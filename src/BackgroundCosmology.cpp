#include "BackgroundCosmology.h"


// Spørsmål: Hvor bør PI defineres fra cmath, i Utils.h under // Physical constants, kanskje mer // Mathematical constants:
// - https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
// - https://stackoverflow.com/questions/6563810/m-pi-works-with-math-h-but-not-with-cmath-in-visual-studio
// - _USE_MATH_DEFINES properly with cmath, men jeg er uvant med å bruke header filer
// spør Kjell og co?

// Spør Kjell med access denied, git pull i Ubuntu on Windows

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // Compute OmegaR, OmegaNu, OmegaLambda, H0 (alle her er egentlig 0)
  //=============================================================================

  H0 = h * H0_over_h;

  // Check: Sjekk at alle omegene summer opp til 1
  OmegaR = 2 * (pow(M_PI, 2.0) / 30.0) * (pow((k_b * TCMB), 4.0) / (pow(hbar, 3.0) * pow(c, 5.0))) * ((8.0 * M_PI * G) / (3 * pow(H0, 2.0)));
  // Hvor skal dotten være for at den skal skjønne at det ikke er heltallsdivisjon? 8. kontra 8.0
  OmegaNu = Neff * (7.0/8) * pow((4.0/11), (4.0/3)) * OmegaR;
  OmegaK = 0;
  OmegaLambda = 1 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // Set the range of x and the number of points for the splines
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  double H_ini;
  ODESolver ode_eta;
  ODESolver ode_time;


  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // Set the rhs of the detadx ODE
    //=============================================================================

    double Hp = Hp_of_x(x);

    detadx[0] = c / Hp;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================

  H_ini = H_of_x(x_start);
  Vector eta_ic{c / H_ini}; // Spørsmål: Skjønner ikke denne koden her, prøvde eta_ic = c / H_ini med Vector eta_ic sammen med de andre deklarasjonene
  ode_eta.solve(detadx, x_array, eta_ic);
  auto eta_array = ode_eta.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array);

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    //=============================================================================
    // Set the rhs of the detadx ODE
    //=============================================================================

    double H = H_of_x(x);

    dtdx[0] = 1/H;

    return GSL_SUCCESS;
  };

  Vector t_ic{1 / (2.0 * H_of_x(x_start))};
  ode_time.solve(dtdx, x_array, t_ic);
  auto t_array = ode_time.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array);

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  double H;

  H = H0 * sqrt((OmegaR + OmegaNu) * exp(-4.0 * x) + (OmegaB + OmegaCDM) * exp(-3.0 * x) + OmegaK * exp(-2.0 * x) + OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  double H = H_of_x(x);
  double H_p;

  H_p = exp(x) * H;
  return H_p;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double H = H_of_x(x);
  double u;
  double dudx;
  double dHpdx;

  u = (OmegaR + OmegaNu) * exp(-4.0 * x) + (OmegaB + OmegaCDM) * exp(-3.0 * x) + OmegaK * exp(-2.0 * x) + OmegaLambda;
  dudx = -((OmegaR + OmegaNu) * 4.0 * exp(-4.0 * x) + (OmegaB + OmegaCDM) * 3.0 * exp(-3.0 * x) + OmegaK * 2.0 * exp(-2.0 * x));
  dHpdx = exp(x) * (H + H0 * (dudx / (2.0 * sqrt(u))));
  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  double dHpdx = dHpdx_of_x(x);
  double u;
  double dudx;
  double dduddx;
  double ddHpddx;

  u = (OmegaR + OmegaNu) * exp(-4.0 * x) + (OmegaB + OmegaCDM) * exp(-3.0 * x) + OmegaK * exp(-2.0 * x) + OmegaLambda;
  dudx = -((OmegaR + OmegaNu) * 4.0 * exp(-4.0 * x) + (OmegaB + OmegaCDM) * 3.0 * exp(-3.0 * x) + OmegaK * 2.0 * exp(-2.0 * x));
  dduddx = (OmegaR + OmegaNu) * 16.0 * exp(-4.0 * x) + (OmegaB + OmegaCDM) * 9.0 * exp(-3.0 * x) + OmegaK * 4.0 * exp(-2.0 * x);
  ddHpddx = dHpdx + exp(x) * (H0 * (dudx / (2.0 * sqrt(u))) + H0 * ((2 * dduddx * sqrt(u) - dudx * dudx) / (4.0 * u * sqrt(u))));
  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB; // Spørsmål: Forandre til OmegaB0 for alle omega-ene, etter fristen

  double OmegaB_x;
  double H = H_of_x(x);
  
  OmegaB_x = OmegaB / (exp(3.0 * x) * (pow(H, 2.0) / pow(H0, 2.0)));
  return OmegaB_x;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double OmegaR_x;
  double H = H_of_x(x);

  OmegaR_x = OmegaR / (exp(4.0 * x) * (pow(H, 2.0) / pow(H0, 2.0)));
  return OmegaR_x;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double OmegaNu_x;
  double H = H_of_x(x);

  OmegaNu_x = OmegaNu / (exp(4.0 * x) * (pow(H, 2.0) / pow(H0, 2.0)));
  return OmegaNu_x;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double OmegaCDM_x;
  double H = H_of_x(x);

  OmegaCDM_x = OmegaCDM / (exp(3.0 * x) * (pow(H, 2.0) / pow(H0, 2.0)));
  return OmegaCDM_x;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double OmegaLambda_x;
  double H = H_of_x(x);

  OmegaLambda_x = OmegaLambda / (pow(H, 2.0) / pow(H0, 2.0));
  return OmegaLambda_x;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  double OmegaK_x;
  double H = H_of_x(x);
  
  OmegaK_x = OmegaK / (exp(2.0 * x) * (pow(H, 2.0) / pow(H0, 2.0)));
  return OmegaK_x;
}

double BackgroundCosmology::get_OmegaM(double x) const{
  if (x == 0.0) return OmegaB + OmegaCDM;

  return get_OmegaB(x) + get_OmegaCDM(x);
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double chi;
  double eta_0;

  eta_0 = eta_of_x_spline(x_end);
  chi = eta_0 - eta_of_x(x);
  return chi;
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  double r;
  double d_A;
  double arg;
  double OmegaK0 = get_OmegaK();
  double chi = get_comoving_distance_of_x(x);
  double H0 = get_H0();

  arg = sqrt(abs(OmegaK0)) * H0 * chi/c;

  if (OmegaK0 < 0){
    r = chi * sin(arg) / arg;
  }
  
  if (OmegaK0 == 0){
    r = chi;
  }
  
  if (OmegaK0 > 0){
    r = chi * sinh(arg) / arg;
  }

  d_A = exp(x) * r;
  return d_A;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double d_L;
  double d_A = get_angular_distance_of_x(x);

  d_L = d_A / exp(2.0 * x);
  return d_L;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::time_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << time_of_x(x)       << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

