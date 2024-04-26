#include "BackgroundCosmology.h"


// Spørsmål: Hvor bør PI defineres fra cmath, i Utils.h under // Physical constants, kanskje mer // Mathematical constants:
// - https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
// - https://stackoverflow.com/questions/6563810/m-pi-works-with-math-h-but-not-with-cmath-in-visual-studio
// - _USE_MATH_DEFINES properly with cmath, men jeg er uvant med å bruke header filer
// spør Kjell og co?

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

  // Spørsmål: Hva er OmegaM i headeren?
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts);


  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  // Også regn ut the age of the Universe her (se side 17, https://cmb.wintherscoming.no/pdfs/ast5220_milestone_1_2024.pdf)

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

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
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
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
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

