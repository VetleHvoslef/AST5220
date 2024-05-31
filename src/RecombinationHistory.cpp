#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array;
  Vector Xe_arr;
  Vector ne_arr;

  // Spørsmål: Antar dette her?
  x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Xe_arr = Utils::linspace(x_start, x_end, npts_rec_arrays);
  ne_arr = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Something is very wrong here, is not in Saha regmie for i = 0

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
    } else {
      
      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // Set up IC, solve the ODE and fetch the result 
      //=============================================================================

      Vector x_peebles_array{x_array[i-1], x_array[i]};
      Vector Xe_ic{Xe_arr[i-1]};
      peebles_Xe_ode.solve(dXedx, x_peebles_array, Xe_ic);

      double n_H;
      double n_b;
      double a = exp(x_array[i]);

      // Fetch cosmological parameters
      const double OmegaB      = cosmo->get_OmegaB();
      const double H0          = cosmo->get_H0();

      // Physical constants
      const double G           = Constants.G;
      const double m_H         = Constants.m_H;

      n_b = (1 - Yp) * ((3.0 * pow(H0, 2.0) * OmegaB) / (8.0 * M_PI * G * m_H * pow(a, 3.0)));
      n_H = (1 - Yp) * n_b;

      Xe_arr[i] = peebles_Xe_ode.get_final_data_by_component(0);
      ne_arr[i] = Xe_arr[i] * n_H;
    }
  }

  Vector log_Xe_arr = log(Xe_arr);
  Vector log_ne_arr = log(ne_arr);

  // log_Xe_of_x_spline.create(x_array, log(Xe_arr)); // Hvorfor funker ikke dette?
  log_Xe_of_x_spline.create(x_array, log_Xe_arr);
  log_ne_of_x_spline.create(x_array, log_ne_arr);

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double TCMB        = cosmo->get_TCMB();
  const double H0          = cosmo->get_H0();

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  // Variables for the equation
  double Tb;
  double quad_constant;
  double n_b;
  double n_H;


  Tb = TCMB/a;
  n_b = (1 - Yp) * ((3.0 * pow(H0, 2.0) * OmegaB) / (8.0 * M_PI * G * m_H * pow(a, 3.0)));
  quad_constant = (1/n_b) * pow(((m_e * Tb * k_b) / (2.0 * M_PI * pow(hbar, 2.0))), 3.0/2.0) * exp(-epsilon_0/(Tb * k_b));
  Xe = (2.0 * quad_constant)/(quad_constant + sqrt(quad_constant * quad_constant + 4.0 * quad_constant));
  
  n_H = (1 - Yp) * n_b;
  ne = Xe * n_H;
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double TCMB        = cosmo->get_TCMB();
  const double H           = cosmo->H_of_x(x);
  const double H0          = cosmo->get_H0();


  // Variables for the equation
  double Tb;
  double phi_2_Tb;
  double alpha_2_Tb;
  double beta_Tb;
  double beta_2_Tb;
  double n_b;
  double n_H;
  double n_1s;
  double lambda_alpha;
  double Cr_Tb;


  Tb = TCMB/a;
  phi_2_Tb = 0.488 * log(epsilon_0 / (k_b * Tb));
  alpha_2_Tb = (8.0 / sqrt(3.0 * M_PI)) * c * sigma_T * sqrt(epsilon_0 / (k_b * Tb)) * phi_2_Tb;
  beta_Tb = alpha_2_Tb * pow((m_e * k_b * Tb) / (2.0 * M_PI * pow(hbar, 2.0)), 3.0/2.0) * exp(-(epsilon_0 / (k_b * Tb)));
  beta_2_Tb = alpha_2_Tb * pow((m_e * k_b * Tb) / (2.0 * M_PI * pow(hbar, 2.0)), 3.0/2.0) * exp(-(epsilon_0 / (4.0 * k_b * Tb)));
  n_b = (1 - Yp) * ((3.0 * pow(H0, 2.0) * OmegaB) / (8.0 * M_PI * G * m_H * pow(a, 3.0)));
  n_H = (1 - Yp) * n_b;
  n_1s = (1 - X_e) * n_H;
  lambda_alpha = H * (pow((3.0 * epsilon_0), 3.0) / (pow(8.0 * M_PI, 2.0) * pow(c, 3.0) * pow(hbar, 3.0) * n_1s));
  Cr_Tb = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta_2_Tb);
  
  dXedx[0] = ((Cr_Tb) / H) * (beta_Tb * (1 - X_e) - n_H * alpha_2_Tb * X_e * X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_reverse_array = Utils::linspace(x_end, x_start, npts);
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector tau_array;
  Vector g_tilde_array = Utils::linspace(x_start, x_end, npts);
  ODESolver ode;

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Physical constants in SI units
    const double c           = Constants.c;
    const double sigma_T     = Constants.sigma_T;

    // Cosmological parameters
    const double H           = cosmo->H_of_x(x);

    double ne                = ne_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -((c * ne * sigma_T)/H);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  Vector tau_ic{0.0};
  ode.solve(dtaudx, x_reverse_array, tau_ic);
  tau_array = ode.get_data_by_component(0);
  std::reverse(tau_array.begin(), tau_array.end());

  tau_of_x_spline.create(x_array, tau_array);


  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================

  for (int i=0; i < npts; i++){
    g_tilde_array[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_array);

  // Husk å skrive sound horizon
  // TODO: Fortsett her
  // ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
  //   // Physical constants in SI units
  //   const double c           = Constants.c;
  //   const double sigma_T     = Constants.sigma_T;

  //   // Cosmological parameters
  //   const double H           = cosmo->H_of_x(x);

  //   double ne                = ne_of_x(x);

  //   // Set the derivative for photon optical depth
  //   dsdx[0] = -((c * ne * sigma_T)/H);

  //   return GSL_SUCCESS;
  // };

  // Vector tau_ic{0.0};
  // ode.solve(dtaudx, x_reverse_array, tau_ic);
  // tau_array = ode.get_data_by_component(0);
  // std::reverse(tau_array.begin(), tau_array.end());

  // tau_of_x_spline.create(x_array, tau_array);

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::Tb_of_x(double x) const{
  const double TCMB = cosmo->get_TCMB(x);
  return TCMB/exp(x);
}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

