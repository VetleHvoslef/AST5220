#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array;
  // Vector log_k_array = log(k_array);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  k_array = exp(log_k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  int n_z = 30000;
  Vector z_array = Utils::linspace(0, 30000, n_z);
  double z;
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    Vector j_ell_array = Utils::linspace(0, 30000, n_z);
    for (int iz=0; iz < n_z; iz++){
      z = z_array[iz];
      j_ell_array[iz] = Utils::j_ell(ell, z);
    }
    j_ell_splines[i].create(z_array, j_ell_array);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(// Dette kan være feil
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  Vector x_array  = Utils::linspace(x_start, x_end, n_x); 

  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    double k = k_array[ik];
    double eta_0 = cosmo->eta_of_x(x_end);
    double delta_x_ix;
    double x_minus_1;
    double x;
    double f_x_minus_1;
    double f_x;

    for (int iell=0; iell < ells.size(); iell++){
      // Integrate over x (trapizodial rule)
      double F_ell_k = 0.0;

      for (int ix=1; ix < x_array.size(); ix++){
        x = x_array[ix];
        x_minus_1 = x_array[ix - 1];
        delta_x_ix = x - x_minus_1;

        f_x_minus_1 = source_function(x_minus_1, k) * j_ell_splines[iell](k * (eta_0 - cosmo->eta_of_x(x_minus_1)));
        f_x = source_function(x, k) * j_ell_splines[iell](k * (eta_0 - cosmo->eta_of_x(x)));
      
        F_ell_k = F_ell_k + ((f_x_minus_1 + f_x) / 2.0) * delta_x_ix;
      }

      // Store the result for Source_ell(k) in result[ell][ik] (ble litt forvirret)
      // skal vel store theta_l(k)?
      result[iell][ik] = F_ell_k;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  // k_array, men ikke ell array?

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int iell=0; iell < nells; iell++){
    thetaT_ell_of_k_spline[iell].create(k_array, thetaT_ell_of_k[iell]);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  // if(Constants.polarization){

  //   // ...
  //   // ...
  //   // ...
  //   // ...

  // }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell( // C_ell til neste gang :)
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate C_ell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  Vector result = Vector(ells.size());

  for (int iell=0; iell < nells; iell++){
    double C_ell_k = 0.0;

    // Integrate over k (trapizodial rule)
    for (int ik=1; ik < log_k_array.size(); ik++){
      double delta_k_ik;
      double k_minus_1;
      double k;
      double f_k_minus_1;
      double f_k;
      k = exp(log_k_array[ik]);
      k_minus_1 = exp(log_k_array[ik - 1]);
      delta_k_ik = log_k_array[ik] - log_k_array[ik - 1];

      f_k_minus_1 = primordial_power_spectrum(k_minus_1) * f_ell_spline[iell](k_minus_1) * g_ell_spline[iell](k_minus_1);
      f_k = primordial_power_spectrum(k) * f_ell_spline[iell](k) * g_ell_spline[iell](k);
      
      C_ell_k = C_ell_k + ((f_k_minus_1 + f_k) / 2.0) * delta_k_ik;
    }

    result[iell] = 4.0 * M_PI * C_ell_k;
  }
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;
  double c = Constants.c;
  double OmegaM0 = cosmo->get_OmegaM();
  double k = k_mpc / Constants.Mpc;
  double Phi = pert->get_Phi(x, k);
  double H0 = cosmo->get_H0();
  double a = exp(x);
  double delta_M;

  //=============================================================================
  // Compute the matter power spectrum
  //=============================================================================
  delta_M = (c * c * k * k * Phi) / ((3.0 / 2.0) * OmegaM0 * pow(a, -1.0) * H0 * H0);

  pofk = pow(std::abs(delta_M), 2.0) * primordial_power_spectrum(k) * 2.0 * M_PI * M_PI / pow(k, 3.0);
  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2); // Dette her blir ikke brukt?
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI); // Dette her blir ikke brukt?
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    // if(Constants.polarization){
    //   fp << cell_EE_spline( ell ) * normfactor  << " ";
    //   fp << cell_TE_spline( ell ) * normfactor  << " ";
    // }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);

  std::ofstream fp_k("functions_of_k.txt");
  auto log_k_values = Utils::linspace(log(k_min), log(k_max), n_k);
  auto kvalues = exp(log_k_values);
  auto print_k_data = [&] (const double k) {
    fp_k << k                                                   << " ";
    fp_k << thetaT_ell_of_k_spline[0](k)                        << " "; // Indeks 0, ell 2
    fp_k << thetaT_ell_of_k_spline[10](k)                       << " "; // Indeks 10, ell 20
    fp_k << thetaT_ell_of_k_spline[20](k)                       << " "; // Indeks 20, ell 120
    fp_k << get_matter_power_spectrum(0.0, k * Constants.Mpc)   << " ";
    fp_k << "\n";
  };
  std::for_each(kvalues.begin(), kvalues.end(), print_k_data);
}
