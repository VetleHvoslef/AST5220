#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

// TENK LITT PÅ DETTE: Testing av kode
/*
using gdb to find all the variables that are created at an iteration in a "for loop"
kan sette en verdi i starten av for løkken som jeg skal teste mot
file:///C:/Users/vetle/Downloads/The.Art.of.Debugging.with.GDB.DDD.and.Eclipse.pdf (åpnes i nettleseren)
- gui gdb:
  - https://www.gdbgui.com/
  - http://cgdb.github.io/
*/

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k);
  Vector k_array_lin_log = Utils::linspace(log(k_min), log(k_max), n_k);
  k_array = exp(k_array_lin_log);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Arrays for splining
  Vector2D delta_cdm_array(n_x, Vector(n_k, 0.0));
  Vector2D delta_b_array(n_x, Vector(n_k, 0.0));
  Vector2D v_cdm_array(n_x, Vector(n_k, 0.0));
  Vector2D v_b_array(n_x, Vector(n_k, 0.0));
  Vector2D Phi_array(n_x, Vector(n_k, 0.0));
  Vector2D PI_array(n_x, Vector(n_k, 0.0));
  Vector2D Psi_array(n_x, Vector(n_k, 0.0));

  Vector2D Theta_0_array(n_x, Vector(n_k, 0.0));
  Vector2D Theta_1_array(n_x, Vector(n_k, 0.0));
  Vector2D Theta_2_array(n_x, Vector(n_k, 0.0));

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    auto x_end_tight_data = get_tight_coupling_time(k);

    double x_end_tight = x_end_tight_data.first;
    int x_end_tight_ind = x_end_tight_data.second;

    int n_x_tight = x_end_tight_ind + 1;
    int n_x_full = (n_x - n_x_tight) + 1;

    Vector x_array_tight_system = Utils::linspace(x_start, x_end_tight, n_x_tight);
    Vector x_array_full_system = Utils::linspace(x_end_tight, x_end, n_x_full);

    //===================================================================
    // Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODESolver tight_coupling_ode;
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    tight_coupling_ode.solve(dydx_tight_coupling, x_array_tight_system, y_tight_coupling_ini);
    auto y_tc = tight_coupling_ode.get_data();

    //===================================================================
    // Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_tight_coupling = tight_coupling_ode.get_final_data();
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODESolver full_system_ode;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end;
    full_system_ode.solve(dydx_full, x_array_full_system, y_full_ini);
    auto y = full_system_ode.get_data();

    for (int ix=0; ix < n_x; ix++){
      // Constants and such
      double c                = Constants.c;
      double x                = x_array[ix];
      double a                = exp(x);
      double Hp               = cosmo->Hp_of_x(x);
      double H0               = cosmo->get_H0();
      double OmegaR0          = cosmo->get_OmegaR();
      double OmegaNu0         = cosmo->get_OmegaNu();

      double dtaudx           = rec->dtaudx_of_x(x);

      double Nu_2             = 0;
      double Theta_p_0        = 0;
      double Theta_p_2        = 0;

      // Tight coupling
      if (ix <= x_end_tight_ind){
        delta_cdm_array[ix][ik] = y_tc[ix][Constants.ind_deltacdm_tc];
        delta_b_array[ix][ik]   = y_tc[ix][Constants.ind_deltab_tc];
        v_cdm_array[ix][ik]     = y_tc[ix][Constants.ind_vcdm_tc];
        v_b_array[ix][ik]       = y_tc[ix][Constants.ind_vb_tc];
        Phi_array[ix][ik]       = y_tc[ix][Constants.ind_Phi_tc];

        Theta_0_array[ix][ik]   = y_tc[ix][Constants.ind_start_theta_tc];
        Theta_1_array[ix][ik]   = y_tc[ix][Constants.ind_start_theta_tc + 1];
        Theta_2_array[ix][ik]   = -((20.0 * c * k) / (45.0 * Hp * dtaudx)) * Theta_1_array[ix][ik];
        
        Psi_array[ix][ik]       = -Phi_array[ix][ik] - ((12.0 * H0 * H0) / (pow(c * k * a, 2.0))) * (OmegaR0 * Theta_2_array[ix][ik] + OmegaNu0 * Nu_2);
        PI_array[ix][ik]        = Theta_2_array[ix][ik] + Theta_p_0 + Theta_p_2;
      }
      else{
        delta_cdm_array[ix][ik] = y[ix - x_end_tight_ind][Constants.ind_deltacdm];
        delta_b_array[ix][ik]   = y[ix - x_end_tight_ind][Constants.ind_deltab];
        v_cdm_array[ix][ik]     = y[ix - x_end_tight_ind][Constants.ind_vcdm];
        v_b_array[ix][ik]       = y[ix - x_end_tight_ind][Constants.ind_vb];
        Phi_array[ix][ik]       = y[ix - x_end_tight_ind][Constants.ind_Phi];

        Theta_0_array[ix][ik]   = y[ix - x_end_tight_ind][Constants.ind_start_theta];
        Theta_1_array[ix][ik]   = y[ix - x_end_tight_ind][Constants.ind_start_theta + 1];
        Theta_2_array[ix][ik]   = y[ix - x_end_tight_ind][Constants.ind_start_theta + 2];

        Psi_array[ix][ik]       = -Phi_array[ix][ik] - ((12.0 * H0 * H0) / (pow(c * k * a, 2.0))) * (OmegaR0 * Theta_2_array[ix][ik] + OmegaNu0 * Nu_2);
        PI_array[ix][ik]        = Theta_2_array[ix][ik] + Theta_p_0 + Theta_p_2;
      }
    }

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================

  delta_cdm_spline.create(x_array, k_array, delta_cdm_array);
  delta_b_spline.create(x_array, k_array, delta_b_array);
  v_cdm_spline.create(x_array, k_array, v_cdm_array);
  v_b_spline.create(x_array, k_array, v_b_array);
  Phi_spline.create(x_array, k_array, Phi_array);
  Psi_spline.create(x_array, k_array, Psi_array);

  Theta_0_spline.create(x_array, k_array, Theta_0_array);
  Theta_1_spline.create(x_array, k_array, Theta_1_array);
  Theta_2_spline.create(x_array, k_array, Theta_2_array);
  Theta_spline = {Theta_0_spline, Theta_1_spline, Theta_2_spline};

  PI_spline.create(x_array, k_array, PI_array);
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Nesecary constants and variables: f_v = 0 (no neutrinos)
  const double Psi         = -2.0/3.0;
  const double c           = Constants.c;  
  const double Hp          = cosmo->Hp_of_x(x);
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================


  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  Phi = -Psi;
  delta_cdm = -(3.0/2.0) * Psi;
  delta_b = delta_cdm;
  v_cdm = -((c * k) / (2.0 * Hp)) * Psi;
  v_b = v_cdm;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -(1.0/2.0) * Psi;
  Theta[1] = ((c * k) / (6.0 * Hp)) * Psi;


  // // SET: Neutrino perturbations (N_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  // }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  // Nesecary constants and variables
  const double c          = Constants.c;

  double Hp               = cosmo->Hp_of_x(x);
  double dtaudx           = rec->dtaudx_of_x(x);

  //=============================================================================
  // Fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = -((20.0 * c * k) / (45.0 * Hp * dtaudx)) * Theta[1];

  for (int il=3; il < n_ell_theta; il++){
    Theta[il] = -(il / (2.0 * il + 1)) * ((c * k) / (Hp * dtaudx)) * Theta[il - 1];
  }

  // TODO: Hvis setter Constants.polarization og Constants.neutrinos til False i Utils header filen så 
  // kræsjer koden, så setter dem til True og kommenterer ut alle if testene

  // // SET: Photon polarization perturbations (Theta_p_ell)
  // if(polarization){
  //   // ...
  //   // ...
  // }

  // // SET: Neutrino perturbations (N_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  // }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double,int>  Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;
  int x_tight_coupling_end_ind = 0;
  double c = Constants.c;
  double dtaudx;
  double min_var;
  double Hp;
  double x;

  //=============================================================================
  // Compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  for (int ix=0; ix < n_x; ix++){
    x = x_array[ix];
    Hp = cosmo->Hp_of_x(x);

    min_var = std::min(1.0, (c * k) / Hp);
    dtaudx = rec->tau_of_x(x);
    x_tight_coupling_end = x;
    x_tight_coupling_end_ind = ix;

    if (std::abs(dtaudx) <= 10 * min_var || x >= -8.3){
      break;
    }
  }

  return std::pair<double,int>(x_tight_coupling_end, x_tight_coupling_end_ind);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array;
  Vector k_array_lin_log = Utils::linspace(log(k_min), log(k_max), n_k);
  k_array = exp(k_array_lin_log);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      const double Hp        = cosmo->Hp_of_x(x);
      const double tau       = rec->tau_of_x(x);
      const double g_tilde   = rec->g_tilde_of_x(x);
      const double Theta_0    = get_Theta(x, k, 0);
      const double Psi       = get_Psi(x, k);
      const double v_b       = get_v_b(x, k);
      const double PI_eq     = get_PI(x, k);
      const double c         = Constants.c;

      const double dPsidx    = Psi_spline.deriv_x(x, k);
      const double dPhidx    = Phi_spline.deriv_x(x, k);
      const double dgdx_tilde = rec->dgdx_tilde_of_x(x);
      const double ddgddx_tilde = rec->ddgddx_tilde_of_x(x);
      const double dHpdx = cosmo->dHpdx_of_x(x);
      const double dv_bdx = v_b_spline.deriv_x(x, k);
      const double ddHpddx = cosmo->ddHpddx_of_x(x);
      const double dPI_eqdx = PI_spline.deriv_x(x, k);
      const double ddPI_eqddx = PI_spline.deriv_xx(x, k);

      const double dHpg_tilde_v_bdx = dHpdx * g_tilde * v_b + Hp * (dgdx_tilde * v_b + g_tilde * dv_bdx);
      const double dHpdHpg_tildePIddx = (dHpdx * dHpdx + Hp * ddHpddx) * g_tilde * PI_eq +
                                        3.0 * Hp * dHpdx * (dgdx_tilde * PI_eq + g_tilde * dPI_eqdx) +
                                        Hp * Hp * (ddgddx_tilde * PI_eq + 2.0 * dgdx_tilde * dPI_eqdx + g_tilde * ddPI_eqddx);

      // Temperatur source
      ST_array[index] = g_tilde * (Theta_0 + Psi + (1.0/4.0) * PI_eq) + 
                        exp(-tau) * (dPsidx - dPhidx) -
                        (1.0 / (c * k)) * dHpg_tilde_v_bdx + 
                        (3.0 / (4.0 * c * c * k * k)) * dHpdHpg_tildePIddx;

      // // Polarization source
      // if(Constants.polarization){
      //   SE_array[index] = 0.0;
      // }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  // if(Constants.polarization){
  //   SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
  // }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  // Nesecary constants and variables
  const double OmegaR0     = cosmo->get_OmegaR();
  const double OmegaB0     = cosmo->get_OmegaB();
  const double OmegaNu0    = cosmo->get_OmegaNu();
  const double OmegaCDM0   = cosmo->get_OmegaCDM();
  const double H0          = cosmo->get_H0();
  const double c           = Constants.c;

  double a           = exp(x);
  double Hp          = cosmo->Hp_of_x(x);
  double dHpdx       = cosmo->dHpdx_of_x(x);
  double dtaudx      = rec->dtaudx_of_x(x);
  double ddtauddx    = rec->ddtauddx_of_x(x);
  double Nu_0         = 0;
  double Nu_2         = 0;
  double R;
  double Theta2;
  double top;
  double q;
  double Psi;
  double Omega_const;

  //=============================================================================
  // Fill in the expressions for all the derivatives
  //=============================================================================

  // Calculated variables
  Theta2 = -((20.0 * c * k) / (45.0 * Hp * dtaudx)) * Theta[1];
  R = (4.0 * OmegaR0) / (3.0 * OmegaB0 * a);

  // SET: Scalar quantities (Phi, delta, v, ...)
  Psi = -Phi - ((12.0 * H0 * H0) / (pow(c * k * a, 2.0))) * (OmegaR0 * Theta2 + OmegaNu0 * Nu_2);
  Omega_const = OmegaCDM0 * pow(a, -1.0) * delta_cdm + OmegaB0 * pow(a, -1.0) * delta_b + 4.0 * OmegaR0 * pow(a, -2.0) * Theta[0] + 4.0 * OmegaNu0 * pow(a, -2.0) * Nu_0;
  dPhidx = Psi - ((pow(c * k, 2.0)) / (3.0 * Hp * Hp)) * Phi + ((H0 * H0) / (2.0 * Hp * Hp)) * Omega_const;
  ddelta_cdmdx = ((c * k)/Hp) * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ((c * k)/Hp) * Psi;
  ddelta_bdx = ((c * k)/Hp) * v_b - 3.0 * dPhidx;

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -((c * k)/Hp) * Theta[1] - dPhidx;

  top = -((1 - R) * dtaudx + (1 + R) * ddtauddx) * (3.0 * Theta[1] + v_b) - ((c * k)/Hp) * Psi + (1 - (dHpdx/Hp)) * ((c * k)/Hp) * (-Theta[0] + 2.0 * Theta2) - ((c * k)/Hp) * dThetadx[0];
  q = top / ((1 + R) * dtaudx + dHpdx/Hp - 1);

  dv_bdx = (1 / (1 + R)) * (-v_b - ((c * k)/Hp) * Psi + R * (q + ((c * k)/Hp) * (-Theta[0] + 2.0 * Theta2) - ((c * k)/Hp) * Psi));
  dThetadx[1] = 1.0/3.0 * (q - dv_bdx);
  // // SET: Neutrino mutlipoles (Nu_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  //   // ...
  // }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Nesseraies constants
  const double c          = Constants.c;

  // Cosmological parameters and variables
  double a                = exp(x);
  double Hp               = cosmo->Hp_of_x(x);
  double eta              = cosmo->eta_of_x(x);
  double H0               = cosmo->get_H0();
  double OmegaR0          = cosmo->get_OmegaR();
  double OmegaNu0         = cosmo->get_OmegaNu();
  double OmegaCDM0        = cosmo->get_OmegaCDM();
  double OmegaB0          = cosmo->get_OmegaB();

  // Recombination variables
  double dtaudx           = rec->dtaudx_of_x(x);

  double Nu_0             = 0;
  double Nu_2             = 0;
  double Theta_p_0        = 0;
  double Theta_p_2        = 0;
  double PI_eq;
  double R;
  double Omega_const;
  double Psi;

  // if(k*Constants.Mpc >  0.124592-0.0002 and k*Constants.Mpc <  0.124592 + 0.0002 and rand() % 100 == 0){
  //   std::cout << "x = " << x << " " << Phi << " " << delta_b << std::endl;
  //   if(std::abs(Phi) > 2.0) exit(1);
  // }


  //=============================================================================
  // Fill in the expressions for all the derivatives
  //=============================================================================

  R = (4.0 * OmegaR0) / (3.0 * OmegaB0 * a);

  // SET: Scalar quantities (Phi, delta, v, ...)
  Psi = -Phi - ((12.0 * H0 * H0) / (pow(c * k * a, 2.0))) * (OmegaR0 * Theta[2] + OmegaNu0 * Nu_2);
  Omega_const = OmegaCDM0 * pow(a, -1.0) * delta_cdm + OmegaB0 * pow(a, -1.0) * delta_b + 4.0 * OmegaR0 * pow(a, -2.0) * Theta[0] + 4.0 * OmegaNu0 * pow(a, -2.0) * Nu_0;
  dPhidx = Psi - ((pow(c * k, 2.0)) / (3.0 * Hp * Hp)) * Phi + ((H0 * H0) / (2.0 * Hp * Hp)) * Omega_const;
  ddelta_cdmdx = ((c * k)/Hp) * v_cdm - 3.0 * dPhidx;
  dv_cdmdx = -v_cdm - ((c * k)/Hp) * Psi;
  ddelta_bdx = ((c * k)/Hp) * v_b - 3.0 * dPhidx;
  dv_bdx = -v_b - ((c * k)/Hp) * Psi + dtaudx * R * (3.0 * Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -((c * k)/Hp) * Theta[1] - dPhidx;
  dThetadx[1] = ((c * k) / (3.0 * Hp)) * Theta[0] - ((2.0 * c * k) / (3.0 * Hp)) * Theta[2] + ((c * k) / (3.0 * Hp)) * Psi + dtaudx * (Theta[1] + (1.0/3.0) * v_b);

  for (int il=2; il < n_ell_theta - 1; il++){
    if (il == 2){
      PI_eq = Theta[2] + Theta_p_0 + Theta_p_2;
    } 
    else{
      PI_eq = 0;
    }
    dThetadx[il] = ((il * c * k) / ((2.0 * il + 1) * Hp)) * Theta[il - 1] - (((il + 1) * c * k) / ((2.0 * il + 1) * Hp)) * Theta[il + 1] + dtaudx * (Theta[il] - (1.0/10.0) * PI_eq);
  }

  dThetadx[n_ell_theta - 1] = ((c * k)/Hp) * Theta[n_ell_theta - 2] - c * ((n_ell_theta + 1) / (Hp * eta)) * Theta[n_ell_theta - 1] + dtaudx * Theta[n_ell_theta - 1];

  // // SET: Photon polarization multipoles (Theta_p_ell)
  // if(polarization){
  //   // ...
  //   // ...
  //   // ...
  // }

  // // SET: Neutrino mutlipoles (Nu_ell)
  // if(neutrinos){
  //   // ...
  //   // ...
  //   // ...
  // }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_PI(const double x, const double k) const{
  return PI_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
// double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
//   return Theta_p_spline[ell](x,k);
// }
// double Perturbations::get_Nu(const double x, const double k, const int ell) const{
//   return Nu_spline[ell](x,k);
// }

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_PI(x,k)        << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << " | ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

