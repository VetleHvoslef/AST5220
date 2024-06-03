#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.0; // Since I am a master student

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  // Best fit parameters, h, OmegaK, OmegaM (Omega B konst, endre på OmegaCDM = OmegaM - OmegaB)
  // double h           = 0.67;
  // double OmegaB      = 0.05;
  // double OmegaCDM    = 0.267;
  // double OmegaK      = 0.0;
  // double Neff        = 3.046;
  // double TCMB        = 2.7255;
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  // cosmo.output("cosmology.txt");

  // best_params_fit_cosmology.txt
  // return 0.0;

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  // double kvalue_0 = 0.001 / Constants.Mpc;
  double kvalue_1 = 0.01 / Constants.Mpc;
  // double kvalue_2 = 0.1 / Constants.Mpc;
  // pert.output(kvalue_0, "perturbations_k0.001.txt");
  pert.output(kvalue_1, "perturbations_k0.01.txt");
  // pert.output(kvalue_2, "perturbations_k0.1.txt");
  
  // // Remove when module is completed
  // return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
