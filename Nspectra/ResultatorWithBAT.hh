//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef ResultatorWithBATWITHBAT_hh
#define ResultatorWithBATWITHBAT_hh


#include<vector>
#include<string>
#include<map>

#include "ModelConfiguratorZprime.hh"
#include "DataBox.hh"
#include "DataPruner.hh"
#include "Resultator.hh"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCInterval.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

class ResultatorWithBAT : public Resultator {

   public: 

      // constructor
      ResultatorWithBAT(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner);
      // destructor
      ~ResultatorWithBAT();

      //calculate Bayesian limits using MCMCCalculator for the numerical integration
      void calculateMCMClimit( 
            UInt_t mcmc_iter = 100000,        // number of MCMC iterations
            UInt_t mcmc_burnin = 100,         // number of MCMC burn in steps to be discarded
            UInt_t ntoys = 0,                 // number of pseudoexperiments for expected limit
            std::string suffix = "",         // suffix for output file names
            std::string mode = ""    // obsereved, expected, mass limit (extra k-factor uncertainty)
      );

      MCMCInterval * GetMcmcInterval(ModelConfig mc, 
               RooDataSet * data,
               double conf_level,
               int n_iter,
               int n_burn,
               double left_side_tail_fraction,
               int n_bins);

      double printMcmcUpperLimit( double peak, std::string filename );

   private:

      std::string sublegend;

};
#endif