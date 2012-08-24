//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef RESULTATOR_hh
#define RESULTATOR_hh


#include<vector>
#include<string>
#include<map>

#include "TFile.h"

#include "ModelConfiguratorZprime.hh"
#include "DataBox.hh"
#include "DataPruner.hh"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCInterval.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

class Resultator {

   public: 

      // constructor
      Resultator(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner, std::string plotfile = "");
      // destructor
      ~Resultator();

      //calculate Bayesian limits using MCMCCalculator for the numerical integration
      void calculateMCMClimit( 
            UInt_t mcmc_iter = 100000,        // number of MCMC iterations
            UInt_t mcmc_burnin = 100,         // number of MCMC burn in steps to be discarded
            UInt_t ntoys = 0,                 // number of pseudoexperiments for expected limit
            std::string suffix = "",         // suffix for output file names
            std::string mode = ""    // obsereved, expected, mass limit (extra k-factor uncertainty)
      );



      void SetPoiUpperLimitByHand(double poiUpperRange){_poiUpperRange = poiUpperRange;}; //FIXME Makeshift function for first tests

      MCMCInterval * GetMcmcInterval(ModelConfig mc, 
               RooDataSet * data,
               double conf_level,
               int n_iter,
               int n_burn,
               double left_side_tail_fraction,
               int n_bins);

      double printMcmcUpperLimit( double peak, std::string filename );
      void printDileptonRatioSignificance( std::string filename, Double_t peak, Double_t pllr);
      
      void setNbinsPosterior(int nbinsPosterior){_nbinsPosterior = nbinsPosterior;};
      int getNbinsPosterior(int nbinsPosterior){return _nbinsPosterior;};

      void makeMcmcPosteriorPlot( MCMCInterval * mcInt, std::string filename);

      Double_t calculateRatioSignificance( std::string mode,
                    Double_t mass_low,
                    Double_t mass_high,
                    Double_t mass_step,
                    Int_t ntoys,
                    std::string suffix,    // suffix for output file names
                    Double_t fit_range_low,
                    Double_t fit_range_high,
                    Double_t width_factor,
                    bool makeLlrPlot // plot option 
                    );
      std::pair<Double_t, Double_t> get_pllr_max( Double_t mass_low,
                       Double_t mass_high,
                       Double_t mass_step,
                       std::string fit_plot_file_name,
                       std::string null_parameters_name,
                       RooDataSet * data );
      Double_t get_pllr( std::string nullParametersName, RooDataSet * data );

      DataBox * getDataBox() { return _databox;};
      double getPoiRangefactor() { return _poirangefactor;};
      void setPoiRangeFactor(double poirangefactor) { _poirangefactor = poirangefactor;};

   protected:

      const std::string legend;
      ModelConfiguratorZprime * _configurator; // pointer to ModelConfigurator object which contains the combined Workspace and ModelConfig
      RooWorkspace * _myWS;  // pointer to the combined workspace
      ModelConfig * _myModelConfig; // pointer to the combined ModelConfig
      DataPruner * _myDataPruner;
      DataBox * _databox; // pointer to DataBox object (returns RooDataSets for the data (observed or PE) )
      double _poiUpperRange; // upper range limit for the POI
      MCMCInterval * _mcInt; 
      bool _bMcmcConverged;
      int _nbinsPosterior; // number of bins for posterior histogram
      bool _makeLlrPlot; // optional plot for significance estimation procedure
      std::string _plotfile;
      double _poirangefactor;
      
};
#endif
