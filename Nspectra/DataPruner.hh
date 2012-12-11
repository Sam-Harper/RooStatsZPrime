//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef DATAPRUNER_hh
#define DATAPRUNER_hh


#include<vector>
#include<string>
#include<map>

#include "ModelConfiguratorZprime.hh"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

class DataPruner {

   public: 

      // constructor: adjust lower range for set of given channels
      DataPruner(std::map<std::string , double> Rangemap, ModelConfiguratorZprime * configurator);
      // constructor: adjust range of all channels according to signal hypothesis
      DataPruner(double n_width, ModelConfiguratorZprime * configurator);
      // destructor
      ~DataPruner();

      void Prune(std::map<string, RooDataSet*> * _Datamap);

   private:

      void Prune(std::string channelname, double massThreshold);
      void Prune(std::string channelname, double massThreshLow, double massThreshHigh);

      
      const std::string legend;
      std::string option;
      double _n_width;
      ModelConfiguratorZprime * _configurator; // pointer to ModelConfigurator object which contains the combined Workspace and ModelConfig
      std::vector<std::string> _channelnames;
      std::map<std::string , double> _Rangemap; // stores limits to masses in the datasets of the channels
      std::map<string, RooDataSet*> * _Datamap;

};
#endif
