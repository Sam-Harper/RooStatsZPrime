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

      // constructor
      DataPruner(std::map<std::string , double> Rangemap);
      // destructor
      ~DataPruner();

      void Prune(std::map<string, RooDataSet*> * _Datamap);

   private:

      void Prune(std::string channelname, double massThreshold);

      const std::string legend;
      //ModelConfiguratorZprime * _configurator; // pointer to ModelConfigurator object which contains the combined Workspace and ModelConfig
      std::vector<std::string> _channelnames;
      std::map<std::string , double> _Rangemap; // stores limits to masses in the datasets of the channels
      std::map<string, RooDataSet*> * _Datamap;

};
#endif