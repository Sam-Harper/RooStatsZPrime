//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef DATABOX_hh
#define DATABOX_hh


#include<vector>
#include<string>
#include<map>

#include "ModelConfiguratorZprime.hh"
#include "DataPruner.hh"

#include "TRandom3.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"


using namespace RooFit;
using namespace RooStats;
using namespace std;

class DataBox {

   public: 

      // constructor
      DataBox(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner);
      // destructor
      ~DataBox();

      RooDataSet * createObservedData();
      RooDataSet * createToyMc(); //generates toy data for the expected limit
      RooDataSet * createPeData(); //generates toy data for the evaluation of significance

   private:

      Double_t GetRandom( std::string pdf, std::string var );
      void SetNeventsData(); // store the total number of events in each channel in _neventsDataMap for usage in createPeData(); called in the class constructor

      const std::string legend;
      ModelConfiguratorZprime * _configurator;
      DataPruner * _myDataPruner;
      bool _pruned;
      RooWorkspace * _myWS;
      ModelConfig * _myModelConfig;
      std::map<string, RooDataSet*> _DatamapToy; // contains datasets for the current toy experiment
      std::vector<std::string> _channelnames;
      TRandom3 * _myrandgen;
      std::map<std::string, int> _neventsDataMap; // total number of events in each channel for usage in createPeData()

};
#endif
