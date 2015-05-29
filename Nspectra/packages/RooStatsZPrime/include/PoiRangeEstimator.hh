//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef POIRANGEESTIMATOR_hh
#define POIRANGEESTIMATOR_hh


#include<vector>
#include<string>
#include<map>

#include "ModelConfiguratorZprime.hh"
#include "DataBox.hh"
#include "Resultator.hh"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

class PoiRangeEstimator {

   public: 

      // constructor
      PoiRangeEstimator(ModelConfiguratorZprime* configurator, Resultator * myResultator);
      PoiRangeEstimator(ModelConfiguratorZprime* configurator, RooDataSet * inputdata);
      // destructor
      ~PoiRangeEstimator();

      Double_t GetPoiUpper(std::string channel, Double_t peak);
      double GetPoiUpperSimple(std::string channel, Double_t peak);
      LikelihoodInterval * GetPlrInterval( double conf_level );

   private:

      const std::string legend;
      ModelConfiguratorZprime * _configurator; // pointer to ModelConfigurator object which contains the combined Workspace and ModelConfig
      RooWorkspace * _myWS;  // pointer to the combined workspace
      ModelConfig * _myModelConfigCopy; // pointer to the combined ModelConfig
      Resultator * _myResultator;

      LikelihoodInterval * mpPlrInt;
      RooDataSet * data;

};
#endif
