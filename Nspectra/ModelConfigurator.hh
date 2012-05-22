#ifndef MODELCONFIGURATOR_hh
//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#define MODELCONFIGURATOR_hh

#include "Pixie.hh"

#include<vector>
#include<string>
#include<map>

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooCategory.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

class ModelConfigurator {

   public: 
      // constructor
      ModelConfigurator(std::map<string, RooWorkspace*>  WSmap);
      // destructor
      ~ModelConfigurator();

      // creates combined workspace and ModelConfig for the analysis
      virtual void Setup() = 0; 

      RooWorkspace * getCombinedWS(){ return _CombinedWS;}  // get function for the combined WS
      ModelConfig * getCombinedModelConfig() { return _CombinedModelConfig;} // get function for the combined ModelConfig
      void setSharedVars(std::string sharedVarsString);
      const std::string getSharedVars() {return _sharedVarsString;}
      void setPoiString(std::string poiString);
      const std::string getPoiString() {return _poiString;}

      void setNuisanceParams(std::string channelname, std::string nuisanceParamsString);
      void setGlobalObs(std::string channelname, std::string GlobalObsParamsString);
      void setObservables(std::string channelname, std::string ObservableParamsString);

      RooCategory * get_channel_cat(){return channel_cat;};
      std::map<string, RooDataSet*> getChannelDataMap(){return _Datamap;};
      std::map<string, RooDataSet*> * getChannelDataMapPointer (){return &(_Datamap);};
      std::vector<std::string> getChannelNames(){ return _channelnames;};

      void setChannelData(std::string channelname, RooDataSet* channeldata);

      void setVar(std::string name, double value);
      void setVarRange(std::string name, double valuelow, double valuehigh);

      void WriteCombinedWS(std::string filename);
 
   protected:

      std::map<string, RooWorkspace*> _WSmap;
      std::map<string, RooDataSet*> _Datamap;
      const std::string legend;

      RooCategory * channel_cat; // COMMENT : check this

      RooWorkspace * _CombinedWS;
      ModelConfig * _CombinedModelConfig;
      std::string _sharedVarsString;
      std::vector<std::string> _sharedVarsStringVec;
      std::string _poiString; // string for the name of the parameter of interest
      std::vector<std::string> _channelnames;
      std::map<std::string, std::vector<std::string> > _nuisanceStringMap;
      std::map<std::string, std::vector<std::string> > _globalObsStringMap;
      std::map<std::string, std::vector<std::string> > _observablesStringMap;
      Pixie _myPixie;


};
#endif