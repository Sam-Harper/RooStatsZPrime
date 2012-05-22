//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "ModelConfiguratorZprime.hh"

#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "RooAbsData.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

using namespace std;
using namespace RooFit;
//using namespace RooStats;


//Constructor
ModelConfiguratorZprime::ModelConfiguratorZprime( std::map<string, RooWorkspace*> WSmap ) : sublegend("Zprime - "), ModelConfigurator(WSmap) {
   _peak = 1000.;
}


//Destructor
ModelConfiguratorZprime::~ModelConfiguratorZprime(){
;
}

// creates combined workspace and ModelConfig for the analysis
void ModelConfiguratorZprime::Setup() {
   std::string funclegend = "Setup() - ";
   cout << legend << sublegend << funclegend << endl;

   if(_poiString == ""){
      cerr << legend << sublegend << funclegend << "ERROR: poi has not been defined, use ModelConfigurator::setPoiString(std::string poiString)" << endl;
   }

   //iterate over channels
   for(std::map<string, RooWorkspace*>::iterator channelit = _WSmap.begin(); channelit != _WSmap.end() ; channelit++){
      Setup_channel( channelit->first );

   }

   //object for the combined likelihood from all channels
   RooSimultaneous model("model", "model", *channel_cat); 
   //construct combined likelihood from all channels
   for(std::map<string, RooWorkspace*>::iterator channelit = _WSmap.begin(); channelit != _WSmap.end() ; channelit++){
      std::string channel_model_string = "model_";
      channel_model_string += channelit->first; //in Setup_channel the channels's object names have been extended with the channel name
      model.addPdf(*_CombinedWS->pdf(channel_model_string.c_str()), (channelit->first).c_str()); // add  channel's likelihood to the combined model likelihood
   }
   _CombinedWS->import(model, RecycleConflictNodes()); //COMMENT: don't know what RecycleConflictNodes() does

   // construct the prior function for the parameter of interest
   // CONVENTION: poi needs to have the same name in all single channel workspaces
   std::string priorstring = "Uniform::prior(" + _poiString + ")";
   _CombinedWS->factory( priorstring.c_str() );  

   // COMMENT: at this point twobody.C calls Int_t TwoBody::SetDileptonRatioParameters( std::string channel,std::string option ), but all functionality in it is commented out -> ask Gena if this is supposed to be this way (if yes, this comment can be deleted)

   Setup_ModelConfig();

}

void ModelConfiguratorZprime::Setup_channel(std::string channelname) {
   std::string funclegend = "Setup_channel() - ";
   funclegend += channelname;
   cout << legend << sublegend << funclegend << " - begin ..." << endl;

   //add channel name to string vector with channel anmes
   _channelnames.push_back(channelname);

   //handle 
   // COMMENT: check what this is actually doing
   channel_cat->defineType(channelname.c_str());

   //get pointer to the channels's workspace
   RooWorkspace * channelWS = (_WSmap.find(channelname))->second;

   // get pointer to the channel's likelihood function
   // CONVENTION: the channel's likelihood has to named "model" 
   RooAbsPdf * channelmodel = channelWS->pdf("model");

   // check for "shared" in variable name for automatic detection of shared variables
   RooArgSet allVariables = channelWS->allVars();
   cout << "-#-#-#-#-#-#-#-#-" << endl;
   allVariables.Print();
   RooLinkedListIter setIt = allVariables.iterator() ;

   while(setIt.Next()){ 
   //while((*setIt) != 0){ // doesn't wok
      std::string nameString = (*setIt)->GetName();
      cout << legend << funclegend << "checking if variable " << nameString << " is shared" << endl;
      if (nameString.find("shared")!=std::string::npos){
         cout << legend << funclegend << "identified variable " << nameString << " as shared" << endl;
         _sharedVarsString += ",";
         _sharedVarsString += nameString;
      }
    }
   //for ( RooLinkedListIter ){
   //   cout << legend << funclegend << "checking if variable " << (*setIt)->GetName() << " is shared" << endl;
   //}

   // import the channel likelihood into the combined workspace
   // add prefix channelname to all nodes except shared_vars
   // COMMENT: in twobody.C the list of shared variables depends on the number of channels (1 vs. 2), will it cause problems if the list longer or contains undefined variables?
    _CombinedWS->import( *channelmodel,
          RooFit::RenameAllNodes( channelname.c_str() ),
          RooFit::RenameAllVariablesExcept(channelname.c_str(), _sharedVarsString.c_str()) );


   // load the data associated with the channel's workspace
   RooAbsData * channelRooAbsData = channelWS->data("data");
   RooDataSet * channelRooDataSet = new RooDataSet( *(RooDataSet *) channelRooAbsData );
   //insert the data associated with the channel's workspace into map relating it to the name of the channel
   _Datamap.insert( pair<string, RooDataSet *>(channelname, channelRooDataSet) );

   cout << legend << sublegend << funclegend << " - end ..." << endl;
}

void ModelConfiguratorZprime::Setup_ModelConfig() {
   //start from line 1766 in twobody.C 
   std::string funclegend = "Setup_ModelConfig() - ";
   cout << legend << sublegend << funclegend << endl;


   _CombinedModelConfig = new ModelConfig();

   // construct RooArgSet with the parameter of interest
   RooArgSet sPoi( *(_CombinedWS->var(_poiString.c_str())) );
   // construct RooArgSet with the nuisance parameters
   RooArgSet sNuis;
   for( std::map<std::string, std::vector<std::string> >::iterator mapIt = _nuisanceStringMap.begin(); mapIt != _nuisanceStringMap.end(); mapIt++){
      for(std::vector<std::string>::iterator vecIt = ((*mapIt).second).begin(); vecIt != ((*mapIt).second).end(); vecIt++){
         cout << legend << sublegend << funclegend << " add nuisance parameter: " << (*vecIt).c_str() << endl;
         sNuis.add( *(_CombinedWS->var( (*vecIt).c_str() )) );
      }
   }
   // construct RooArgSet with the global observables
   RooArgSet sGlobalObs;
   for( std::map<std::string, std::vector<std::string> >::iterator mapIt = _globalObsStringMap.begin(); mapIt != _globalObsStringMap.end(); mapIt++){
      for(std::vector<std::string>::iterator vecIt = ((*mapIt).second).begin(); vecIt != ((*mapIt).second).end(); vecIt++){
         cout << legend << sublegend << funclegend << " add global observable: " << (*vecIt).c_str() << endl;
         sGlobalObs.add( *(_CombinedWS->var( (*vecIt).c_str() )) );
      }
   }
   // construct RooArgSet with the observables
   //RooArgSet sObs( *(_CombinedWS->var("mass")) );
   RooArgSet sObs;
   for( std::map<std::string, std::vector<std::string> >::iterator mapIt = _observablesStringMap.begin(); mapIt != _observablesStringMap.end(); mapIt++){
      for(std::vector<std::string>::iterator vecIt = ((*mapIt).second).begin(); vecIt != ((*mapIt).second).end(); vecIt++){
         cout << legend << sublegend << funclegend << " add observable: " << (*vecIt).c_str() << endl;
         sObs.add( *(_CombinedWS->var( (*vecIt).c_str() )) );
      }
   }

   //now add the respective variables from all channels ...
   
   // construct the ModelConfig -> associates the elements of the workspace with their statistical meaning 
   _CombinedModelConfig->SetWorkspace(*_CombinedWS);
   _CombinedModelConfig->SetPdf(*(_CombinedWS->pdf("model")));
   _CombinedModelConfig->SetParametersOfInterest( sPoi );
   _CombinedModelConfig->SetPriorPdf( *(_CombinedWS->pdf("prior")) );
   _CombinedModelConfig->SetNuisanceParameters( sNuis );
   _CombinedModelConfig->SetObservables( sObs );
   _CombinedModelConfig->SetGlobalObservables( sGlobalObs );

   cout << legend << sublegend << funclegend << "importing combined ModelConfig" << _CombinedModelConfig->GetName() << " into workspace" << endl;
   
   //load the ModelConfig objects for the combined model into the combined workspace
   _CombinedWS->import(*_CombinedModelConfig);

   cout << legend << sublegend << funclegend << "setup for combined ModelConfig" << _CombinedModelConfig->GetName() << " complete" << endl;

}

void ModelConfiguratorZprime::setMassHypothesis(double peak){
   std::string funclegend = "setMassHypothesis(double peak) - ";
   cout << legend << sublegend << funclegend << "adjusting Z' mass hypothesis from "<< _peak << " to " << peak << endl;
   //CONVENTION: variable for Z' mass must be named "peak"
   _CombinedWS->var("peak")->setVal(peak);
}

double ModelConfiguratorZprime::getMassHypothesis(){
   return _CombinedWS->var("peak")->getVal();
}
