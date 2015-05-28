//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "ModelConfigurator.hh"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include "RooRealVar.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;


//--------------------Constructor
ModelConfigurator::ModelConfigurator(std::map<string, RooWorkspace*> WSmap) : legend("ModelConfigurator - "), _myPixie(),verbose_(false) {

   _WSmap = WSmap;
   _CombinedWS = new RooWorkspace("CombinedWS");
   _sharedVarsString = "";
   _poiString = "";
   channel_cat = new RooCategory("channels", "channels"); //COMMENT: check why/ what for we need this
   if(verbose_) cout << legend << "constructed" << endl;
}

ModelConfigurator::~ModelConfigurator() {
   //FIXME: delete the map entries of _WSmap ?!
   //FIXME: delete the map entries of _Datamap ?!
   delete _CombinedModelConfig;
   delete _CombinedWS;
   delete channel_cat;
   if(verbose_) cout << legend << "deleted" << endl;
}

void ModelConfigurator::setSharedVars(std::string sharedVarsString) {
   if (_sharedVarsString == ""){
      _sharedVarsString = sharedVarsString;
      _sharedVarsStringVec = _myPixie.commaSepStringToStringVec(sharedVarsString);
   }
   else{
      if(verbose_) cout << legend << "cannot change shared variables once they are defined" << endl;
      if(verbose_) cout << legend << "shared variables are: " << _sharedVarsString << endl;
   }
}

void ModelConfigurator::setPoiString(std::string poiString) {
   if (_poiString == ""){
      _poiString = poiString;
      //FIXME: issue warning if poiString does not correspond to a variable in the workspaces
   }
   else{
      if(verbose_) cout << legend << "cannot change poi once it is defined" << endl;
      if(verbose_) cout << legend << "poi is: " << _poiString << endl;
   }
}

void ModelConfigurator::setNuisanceParams(std::string channelname, std::string nuisanceParamsString){

   std::string funclegend = "setNuisanceParams(std::string channelname, std::string nuisanceParamsString) - ";
   funclegend += channelname;

   if(_nuisanceStringMap.find(channelname) == _nuisanceStringMap.end()){
      std::vector<std::string> nuisanceStringVec = _myPixie.commaSepStringToStringForRooArgSet(channelname, nuisanceParamsString, _sharedVarsStringVec);
      _nuisanceStringMap.insert( pair<string, std::vector<std::string> >(channelname, nuisanceStringVec) );
   }
   else{
      if(verbose_) cout << legend << funclegend << " cannot change these parameters once they have been defined " << endl;
   }
}

void ModelConfigurator::setGlobalObs(std::string channelname, std::string GlobalObsParamsString){

   std::string funclegend = "setGlobalObs(std::string channelname, std::string GlobalObsParamsString) - ";
   funclegend += channelname;

   if(_globalObsStringMap.find(channelname) == _globalObsStringMap.end()){
      std::vector<std::string> globalObsStringVec = _myPixie.commaSepStringToStringForRooArgSet(channelname, GlobalObsParamsString, _sharedVarsStringVec);
      _globalObsStringMap.insert( pair<string, std::vector<std::string> >(channelname, globalObsStringVec) );
   }
   else{
      if(verbose_) cout << legend << funclegend << " cannot change these parameters once they have been defined " << endl;
   }
}

void ModelConfigurator::setObservables(std::string channelname, std::string ObservableParamsString){

   std::string funclegend = "setObservables(std::string channelname, std::string ObservableParamsString) - ";
   funclegend += channelname;

   if(_observablesStringMap.find(channelname) == _observablesStringMap.end()){
      std::vector<std::string> observablesStringVec = _myPixie.commaSepStringToStringForRooArgSet(channelname, ObservableParamsString, _sharedVarsStringVec);
      _observablesStringMap.insert( pair<string, std::vector<std::string> >(channelname, observablesStringVec) );
   }
   else{
      if(verbose_) cout << legend << funclegend << " cannot change these parameters once they have been defined " << endl;
   }
}

void ModelConfigurator::WriteCombinedWS(std::string filename){

   std::string funclegend = "WriteCombinedWS(std::string filename) - ";

   TFile* file = new TFile( filename.c_str() ,"RECREATE"); 
   if ( file->IsOpen() ) if(verbose_) cout << legend << funclegend << "File "<< filename <<" opened successfully" << endl;  
   file->cd();

   //make sure that user defined RooFit functions are saved
   _CombinedWS->importClassCode();

   file->WriteTObject(_CombinedWS);
   file->Close();
   if ( file->IsOpen() ) if(verbose_) cout << legend << funclegend << "combined workspace has been written to file " << filename << endl;

}

void ModelConfigurator::setChannelData(std::string channelname, RooDataSet* channeldata){

   std::string funclegend = "setChannelData(std::string channelname, RooDataSet* _Datamap) - ";

   if ( find(_channelnames.begin(), _channelnames.end(), channelname) != _channelnames.end() ){
      if(verbose_) cout << legend << funclegend << "WARNING: replacing dataset for channel " << channelname << endl;
      if(verbose_) cout << legend << funclegend << "WARNING: dataset for channel will be used for calculations but will not be saved to the combined WorkSpace"<< endl; 
      //FIXME: not saving to the WS might not be the desired solution?
      RooDataSet * temppointer = _Datamap.find(channelname)->second;
      _Datamap.find(channelname)->second = channeldata;
      delete temppointer;
   }
   else{
      if(verbose_) cout << legend << funclegend << "ERROR: channel " << channelname << " not contained in the model" << endl;
   }
};

void ModelConfigurator::setVar(std::string name, double value){
   std::string funclegend = "ModelConfigurator::setVar(std:::string name, double value) - ";
   if(verbose_) cout << legend << funclegend << "WARNING: manually adjusting value of variable "<< name << " from " << _CombinedWS->var(name.c_str())->getVal() << " to " << value << endl;
   _CombinedWS->var(name.c_str())->setVal(value);
}

void ModelConfigurator::setVarRange(std::string name, double valuelow, double valuehigh){
   std::string funclegend = "ModelConfigurator::setVarRange(std::string name, double valuelow, double valuehigh) - ";
   if(verbose_) cout << legend << funclegend << "WARNING: manually adjusting range of variable " << name << " to " << valuelow << " .. " << valuehigh << endl;
   _CombinedWS->Print();
   std::cout <<"getting var"<<std::endl;
   _CombinedWS->var(name.c_str())->setRange(valuelow, valuehigh);
}
