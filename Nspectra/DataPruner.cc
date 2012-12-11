//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "DataPruner.hh"
#include "TStopwatch.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <boost/lexical_cast.hpp>

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "TH1F.h"
#include "RooRealVar.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "TTree.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace boost;


//--------------------Constructor
DataPruner::DataPruner( std::map<std::string , double> Rangemap, ModelConfiguratorZprime * configurator) : legend("DataPruner- ") {
   cout << legend << "( std::map<std::string , double> Rangemap) constructed" << endl;

   _configurator = configurator;
   //_channelnames = _configurator->getChannelNames();
   _Rangemap = Rangemap;
   option = "lowerrange_customized";
 
}

//--------------------Constructor
DataPruner::DataPruner( double n_width, ModelConfiguratorZprime * configurator) : legend("DataPruner- ") {
   cout << legend << "( double n_width, ModelConfiguratorZprime * configurator) constructed" << endl;

   _n_width = n_width;
   _configurator = configurator;
   _channelnames = _configurator->getChannelNames();
   option = "range_by_signal_width";

}

//--------------------Destructor
DataPruner::~DataPruner() {
 ;
}

void DataPruner::Prune( std::map<string, RooDataSet*> * Datamap){
   cout << legend << "Prune( std::map<string, RooDataSet*> * Datamap)" << endl;
   _Datamap = Datamap;
   
   if (option == "lowerrange_customized"){
      for(std::map<std::string , double>::iterator mapIt = _Rangemap.begin(); mapIt != _Rangemap.end(); mapIt++ ){
         Prune((*mapIt).first, (*mapIt).second);
      }
   }

   if (option == "range_by_signal_width"){

      double massThreshLow = 20000.; 
      double massThreshHigh = 0.; 

      for(std::map<std::string , RooDataSet*>::iterator mapIt = _Datamap->begin(); mapIt != _Datamap->end(); mapIt++ ){

         //WARNING: ugly reproduction from WS, because I have no idea how to evaluate RooProduct or RooAddition
         std::string paramname1 = "sigma_";
         paramname1 += (*mapIt).first;
         cout << "getting parameter: " << paramname1 << endl;
         double param1 = ((RooProduct*)_configurator->getCombinedWS()->arg(paramname1.c_str()))->getVal();
         cout << "param1: " << param1 << endl;

         std::string paramname2 = "width_";
         paramname2 += (*mapIt).first;
         cout << "getting parameter: " << paramname2 << endl;
         double param2 = ((RooAddition*)_configurator->getCombinedWS()->arg(paramname2.c_str()))->getVal();
         cout << "param2: " << param2  << endl;

         double rangeparam = param1 + param2; // temporary solution; try to find s.th. better (e.g. on FWHM approximation of the Voigt distribution?)

         double massThreshLowTemp = _configurator->getMassHypothesis() - _n_width * rangeparam;
         double massThreshHighTemp = _configurator->getMassHypothesis() + _n_width * rangeparam;

         if ((massThreshLowTemp < massThreshLow) && (massThreshLowTemp > 200.)){ //CONVENTION: never go below 200 in Zprime search
            massThreshLow = massThreshLowTemp;
         }
         if ((massThreshHighTemp > massThreshHigh) && (massThreshHighTemp < 20000.)){
            massThreshHigh = massThreshHighTemp;
         }
      }

      for(std::map<std::string , RooDataSet*>::iterator mapIt = _Datamap->begin(); mapIt != _Datamap->end(); mapIt++ ){
         Prune((*mapIt).first, massThreshLow, massThreshHigh);
      }

      _configurator->setVarRange("mass", massThreshLow, massThreshHigh);  //not really clear to me why this should be necessary
   }

}

void DataPruner::Prune(std::string channelname, double massThreshold){

   std::string funclegend = " Prune(std::string channelname, double massThreshold) ";

   //create constrained dataset
   //const TTree * mytree = (_configurator->getChannelDataMap().find(channelname)->second->tree());
   //TTree * thetree = (TTree*) mytree->Clone(mytree->GetName());
   //TTree thetree = *mytree;

   //CONVENTION: observable "mass" is hardcoded
   std::string optionstring = "mass>";
   optionstring += lexical_cast<std::string>(massThreshold);

   cout << legend << funclegend << "events in channel channel" << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " before pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl ;
   cout << legend << funclegend << "pruning channel" << channelname << " with restriction " << optionstring << endl;

   RooDataSet* prunedData = (RooDataSet *) _Datamap->find(channelname)->second->reduce(optionstring.c_str()) ; //COMMENT: is this cast evil? (according to the Roofit manual 2.91 p. 80 it should not be necessary)

   //FIXME: potential memory leak in lines below, check if old dataset needs to be deleted

   //pointer to old dataset
   RooDataSet * pointerfordelete = _Datamap->find(channelname)->second;
   //store new dataset
   _Datamap->find(channelname)->second = prunedData;
   cout << legend << funclegend << "events in channel " << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " after pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl; ;

   //delete pointer to old dataset
   delete pointerfordelete;

}

void DataPruner::Prune(std::string channelname, double massThreshLow, double massThreshHigh){

   std::string funclegend = " Prune(std::string channelname, double massThreshold) ";

   cout << legend << funclegend << "events in channel channel" << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " before pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl ;

   //CONVENTION: observable "mass" is hardcoded
   std::string optionstring = "mass>";
   optionstring += lexical_cast<std::string>(massThreshLow);
   optionstring += " && mass<";
   optionstring += lexical_cast<std::string>(massThreshHigh);

   cout << legend << funclegend << "pruning channel" << channelname << " with restriction " << optionstring << endl;

   RooDataSet* prunedData = (RooDataSet *) _Datamap->find(channelname)->second->reduce(optionstring.c_str()) ; //COMMENT: is this cast evil? (according to the Roofit manual 2.91 p. 80 it should not be necessary)

   //FIXME: potential memory leak in lines below, check if old dataset needs to be deleted

   //pointer to old dataset
   RooDataSet * pointerfordelete = _Datamap->find(channelname)->second;
   //store new dataset
   _Datamap->find(channelname)->second = prunedData;
   cout << legend << funclegend << "events in channel " << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " after pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl; ;

   //CONVENTION: the following line relies on a hardcoded parameter name
   std::string varstring = "nbkg_est_";
   varstring += channelname;
   _configurator->setVar(varstring.c_str(), _Datamap->find(channelname)->second->sumEntries());

   //delete pointer to old dataset
   delete pointerfordelete;

}

