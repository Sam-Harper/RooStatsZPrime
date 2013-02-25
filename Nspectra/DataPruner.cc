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

   _eventnum = -1.;
   _configurator = configurator;
   //_channelnames = _configurator->getChannelNames();
   _Rangemap = Rangemap;
   _option = "lowerrange_customized";
 
}

//--------------------Constructor
DataPruner::DataPruner( double n_width, ModelConfiguratorZprime * configurator, double thresholdLow_min, double eventnum) : legend("DataPruner- ") {
   cout << legend << "( double n_width, ModelConfiguratorZprime * configurator) constructed" << endl;

   _n_width = n_width;
   _configurator = configurator;
   _channelnames = _configurator->getChannelNames();
   _option = "range_by_signal_width";
   _thresholdLow_min = thresholdLow_min;
   _eventnum = eventnum;
}

//--------------------Destructor
DataPruner::~DataPruner() {
 ;
}

void DataPruner::Prune( std::map<string, RooDataSet*> * Datamap){
   cout << legend << "Prune( std::map<string, RooDataSet*> * Datamap)" << endl;
   _Datamap = Datamap;
   
   if (_option == "lowerrange_customized"){
      for(std::map<std::string , double>::iterator mapIt = _Rangemap.begin(); mapIt != _Rangemap.end(); mapIt++ ){
         Prune((*mapIt).first, (*mapIt).second);
      }
   }

   if (_option == "range_by_signal_width"){

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
         cout << "minimum lower mass threshold is minimum of: " << massThreshLow  << " , " << _thresholdLow_min  << endl;
         massThreshLow = std::min( massThreshLow, _thresholdLow_min );
         cout << "minimum lower mass threshold is: " << massThreshLow  << endl;
         //massThreshLow = std::min( massThreshLow, _thresholdLow_min_map.find((*mapIt).first)->second );
         //cout << "minimum lower mass threshold for channel: " << (*mapIt).first << " is: " <<  _thresholdLow_min_map.find((*mapIt).first)->second << endl;
      }



      for(std::map<std::string , RooDataSet*>::iterator mapIt = _Datamap->begin(); mapIt != _Datamap->end(); mapIt++ ){
         Prune((*mapIt).first, massThreshLow, massThreshHigh);
      }

      _configurator->setVarRange("mass", massThreshLow, massThreshHigh);  //not fully clear to me why this is necessary
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

   cout << legend << funclegend << "events in channel " << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " before pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl ;

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

//    //possibility 1: normalize to number of observed events in the considered mass range
//    //CONVENTION: the following lines rely on a hardcoded parameter name
//    std::string varstring = "nbkg_est_";
//    varstring += channelname;
//    _configurator->setVar(varstring.c_str(), _Datamap->find(channelname)->second->sumEntries());

   //possibility 2: normalize to number of observed events in the considered mass range
   //CONVENTION: the following lines rely on a hardcoded name for the background pdf and several parameters
   std::string varstring = "nbkg_est_";
   varstring += channelname;
   std::string pdfstring = "bkgpdf_";
   pdfstring += channelname;
   double varstring_value = _configurator->getCombinedWS()->var(varstring.c_str())->getVal();
   RooRealVar* x = _configurator->getCombinedWS()->var("mass");
   x->setRange("pruned_range", massThreshLow, massThreshHigh); // create range to integrate over
   RooAbsReal* i = _configurator->getCombinedWS()->pdf(pdfstring.c_str())->createIntegral(*x, RooFit::NormSet(*x), RooFit::Range("pruned_range"));
   double pdfstring_value = i->getVal();
   std::cout << "Integral value of " << pdfstring << " in pruned range: " << pdfstring_value << std::endl;
   _configurator->setVar(varstring.c_str(), (varstring_value*pdfstring_value) );
   //_configurator->setVar(varstring.c_str(), (varstring_value*pdfstring_value*1.2) ); //FIXME: factor 1.2 only for testing

   //delete pointer to old dataset
   delete pointerfordelete;

}

void DataPruner::Norm(std::map<string, RooDataSet*> * Datamap){

   if(_eventnum < 0.){return;}

   std::string funclegend = " Norm(double eventnum) ";

   for(std::map<std::string , RooDataSet*>::iterator mapIt = Datamap->begin(); mapIt != Datamap->end(); mapIt++ ){
         std::string channelname = (*mapIt).first;
         cout << legend << funclegend << "events in channel " << channelname << " at address:" << &(_Datamap->find(channelname)->second) << Datamap->find(channelname)->second->sumEntries() << endl ;
         RooDataSet* myRooDataSet = Datamap->find(channelname)->second;

         //find current mass range
         //CONVENTION: name "mass" hardcoded
         RooRealVar* x = _configurator->getCombinedWS()->var("mass");
         double mass_window_lower_border = x->getMin();
         double mass_window_upper_border = x->getMax();

         double norm_window_lower_border = mass_window_lower_border;
         double norm_window_upper_border = mass_window_upper_border;

         //find mass range with 400 events
        double events_in_norm_window = -1;
         while (true){

         std::string optionstring = "mass>";
         optionstring += lexical_cast<std::string>(norm_window_lower_border);
         optionstring += " && mass<";
         optionstring += lexical_cast<std::string>(norm_window_upper_border);
         
         events_in_norm_window = myRooDataSet->sumEntries(optionstring.c_str()); 
         if (events_in_norm_window < _eventnum ){
            break;
         }
         else{
            cout << legend << funclegend << "norm channel" << channelname << " with restriction " << optionstring << " : "<< events_in_norm_window <<endl;
            norm_window_lower_border = norm_window_lower_border + 0.1; //0. GeV steps -> may bot be ok for all applications
            norm_window_upper_border = norm_window_upper_border - 0.1;
         }
      }

      //prepare readjusting normalization

      std::string varstring = "nbkg_est_";
      varstring += channelname;
      std::string pdfstring = "bkgpdf_";
      pdfstring += channelname;

      double varstring_value = _configurator->getCombinedWS()->var(varstring.c_str())->getVal();

      //find bkg norm in mass window

      RooRealVar* x1 = _configurator->getCombinedWS()->var("mass");
      x1->setRange("full_range", mass_window_lower_border, mass_window_upper_border); // create range to integrate over
      RooAbsReal* i1 = _configurator->getCombinedWS()->pdf(pdfstring.c_str())->createIntegral(*x1, RooFit::NormSet(*x1), RooFit::Range("full_range"));
      double norm_mass_window = i1->getVal();
      std::cout << "Integral value of " << pdfstring << " in full range: " << norm_mass_window << std::endl;
   
      //find bkg norm. in norm window

      RooRealVar* x2 = _configurator->getCombinedWS()->var("mass");
      x2->setRange("norm_range", norm_window_lower_border, norm_window_upper_border); // create range to integrate over
      RooAbsReal* i2 = _configurator->getCombinedWS()->pdf(pdfstring.c_str())->createIntegral(*x2, RooFit::NormSet(*x2), RooFit::Range("norm_range"));
      double norm_norm_window = i2->getVal();
      std::cout << "Integral value of " << pdfstring << " in norm range: " << norm_norm_window << std::endl;

      // readjust normalization

      double new_norm = varstring_value * (events_in_norm_window / ((norm_norm_window / norm_mass_window ) * varstring_value ) );
      cout << "old norm: " << varstring_value << " new norm: " << new_norm << endl;
      _configurator->setVar(varstring.c_str(), new_norm );
      cout << "parameter " << varstring << "has been readjusted to: " << _configurator->getCombinedWS()->var(varstring.c_str())->getVal() << endl;
   }


}

