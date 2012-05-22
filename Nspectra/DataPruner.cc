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
#include "TTree.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace boost;


//--------------------Constructor
DataPruner::DataPruner( std::map<std::string , double> Rangemap) : legend("DataPruner- ") {
   cout << legend << " constructed" << endl;

   //_configurator = configurator;
   //_channelnames = _configurator->getChannelNames();
   _Rangemap = Rangemap;
 
}

DataPruner::~DataPruner() {
 ;
}

void DataPruner::Prune( std::map<string, RooDataSet*> * Datamap){
   _Datamap = Datamap;
   for(std::map<std::string , double>::iterator mapIt = _Rangemap.begin(); mapIt != _Rangemap .end(); mapIt++ ){
      Prune((*mapIt).first, (*mapIt).second);
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
   cout << legend << funclegend << "events in channel channel" << channelname << " at address:" << &(_Datamap->find(channelname)->second) << " after pruning: " << _Datamap->find(channelname)->second->sumEntries() << endl; ;

   //delete pointer to old dataset
   delete pointerfordelete;

}

