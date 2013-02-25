//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "DataBox.hh"

#include "TStopwatch.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>

#include "RooStats/SimpleInterval.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/SequentialProposal.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooRandom.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;


//--------------------Constructor
DataBox::DataBox(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner) : legend("DataBox - ")  {
   cout << legend << " constructed" << endl;
   _configurator = configurator;
   _myDataPruner = myDataPruner;
   _myWS = _configurator->getCombinedWS() ;
   _myModelConfig = _configurator->getCombinedModelConfig();
   _channelnames = _configurator->getChannelNames();
   _pruned = false;

   SetNeventsData();

   _myrandgen = new TRandom3(0);
   RooRandom::randomGenerator()->SetSeed(0); //otherwise we will get the same pseudo experiments in all jobs
}

DataBox::~DataBox() {
;
}

RooDataSet * DataBox::createObservedData(){
   
   // create combined dataset
   RooRealVar * mass = _myWS->var("mass"); //FIXME: can be done without hardcoding the parameter name
   
   //prune data if requested
   if(!_pruned){
   _myDataPruner->Prune(_configurator->getChannelDataMapPointer());
   _myDataPruner->Norm(_configurator->getChannelDataMapPointer());
   _pruned = true; //avoid multiple pruning leading to wrong normalization
   }
   
   cout << "creating combined dataset ... " << endl;
   
   RooDataSet * data = new RooDataSet( "combdata", "combined data",
            RooArgSet((const RooRealVar)(*mass) ),
            Index(*(_configurator->get_channel_cat())),
            Import(_configurator->getChannelDataMap()) );
   
   data->Print();
   cout << "... created combined dataset " << endl;
   
   return data;
}

// toy experiments for limit calcualtion
RooDataSet * DataBox::createToyMc(){
   std::string funclegend = "createToyMc() - ";
   cout << legend << funclegend << endl;

   // check if previous toy data needs to be cleared -> if yes, do it properly
   if(_DatamapToy.empty() == false){
      for(std::map<string, RooDataSet*>::iterator mapit = _DatamapToy.begin(); mapit !=  _DatamapToy.end() ; mapit++){
         delete mapit->second;
      }
      _DatamapToy.clear();
   }

   //   loop over channels
   for(std::vector<std::string>::iterator nameit = _channelnames.begin(); nameit != _channelnames.end(); nameit++ ){

      //CONVENTION: variable names in the workspace must  be conform to the following few lines
      std::string syst_nbkg_channelname = "syst_nbkg_" + (*nameit);
      std::string beta_nbkg_channelname = "beta_nbkg_" + (*nameit);
      std::string nbkg_kappa_channelname = "nbkg_kappa_" + (*nameit);
      std::string nbkg_est_channelname = "nbkg_est_" + (*nameit);

      Double_t _ntoy = -1.;
      Double_t _beta = -1.;
      Double_t _kappa = -1.;
      Double_t _nbkg_est = -1.;

      if(((*nameit) == "dielectron_ebeb") || ((*nameit) == "dielectron_ebee")) {

         std::string nbkg_kappa_b_channelname = "nbkg_kappa_b_" + (*nameit);
         std::string nbkg_kappa_m_channelname = "nbkg_kappa_m_" + (*nameit);
         std::string nbkg_kappa_msquared_channelname = "nbkg_kappa_msquared_" + (*nameit);
         std::string peak = "peak";
         Double_t nbkg_kappa_b_channelname_val = _myWS->var(nbkg_kappa_b_channelname.c_str())->getVal();
         Double_t nbkg_kappa_m_channelname_val = _myWS->var(nbkg_kappa_m_channelname.c_str())->getVal();
         Double_t nbkg_kappa_msquared_channelname_val = _myWS->var(nbkg_kappa_msquared_channelname.c_str())->getVal();
         Double_t peak_val = _myWS->var(peak.c_str())->getVal();

         _beta = GetRandom(syst_nbkg_channelname , beta_nbkg_channelname);
         _kappa = nbkg_kappa_b_channelname_val + nbkg_kappa_m_channelname_val * peak_val + nbkg_kappa_msquared_channelname_val * peak_val * peak_val;
         _nbkg_est = _myWS->var(nbkg_est_channelname.c_str())->getVal();
         _ntoy = pow(_kappa,_beta) * _nbkg_est;
      }
      else{
         _beta = GetRandom(syst_nbkg_channelname , beta_nbkg_channelname);
         _kappa = _myWS->var(nbkg_kappa_channelname.c_str())->getVal();
         _nbkg_est = _myWS->var(nbkg_est_channelname.c_str())->getVal();
         _ntoy = pow(_kappa,_beta) * _nbkg_est;
      }
      
      Int_t _n = _myrandgen->Poisson(_ntoy);
      
      // all nuisance parameters:
      //   beta_nsig_dimuon, 
      //   beta_nbkg_dimuon,
      //   lumi_nuis
      
      // create dataset
      RooRealVar * _mass = _myWS->var("mass"); //FIXME: can be done without hardcoding the parameter name
      RooArgSet _vars(*_mass);
      
      //CONVENTION: pdf name for the background description in the workspace must be conform to the following line
      std::string bkgpdf_channelname = "bkgpdf_" + (*nameit);
      RooAbsPdf * _pdf = _myWS->pdf(bkgpdf_channelname.c_str());
      
      RooAbsPdf::GenSpec * _spec = _pdf->prepareMultiGen(_vars,
                           Name("toys"),
                           NumEvents(_n),
                           Extended(kFALSE),
                           Verbose(kTRUE));

      RooDataSet * toydata = _pdf->generate(*_spec); 
      
      delete _spec;
      
      _DatamapToy.insert( pair<string, RooDataSet *>((*nameit), toydata) );
      
      Int_t n_generated_entries = (Int_t)(toydata->sumEntries());
      
      // debug output
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! _beta = " << _beta << std::endl;
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! _kappa = " << _kappa << std::endl;
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! _nbkg_est = " << _nbkg_est << std::endl;
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! _ntoy     = " << _ntoy << std::endl;
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! _n        = " << _n    << std::endl;
      std::cout << legend << funclegend << (*nameit) << " - !!!!!!!!!!!!!! n_generated_entries = " << n_generated_entries    << std::endl;
      
   }
   
   RooRealVar * mass = _myWS->var("mass"); //FIXME: can be done without hardcoding the parameter name

   //Prune data (remove all events below specified mass threshold) if requested in the setup
   if(!_pruned){
   _myDataPruner->Prune(&_DatamapToy);
   _myDataPruner->Norm(&_DatamapToy);
   _pruned = true; //avoid multiple pruning leading to wrong normalization
   }

   // create combined dataset
   RooDataSet * data = new RooDataSet( "combdata", "combined data",
            RooArgSet((const RooRealVar)(*mass) ),
            Index(*(_configurator->get_channel_cat())),
            Import(_DatamapToy) );
   
   data->Print();
   
   return data;
}


//COMMENT: is this the best way to do this?; move to Pixie (extended with workspace as additional argument)?
Double_t DataBox::GetRandom( std::string pdf, std::string var ){
  //
  // generates a random number using a pdf in the workspace
  //
  
  // generate a dataset with one entry
  RooDataSet * _ds = _myWS->pdf(pdf.c_str())->generate(*_myWS->var(var.c_str()), 1);

  Double_t _result = ((RooRealVar *)(_ds->get(0)->first()))->getVal();
  delete _ds;

  return _result;
}

// toy experiments for significance calculation
//except for the content of the loop over the channels the structure is identical createToyMc(), so there is room for an improved code structure
RooDataSet * DataBox::createPeData(){
  //
  // Creates background-only pseudodata 
  // which mimics the observed data.
  // Number of events exactly matches the observed number
  //

   std::string funclegend = "createPeData() - ";
   cout << legend << funclegend << endl;

   // check if previous toy data needs to be cleared -> if yes, do it properly
   if(_DatamapToy.empty() == false){
      for(std::map<string, RooDataSet*>::iterator mapit = _DatamapToy.begin(); mapit !=  _DatamapToy.end() ; mapit++){
         delete mapit->second;
      }
      _DatamapToy.clear();
   }

   //   loop over channels
   for(std::vector<std::string>::iterator nameit = _channelnames.begin(); nameit != _channelnames.end(); nameit++ ){

      Int_t nchannel = _neventsDataMap.find(*nameit)->second;
      
      // create dataset
      RooRealVar * _mass = _myWS->var("mass"); //FIXME: no need to have the name of the observable hardcoded
      RooArgSet _vars(*_mass);
      
      //CONVENTION: pdf name for the background description in the workspace must be conform to the following line
      std::string bkgpdf_channelname = "bkgpdf_" + (*nameit);
      RooAbsPdf * _pdf = _myWS->pdf(bkgpdf_channelname.c_str());
      
      RooAbsPdf::GenSpec * _spec = _pdf->prepareMultiGen(_vars,
                           Name("toys"),
                           NumEvents(nchannel),
                           Extended(kFALSE),
                           Verbose(kTRUE));

      RooDataSet * toydata = _pdf->generate(*_spec); 
      
      delete _spec;
      
      _DatamapToy.insert( pair<string, RooDataSet *>((*nameit), toydata) );
      
      Int_t n_generated_entries = (Int_t)(toydata->sumEntries());
      
      // debug
      //std::cout << "!!!!!!!!!!!!!! _n      = " << nchannel    << std::endl;
      
      //return n_generated_entries;

   }

   RooRealVar * mass = _myWS->var("mass"); //FIXME: can be done without hardcoding the parameter name

   //Prune data (remove all events below specified mass threshold) if requested in the setup
   _myDataPruner->Prune(&_DatamapToy);
   _myDataPruner->Norm(&_DatamapToy);

   // create combined dataset
   RooDataSet * data = new RooDataSet( "combdata", "combined data",
            RooArgSet((const RooRealVar)(*mass) ),
            Index(*(_configurator->get_channel_cat())),
            Import(_DatamapToy) );
   
   data->Print();

   return data;
}

void DataBox::SetNeventsData(){

   //   loop over channels
   for(std::vector<std::string>::iterator nameit = _channelnames.begin(); nameit != _channelnames.end(); nameit++ ){
       _neventsDataMap.insert( pair<string, int>((*nameit), static_cast<int>((_configurator->getChannelDataMap())[(*nameit)]->sumEntries()) ) );
   }

}
