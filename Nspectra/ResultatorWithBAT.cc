//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "ResultatorWithBAT.hh"
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
#include "BAT/BATCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/SequentialProposal.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

#include "TCanvas.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;


//--------------------Constructor
ResultatorWithBAT::ResultatorWithBAT(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner, std::string plotfile) : sublegend("WithBAT - "), Resultator(configurator, myDataPruner, plotfile) {
   ;
}

ResultatorWithBAT::~ResultatorWithBAT() {
   delete _databox;
}

void ResultatorWithBAT::calculateMCMClimit( UInt_t mcmc_iter, UInt_t mcmc_burnin, UInt_t ntoys, std::string suffix, std::string mode){
   std::string funclegend = "calculateMCMClimit( Int_t mcmc_iter, Int_t mcmc_burnin, Int_t nruns, std::string suffix, std::string mode) - ";
   cout << legend << funclegend << endl;


   // time it
   TStopwatch timer;
   timer.Start();

   unsigned int toycounter = 1;

   while(toycounter <= ntoys){

      //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   
      //Double_t limit = -1.0;  
      
      // for mass limit, add k-factor systematics to the nsig systematics
      // FIXME: this is a correlated part of the uncertainty!!!
      //  - different uncertainties for graviton and Z' models
   //     if ( mode.find("mass_") != std::string::npos ){
   //       
   //       std::cout << _legend << std::endl;
   //       std::cout << _legend << "this a mass limit calculation," << std::endl;
   //       std::cout << _legend << "I will add k-factor uncertainty to the nsig uncertainty" << std::endl;
   //       std::cout << _legend << "I will do it " << ntoys << " times, so one can average. " << pe_counter+1 << " of " << ntoys << std::endl;
   //       std::cout << _legend << std::endl;
   // 
   //       Double_t kfactor_err = GetKfactorUncertainty(peak, mode);
   // 
   //       double nsig_kappa = ws->var("nsig_kappa_dimuon")->getVal();
   //       nsig_kappa = 1.0+sqrt((nsig_kappa-1.0)*(nsig_kappa-1.0) + kfactor_err*kfactor_err);
   //       ws->var("nsig_kappa_dimuon")->setVal(nsig_kappa);
   // 
   //       nsig_kappa = ws->var("nsig_kappa_dielectron")->getVal();
   //       nsig_kappa = 1.0+sqrt((nsig_kappa-1.0)*(nsig_kappa-1.0) + kfactor_err*kfactor_err);
   //       ws->var("nsig_kappa_dielectron")->setVal(nsig_kappa);
   // 
   //       //ntoys = 1;
   //     }
   
   
      // calculate observed ratio limit
   
   
      // change POI range
      //delete mpMc;
      //mpMc = new ModelConfig(mc);
      //mpMc->SetWorkspace(*ws);
      //double poiUpper = GetPoiUpper("dilepton", peak);
      //double poiUpper = _poiUpperRange;
      std::cout << legend << funclegend << "setting POI range to [0; " << _poiUpperRange << "]" << std::endl;
   
      _myWS->var("ratio")->setRange(0.0, _poiUpperRange); //FIXME variable name does not need to be hardcoded
   
      // study of the background normalization
      //StudyBackground();
      RooDataSet * data;
   
      if(mode.find("observed")!=std::string::npos){

         std::cout << legend << funclegend << "calculating an OBSERVED limit..." << std::endl;
         std::cout << legend << funclegend << "For stability reasons we calculate it " << ntoys << " times, so one can average. This is iteration " << toycounter << " of " << ntoys << std::endl;

         data = _databox->createObservedData();
         _mcInt = GetMcmcInterval(*_myModelConfig, data, 0.95, mcmc_iter, mcmc_burnin, 0.0, _nbinsPosterior);
         toycounter++; 
      }
      else{

         std::cout << legend << funclegend << "calculating an EXPECTED limit..." << std::endl;
         std::cout << legend << funclegend << "We calculate it for " << ntoys << " toy experiments. This is iteration " << toycounter << " of " << ntoys << std::endl;

         data = _databox->createToyMc();
         _mcInt = GetMcmcInterval(*_myModelConfig, data, 0.95, mcmc_iter, mcmc_burnin, 0.0, _nbinsPosterior);
         toycounter++; 
      }

      if(_bMcmcConverged){
         std::string _outfile = "nchannels_ratio_mcmc_limit_" + suffix + ".ascii";
         printMcmcUpperLimit( _configurator->getMassHypothesis(), _outfile );
      }
      else{
//FIXME: lines below are not the solution because convergence to the right value is the problem and not the poi range
//         std::cout << legend << funclegend << "WARNING: automatically adjusting POI range from" << _poiUpperRange << " to " << (_poiUpperRange*2.) << std::endl;
//         _poiUpperRange = _poiUpperRange*2.;
           toycounter--;
      }
   
      // make extra plots
      //if (plot_name.size() > 0){
      //  char buf[356];
      //  std::string _outfile = suffix + "_%d_" + plot_name;
      //  sprintf(buf, _outfile.c_str(), pe_counter);
      //  makeMcmcPosteriorPlot( buf );
      //}

      delete data; //COMMENT: that should be the right place to do this, but check if there are memory problems

   }

   // get time
   timer.Print();

}

MCMCInterval * ResultatorWithBAT::GetMcmcInterval(ModelConfig mc, RooDataSet * data,
               double conf_level,
               int n_iter,
               int n_burn,
               double left_side_tail_fraction,
               int n_bins){
  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit

  std::string funclegend = " GetMcmcInterval(...) ";
  cout << legend << funclegend << endl;

  MCMCInterval * mcInt;

        // prepare BATCalulator
      //BATCalculator batcalc(data, modelconfig);
      BATCalculator* batcalc = new BATCalculator(*data, mc, true );
      cout << "setting name for BATCalculator" << endl;
      // give the BATCalculator a unique name (always a good idea in ROOT)
      TString namestring = "mybatc_";
      //namestring += obs;
      batcalc->SetName(namestring);

      cout << "fixing ConfidenceLevel" << endl;
      // fix amount of posterior probability in the calculated interval.
      // FIXME: C.L. hardcoded (use conf_level argument of the function)
      batcalc->SetConfidenceLevel(0.90);
      // fix number of Markov chain elements. (in general: the longer the Markov chain the more 
      // precise will be the results)
      batcalc->SetnMCMC(n_iter);

      batcalc->SetNumBins("ratio", 200);

      // retrieve SimpleInterval object containing the information about the interval (this 
      // triggers the actual calculations)
//      cout << "constructing interval object ..." << endl;
//      SimpleInterval * interval = batcalc->GetInterval1D();
//      double LowerLim = interval->LowerLimit();
//      double UpperLim = interval->UpperLimit();

//      std::cout << "BATCalculator: 90% CL interval: [ " << LowerLim << " - " << UpperLim << " ] or 95% CL upper limit\n";


      mcInt = batcalc->GetInterval();


  cout << "test 4" << endl;

  // check if limit makes sense
  _bMcmcConverged = false; // default
  if (mcInt){
    RooRealVar * p_first_poi = (RooRealVar*) mc.GetParametersOfInterest()->first();
    double poi_limit = mcInt->UpperLimit(*p_first_poi);
    double u_poi_min  = p_first_poi->getMin();
    double u_poi_max  = p_first_poi->getMax();
    double u_poi_gap = (u_poi_max-poi_limit)/(u_poi_max-u_poi_min);
    std::cout << legend << "POI upper limit: " << poi_limit << std::endl;
    std::cout << legend << "POI range: [" << u_poi_min << ", " << u_poi_max << "]" << std::endl;
    std::cout << legend << "POI upper gap (fraction of range): " << u_poi_gap << std::endl;
    if (u_poi_gap<0.2){
      std::cout << legend << funclegend
      << "POI limit too close to the upper boundary, MCMC probably failed!!!" << std::endl;
      std::cout << legend << funclegend
      << "returning interval and setting fail flag" << std::endl;
      _bMcmcConverged = false;
    }
    else{
      _bMcmcConverged = true;
    }
  }
  else std::cout << legend << funclegend << "No interval found!" << std::endl;

  if(_plotfile != ""){
    TFile* outfile = new TFile(_plotfile.c_str(),"RECREATE");
    outfile->WriteTObject(batcalc->GetPosteriorPlot1D());
  }

  batcalc->CleanCalculatorForNewData();
  delete batcalc; //COMMENT: keep an eye on this
  
  return mcInt;

  cout << "test end" << endl;

}

double ResultatorWithBAT::printMcmcUpperLimit( double peak, std::string filename ){
  //
  // print out the upper limit on the first Parameter of Interest
  //

  char buf[1024];

  double _limit = numeric_limits<double>::max();

  if (_mcInt){
    RooRealVar * firstPOI = (RooRealVar*) (*_myModelConfig).GetParametersOfInterest()->first();
    _limit = _mcInt->UpperLimit(*firstPOI);
    std::cout << "\n95% upper limit on " << firstPOI->GetName() << " is : " <<
      _limit << endl;
    
    if (_bMcmcConverged){
      sprintf(buf, "%7.1f   %7.6f", peak, _limit);
    }
    else{
      sprintf(buf, "# %7.1f   %7.6f  # MCMC did not converge", peak, _limit);
    }

  }
  else{
    sprintf(buf, "# MCMC did not converge");
  }

  if (filename.size()!=0){
    
    std::ofstream aFile;
    
    // append to file if exists
    aFile.open(filename.c_str(), std::ios_base::app);
    
    aFile << buf << std::endl;
    
    // close outfile here so it is safe even if subsequent iterations crash
    aFile.close();
    
  }
    
  return _limit;
}

