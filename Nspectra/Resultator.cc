//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "Resultator.hh"
#include "TStopwatch.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <fstream>

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

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "PoiRangeEstimator.hh"

using namespace std;
using namespace RooFit;
using namespace RooStats;


//--------------------Constructor
Resultator::Resultator(ModelConfiguratorZprime* configurator, DataPruner * myDataPruner, std::string plotfile) : legend("Resultator - "), _plotfile(plotfile) {
   cout << legend << " constructed" << endl;
   _configurator = configurator;
   _myDataPruner = myDataPruner;
   _poiUpperRange = -1.;
   _myWS = _configurator->getCombinedWS() ;
   _myModelConfig = _configurator->getCombinedModelConfig();
   _bMcmcConverged = false;
   _low_POI_edge_problem = false;
   _high_POI_edge_problem = false;
   _databox = new DataBox(configurator, myDataPruner);
   _nbinsPosterior = 100;

   _makeLlrPlot = false;
   _poirangefactor = 1.0;
   // study of the background normalization
   //StudyBackground();

   //ModelConfig mc = *_myModelConfig;

}

Resultator::~Resultator() {
   delete _databox;
}

void Resultator::calculateMCMClimit( UInt_t mcmc_iter, UInt_t mcmc_burnin, UInt_t ntoys, std::string suffix, std::string mode){
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
   
   
      std::cout << legend << funclegend << "setting POI range to [0; " << _poiUpperRange << "]" << std::endl;
   
      _myWS->var("ratio")->setRange(0.0, _poiUpperRange); //FIXME variable name does not need to be hardcoded
   
      // study of the background normalization
      //StudyBackground();
      RooDataSet * data;
   
      if(mode.find("observed")!=std::string::npos){
         std::cout << legend << funclegend << "calculating an OBSERVED limit..." << std::endl;
         std::cout << legend << funclegend << "For stability reasons we calculate it " << ntoys << " times, so one can average. This is iteration " << toycounter << " of " << ntoys << std::endl;
         data = _databox->createObservedData();
      }
      else{
         std::cout << legend << funclegend << "calculating an EXPECTED limit..." << std::endl;
         std::cout << legend << funclegend << "We calculate it for " << ntoys << " toy experiments. This is iteration " << toycounter << " of " << ntoys << std::endl;
         data = _databox->createToyMc();
         
      }

      //Estimate reasonable POI range -> do this again here to cover changed data in pseudo experiments
      PoiRangeEstimator * myPoiRangeEstimator = new PoiRangeEstimator(_configurator, data);
      double poiUpperLimitGuess = myPoiRangeEstimator->GetPoiUpper("multi", _configurator->getMassHypothesis());
      cout << "estimate for reasonable upper limit of poi range after data generation: " << poiUpperLimitGuess << endl; 
      delete myPoiRangeEstimator;
      SetPoiUpperLimitByHand(poiUpperLimitGuess*_poirangefactor); // times two is of course just an ad-hoc number

      _mcInt = GetMcmcInterval(*_myModelConfig, data, 0.95, mcmc_iter, mcmc_burnin, 0.0, _nbinsPosterior);
      toycounter++; 

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

MCMCInterval * Resultator::GetMcmcInterval(ModelConfig mc, RooDataSet * data,
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

  /*
  RooFitResult* fit = ws->pdf("model")->fitTo(*data,Save());
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction* pf = ph.GetProposalFunction();
  */

  MCMCInterval * mcInt;

  // FIXME: testing: this proposal function seems fairly robust
  SequentialProposal sp(0.5);

  MCMCCalculator mcmc( *data, mc );
  mcmc.SetConfidenceLevel(conf_level);
  mcmc.SetNumIters(n_iter);          // Metropolis-Hastings algorithm iterations

  // FIXME: testing: different proposal function
  //mcmc.SetProposalFunction(*pf);
  mcmc.SetProposalFunction(sp);

  mcmc.SetNumBurnInSteps(n_burn); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(left_side_tail_fraction);
  mcmc.SetNumBins(n_bins);

  cout << "test 1" << endl;
  
  // FIXME: testing good initial values - don't seem to do anything different
  //ws->var("ratio")->setVal(0.01);
  //ws->var("beta_nsig_dielectron")->setRange(-3.0, 3.0);
  //ws->var("beta_nbkg_dielectron")->setRange(-3.0, 3.0);
  //ws->var("beta_mass_dielectron")->setRange(-3.0, 3.0);

//mcInt = mcmc.GetInterval();
//  try {
//    delete mcInt;
//    cout << "test 2" << endl;
    mcInt = mcmc.GetInterval();
//  } catch ( std::length_error &ex) {
//    cout << "test 3" << endl;
//    mcInt = 0;
//  }

  cout << "test 4" << endl;

  // check if limit makes sense
  _bMcmcConverged = false; // default
  _low_POI_edge_problem = false;
  _high_POI_edge_problem = false;
  if (mcInt){
    //plot posterior info if requested
    if(_plotfile != ""){
      makeMcmcPosteriorPlot(mcInt, _plotfile);
    }
    RooRealVar * p_first_poi = (RooRealVar*) mc.GetParametersOfInterest()->first();
    double poi_limit = mcInt->UpperLimit(*p_first_poi);
    double poi_min  = p_first_poi->getMin();
    double poi_max  = p_first_poi->getMax();
    double u_poi_gap = (poi_max - poi_limit)/(poi_max - poi_min);
    double l_poi_gap = (poi_limit - poi_min)/(poi_max - poi_min);
    std::cout << legend << "POI upper limit: " << poi_limit << std::endl;
    std::cout << legend << "POI range: [" << poi_min << ", " << poi_max << "]" << std::endl;
    std::cout << legend << "POI upper gap (fraction of range): " << u_poi_gap << std::endl;
    std::cout << legend << "POI lower gap (fraction of range): " << l_poi_gap << std::endl;
    if (u_poi_gap<0.2){
      std::cout << legend << funclegend
      << "POI limit too close to the upper boundary, MCMC probably failed!!!" << std::endl;
      std::cout << legend << funclegend
      << "returning interval and setting fail flags" << std::endl;
      _bMcmcConverged = false;
      _high_POI_edge_problem = true;
    }
    else if (l_poi_gap<0.1){
      std::cout << legend << funclegend
      << "POI limit too close to the lower boundary, MCMC probably failed!!!" << std::endl;
      std::cout << legend << funclegend
      << "returning interval and setting fail flags" << std::endl;
      _bMcmcConverged = false;
        _low_POI_edge_problem = false;
    }
    else {
      _bMcmcConverged = true;
    }
  }
  else std::cout << legend << funclegend << "No interval found!" << std::endl;
  
  return mcInt;

  cout << "test end" << endl;

}

double Resultator::printMcmcUpperLimit( double peak, std::string filename ){
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

void Resultator::printDileptonRatioSignificance( std::string filename,
                     Double_t peak,
                     Double_t pllr){
  //
  // print out the test statistic and the corresponding resonance mass
  //

  if (filename.size()!=0){
    
    std::ofstream aFile;

    // append to file if exists
    aFile.open(filename.c_str(), std::ios_base::app);

    char buf[1024];
    sprintf(buf, "%7.1f   %7.6f", peak, pllr);

    aFile << buf << std::endl;

    // close outfile here so it is safe even if subsequent iterations crash
    aFile.close();

  }

}

Double_t Resultator::calculateRatioSignificance( std::string mode,
                    Double_t mass_low,
                    Double_t mass_high,
                    Double_t mass_step,
                    Int_t ntoys,
                    std::string suffix,    // suffix for output file names
                    Double_t fit_range_low,
                    Double_t fit_range_high,
                    Double_t width_factor,
                    bool makeLlrPlot // plot option 
                    ){ 
   //
   // Significance for xsec(Z')/xsec(Z) 
   //

   _makeLlrPlot = makeLlrPlot;

  std::string funclegend = " calculateRatioSignificance(...) ";
  cout << legend << funclegend << endl;

  // settings for adding signal (special study)
  //double _peak = 965.4; // best fit combined
  //double _peak = 961.9; // best fit ee
  //double mu_sig_ee   = 5.2;
  //double mu_sig_ee   = 0.0;
  //double mu_sig_mumu = mu_sig_ee * 1.4;
  //std::string sigPdfDimuonName = "sigpdf_dimuon";
  //std::string sigPdfDielectronName = "sigpdf_dielectron";

  // this is the significance we want to evaluate and return
  Double_t significance = 0.0;

  //PrintWsVars(); //FIXME: maybe migrate from twobody.C later
  //Diagnostics(); //FIXME: maybe migrate from twobody.C later
  //mc.Print();


  // get observed numbers of events
  //Int_t n_ee = -1; // safe default
  //Int_t n_mumu = -1;
  //if (do_dimuon) n_mumu = mData["dimuon"]->sumEntries();
  //if (do_dielectron) n_ee = mData["dielectron"]->sumEntries();

  _myWS->var("mass")->setRange(fit_range_low, fit_range_high); //FIXME: no need to have the observable name hardcoded

  Double_t _width = _myWS->var("width_p1")->getVal(); 
  _myWS->var("width_p1")->setVal(_width*width_factor); //COMMENT: not sure what is happening here; probably used to accomodate different signal models

  std::pair<Double_t, Double_t> _max_pllr;
  Double_t _max_pllr_obs = 0.0;
  Double_t _max_pllr_mass = mass_low;

  // let's time this
  TStopwatch t;
  t.Start();

  // observed
  std::string null_parameters_name = "NullParameters";
  std::string null_parameters_snapshot_name = null_parameters_name+ "Snapshot";
  if (mode.find("data")!=std::string::npos){
    RooDataSet* data = _databox->createObservedData();

    // set parameter snapshot for bg-only hypothesis
    RooAbsReal * pNll = _myModelConfig->GetPdf()->createNLL( *data );
    const RooArgSet * pGlobalObs = _myModelConfig->GetGlobalObservables();
    RooArgSet poiAndGlobalObs("poiAndGlobalObs");
    poiAndGlobalObs.add( *pGlobalObs );
    const RooArgSet * pPoi = _myModelConfig->GetParametersOfInterest();
    poiAndGlobalObs.add( *pPoi );
    RooAbsReal * pProfile = pNll->createProfile( poiAndGlobalObs ); // do not profile global observables
    ((RooRealVar *)pPoi->first())->setVal( 0 );  // set xsec=0 here
    pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
    RooArgSet * pPoiAndNuisance = new RooArgSet("poiAndNuisance");
    pPoiAndNuisance->add(*_myModelConfig->GetNuisanceParameters());
    pPoiAndNuisance->add(*_myModelConfig->GetParametersOfInterest());
    _myModelConfig->SetSnapshot(*pPoiAndNuisance);

    // this snapshot is redundant
    // just to save time and not rename anything
    _myWS->defineSet( null_parameters_name.c_str(), *pPoiAndNuisance );
    _myWS->saveSnapshot(null_parameters_snapshot_name.c_str(), *_myWS->set(null_parameters_name.c_str())); //COMMENT: need to clarify the usage of snapshots

    delete pProfile;
    delete pNll;
    delete pPoiAndNuisance;


    // FIXME: testing
    // fit and get signal yield
    //PlotData("dielectron");
    //ws->var("peak")->setVal(_peak);
    //RooFitResult * _fit = ws->pdf("model")->fitTo(*data,Save());
    //ws->pdf("model")->fitTo(*data,Save());
    //ws->function("nsig_dielectron")->Print();
    //ws->function("nsig_dimuon")->Print();

    _max_pllr = get_pllr_max(mass_low, mass_high, mass_step, "nchannels_ratio_best_fit.pdf", null_parameters_name, data);
    _max_pllr_obs = _max_pllr.first;
    _max_pllr_mass = _max_pllr.second;
    
    // print to file
    //std::string _outfile = "dilepton_ratio_pllr_" + mode + suffix + ".ascii";
    std::string _outfile = "nchannels_ratio_pllr_" + mode + suffix + ".ascii";
    printDileptonRatioSignificance( _outfile, _max_pllr_mass, _max_pllr_obs );

    delete data;
  }
  else if (mode.find("pseudo")!=std::string::npos){

    int pe_counter = 0;
    while (pe_counter < ntoys){
      
      std::cout<< legend << funclegend << "this is pseudoexperiment " << pe_counter+1 << " of " << ntoys << std::endl;
      std::cout << legend << funclegend << "for the significance estimate" << std::endl;
      
      // prepare PE data
      RooDataSet * data = _databox->createPeData();

      //std::cout <<"data "<<*data<<std::endl;
      //  data->Print("v");
   // add signal (optional)
   //AddPseudoData("dimuon", sigPdfDimuonName, mu_sig_mumu, _peak, true);
   //AddPseudoData("dielectron", sigPdfDielectronName, mu_sig_ee, _peak, true);
      //std::cout <<"value of nbkg_est "<<_myWS->var("nbkg_est_diele_13TeV_phys14")->getVal()<<std::endl;
      //      _myWS->var("nbkg_est_diele_13TeV_phys14")->setVal(78000);

    // set parameter snapshot for bg-only hypothesis
    RooAbsReal * pNll = _myModelConfig->GetPdf()->createNLL( *data );
    const RooArgSet * pGlobalObs = _myModelConfig->GetGlobalObservables();
    RooArgSet poiAndGlobalObs("poiAndGlobalObs");
    poiAndGlobalObs.add( *pGlobalObs );
    const RooArgSet * pPoi = _myModelConfig->GetParametersOfInterest();
    poiAndGlobalObs.add( *pPoi );
    RooAbsReal * pProfile = pNll->createProfile( poiAndGlobalObs ); // do not profile global observables
    ((RooRealVar *)pPoi->first())->setVal( 0 );  // set xsec=0 here
    pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
    RooArgSet * pPoiAndNuisance = new RooArgSet("poiAndNuisance");
    pPoiAndNuisance->add(*_myModelConfig->GetNuisanceParameters());
    pPoiAndNuisance->add(*_myModelConfig->GetParametersOfInterest());
    _myModelConfig->SetSnapshot(*pPoiAndNuisance);

    // this snapshot is redundant
    // just to save time and not rename anything
    _myWS->defineSet( null_parameters_name.c_str(), *pPoiAndNuisance );
    _myWS->saveSnapshot(null_parameters_snapshot_name.c_str(), *_myWS->set(null_parameters_name.c_str()));

    // _myWS->Print();

    std::cout <<"value of nbkg_est "<<_myWS->var("nbkg_est_diele_13TeV_phys14")->getVal()<<std::endl;

    delete pProfile;
    delete pNll;
    delete pPoiAndNuisance;

      _max_pllr = get_pllr_max(mass_low, mass_high, mass_step, "", null_parameters_name, data);
      _max_pllr_obs = _max_pllr.first;
      _max_pllr_mass = _max_pllr.second;
      
      std::string _outfile = "dilepton_ratio_pllr_" + mode + suffix + ".ascii";
      printDileptonRatioSignificance( _outfile, _max_pllr_mass, _max_pllr_obs );
      
      ++pe_counter;
      delete data;
    } // end of while
  } // end of if (pseudo)

  t.Print();

  return significance;
}

std::pair<Double_t, Double_t> Resultator::get_pllr_max( Double_t mass_low,
                       Double_t mass_high,
                       Double_t mass_step,
                       std::string fit_plot_file_name,
                       std::string null_parameters_name,
                       RooDataSet * data ){
  //
  // Scan the mass range and find the most extreme test statistic value
  // 
  // If fit_plot_file_name is not "", draw fit projections, save to file
  //

  std::string funclegend = " get_pllr_max(...) " ;
  cout << legend << funclegend << " mass_low: " << mass_low << " mass_high: " << mass_high << " mass_step: " << mass_step << std::endl;

  std::pair<Double_t, Double_t> _res;

  //Double_t _max_pllr = 0.0;
  Double_t _max_pllr = -1.0;
  //Double_t _max_mass = 0.0;
  Double_t _max_mass = mass_low;

  // does llr slope with mass? generally, yes but not in this analysis
  // SetPlotStyle(); //FIXME: maybe migrate this later from twobody.C to the Pixie class
  TH1F h_llr("llr","",(int)((mass_high-mass_low)/mass_step)+1,mass_low, mass_high);
  int _co = 1;
  for (Double_t _p = mass_low; _p < mass_high; _p += mass_step){

    // set resonance peak position
    // CONVENTION: variable "peak" for peak position hardcoded
    _myWS->var("peak")->setVal(_p);
    cout << legend << funclegend << " setting peak to: " << _p << endl;

    Double_t _pllr = get_pllr(null_parameters_name, data);
    if (_pllr > _max_pllr){
      _max_pllr = _pllr;
      _max_mass = _p;
    }

    // FIXME: does llr slopes with mass?
    std::cout << legend << funclegend << "DEBUG: pllr=" << _pllr << std::endl;
    h_llr.SetBinContent(_co,_pllr);
    ++_co;
  }
  if (_makeLlrPlot){  //FIXME: migrate this to the new code later 
    //TCanvas c1;
    //h_llr.Draw();
    //CmsPrel();
    //XLabel("dilepton invariant mass");
    //YLabel("Z score");
    //DrawTitle("Local asymptotic Z score");
    //c1.SaveAs("llr.pdf");
  }

  // try adaptive
  _myWS->var("peak")->setVal(_max_mass);
  Double_t _min_range = _myWS->var("peak")->getMin();
  Double_t _max_range = _myWS->var("peak")->getMax();

  // set range for local minimization
  _myWS->var("peak")->setRange(std::max(_max_mass-mass_step,mass_low),
             std::min(_max_mass+mass_step,mass_high));

  _myWS->var("peak")->setConstant(kFALSE);
  //CONVENTION: likelihood of combined model is ahrdcoded as "model"
  RooFitResult * _fit = _myWS->pdf("model")->fitTo(*data,Save());
  //RooFitResult * _fit = 0;
  _max_mass =_myWS->var("peak")->getVal();

  // put all back
  _myWS->var("peak")->setRange(_min_range, _max_range);
  _myWS->var("peak")->setConstant(kTRUE);
  //_max_pllr = get_pllr(null_parameters_name);

  if (fit_plot_file_name.size()>0){
    // make and save fit plot
  /*

    ws->function("nsig_dimuon")->Print();
    ws->function("nsig_dielectron")->Print();
    RooPlot * _frame1 = ws->var("mass")->frame(100) ;
    RooPlot * _frame2 = ws->var("mass")->frame(100) ;
    *data->plotOn(_frame1, Cut("channels==channels::dimuon"), Name("dimuon")) ;
    ws->pdf("model")->plotOn(_frame1,Slice(*ws->cat("channels"),"dimuon"),ProjWData(*ws->cat("channels"),*data)) ;
    *data->plotOn(_frame2, Cut("channels==channels::dielectron"), Name("dielectron")) ;
    ws->pdf("model")->plotOn(_frame2,Slice(*ws->cat("channels"),"dielectron"),ProjWData(*ws->cat("channels"),*data)) ;
    TCanvas c2("c2", "c2", 425, 550);
    c2.Divide(1,2);
    c2.cd(1);
    _frame1->Draw();
    c2.cd(2);
    _frame2->Draw();
    c2.SaveAs(fit_plot_file_name.c_str());
    delete _frame1;
    delete _frame2;
  */
  }

  _max_pllr = get_pllr(null_parameters_name, data);

  _res.first = _max_pllr;
  _res.second = _max_mass;

  // clean up
  //delete _fit;  //COMMENT: deleting _fit crashes the code as destructor is called later from somewhere else in RooFit/RooStats

  return _res;
}

Double_t Resultator::get_pllr( std::string nullParametersName, RooDataSet * data ){
  //
  // Compute the profile log likelihood ratio test statistic
  // Gets the set of null model parameters from the workspace by name
  //

  std::string funclegend = " get_pllr( std::string nullParametersName ) ";
  cout << legend << funclegend << endl;

  // use ProfileLikelihood
  ProfileLikelihoodCalculator plc(*data, *_myModelConfig);
  plc.SetConfidenceLevel(0.68);

  std::string null_parameters_snapshot_name = nullParametersName+ "Snapshot";
  _myWS->loadSnapshot( null_parameters_snapshot_name.c_str() );

  const RooArgSet * pNullParameters = _myWS->set(nullParametersName.c_str());
  plc.SetNullParameters(*pNullParameters);

  // suppress RooFit verbosity
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  RooStats::HypoTestResult * htr = plc.GetHypoTest();
  Double_t significance = htr->Significance();

  std::cout << legend << funclegend << "null p-value = " << htr->NullPValue() << std::endl;
  std::cout << legend << funclegend << "significance = " << significance << std::endl;
  //logfile << legend << funclegend << "significance = " << significance << std::endl; //COMMENT: Why put it into a separate logfile if it goes to cout anyway
  delete htr;

  // return RooFit verbosity to previous level
  RooMsgService::instance().setGlobalKillBelow(msglevel);

  //std::exit(-1);

  return significance;
  //return 4.0;
}

void Resultator::makeMcmcPosteriorPlot( MCMCInterval * mcInt, std::string filename){
   std::string funclegend = " makeMcmcPosteriorPlot( MCMCInterval * mcInt, std::string filename) ";
   cout << legend << funclegend << endl;

   TFile* outfile = new TFile(_plotfile.c_str(),"RECREATE");
   //outfile->cd();
   TH1D* myhist = new TH1D("testhist","testhist",200, 0.9, 1.1); //histogram for storing information from special tests

   std::string p0name = "width_p0"; //Comment: unfortunately hardcoded
   //p0name += "_";
   //p0name += *vecIt;
   std::string p1name = "width_p1"; //Comment: unfortunately hardcoded

   cout << p0name << " is: " << _myWS->var(p0name.c_str())->getVal() << endl;
   cout << p1name << " is: " << _myWS->var(p1name.c_str())->getVal() << endl;

   cout << "Hei-ho 1" << endl;
//    //make plot of RooVoigtian
//    RooVoigtian * myVoigt = static_cast<RooVoigtian *>(_myWS->pdf("sigpdf_dimuon2012"));
//    RooRealVar * mymass = static_cast<RooRealVar *>(_myWS->var("mass"));
//    RooArgSet * mySet = new RooArgSet(*mymass);
//    RooArgSet * mySet2 = new RooArgSet(*mymass);
//    RooRealVar * my_beta_nsig_dimuon2012 = static_cast<RooRealVar *>(_myWS->var("beta_nsig_dimuon2012"));
//    RooRealVar * my_beta_nbkg_dimuon2012 = static_cast<RooRealVar *>(_myWS->var("beta_nbkg_dimuon2012"));
//    mySet2->add(*my_beta_nsig_dimuon2012);
//    mySet2->add(*my_beta_nbkg_dimuon2012);
//    cout << "signal shape is normalized to: " << myVoigt->getNorm(mySet) << endl;
// 
//    RooRealVar mass ("mass_scaled_dimuon2012","mass_scaled_dimuon2012",200.,8000.);
//    //RooRealVar testvar = = static_cast<RooVoigtian *>(_myWS->var("mass"));
// 
//    RooPlot* frame = mass.frame();
// 
//    myVoigt->plotOn(frame);
// 
//    outfile->WriteTObject(frame);
//    cout << "Hei-ho 6" << endl;
// 
//    RooSimultaneous * mySimultaneous_model  = static_cast<RooSimultaneous*>(_myWS->pdf("model"));
//    cout << "model is normalized to: " << mySimultaneous_model->getNorm(mySet2) << endl;
//    RooSimultaneous * mySimultaneous_model_dimuon2012  = static_cast<RooSimultaneous*>(_myWS->pdf("model_dimuon2012"));
//    cout << "model is normalized to: " << mySimultaneous_model_dimuon2012->getNorm(mySet2) << endl;

   //test likelihood evaluation

   // load and print S+B Model Config
   RooStats::ModelConfig * pSbHypo = _myModelConfig;
   pSbHypo->Print();
   
   // get RooArgSet with poi 
   RooArgSet * poiSet = new RooArgSet();
   poiSet->add(*(_myWS->var("ratio")));
   
   // set parameter snapshot that corresponds to the best fit to data
   RooAbsReal * pNll = pSbHypo->GetPdf()->createNLL( *(_databox->createObservedData()) );
   
   // global fit preparation
   // (RooArgSet here defines the parameters that are NOT fitted)
   //RooAbsReal * pProfile = pNll->createProfile( RooArgSet() ); 
   
   // ...or do not profile global observables
   RooAbsReal * pProfile = pNll->createProfile( *poiSet ); 
   
   // do fit and set POI and nuisance parameters to fitted values
   double _pnll = pProfile->getVal(); 
   
   std::cout << std::endl;
   std::cout << "profiled: -log L = " << _pnll << std::endl;


   //frame->Write();

   if (mcInt && filename.size() > 0){

      TCanvas c1("c1");
      MCMCIntervalPlot plot(*mcInt);
      plot.Draw();
      std::string plot_name = "posterior_" + filename;
      c1.SetName(plot_name.c_str());
      outfile->WriteTObject(&c1);

      TCanvas c2("c2");
      plot.DrawNLLHist();
      plot_name = "nll_" + filename;
      c2.SetName(plot_name.c_str());
      outfile->WriteTObject(&c2);

      std::vector<std::string> _vNames;
      _vNames.push_back("ratio");
      _vNames.push_back("beta_nsig_dimuon2011");
      _vNames.push_back("beta_nbkg_dimuon2011");
      _vNames.push_back("beta_nsig_dielectron2011");
      _vNames.push_back("beta_nbkg_dielectron2011");
      _vNames.push_back("beta_mass_dielectron2011");
      _vNames.push_back("beta_nsig_dimuon2012");
      _vNames.push_back("beta_nbkg_dimuon2012");
      _vNames.push_back("beta_nsig_dielectron2012");
      _vNames.push_back("beta_nbkg_dielectron2012");
      _vNames.push_back("beta_mass_dielectron2012");
      
      plot_name = filename;

      for (std::vector<std::string>::const_iterator _name = _vNames.begin(); _name != _vNames.end(); ++_name){

         RooRealVar * _nuis = _myWS->var( _name->c_str() );
         if (_nuis){
            TCanvas c3("c3");
            plot.DrawChainScatter(*_myWS->var("ratio"), *_nuis);
            std::string _fname( _name->c_str() );
            _fname += "_vs_poi_"+plot_name;
            c3.SetName(_fname.c_str());
            outfile->WriteTObject(&c3);
         }
      }

      // test scatter
      //TCanvas c2("c2");
      //plot.DrawChainScatter(*ws->var("ratio"),*ws->var("beta_nsig_dielectron"));
      //c2.SaveAs("scatter_mcmc_e.png");
      // 
      //TCanvas c3("c3");
      //plot.DrawChainScatter(*ws->var("ratio"),*ws->var("beta_nsig_dimuon"));
      //c3.SaveAs("scatter_mcmc_mu.png");
   }
   outfile->Close();

   bool debugflag = false;
   if(debugflag){
      //dump Markov Chain
      cout << "start: dumping Markov Chain ..." << endl;
      cout << "mcInt->GetDimension() : " << mcInt->GetDimension() << endl;
      const MarkovChain * mychain = mcInt->GetChain();
      int chainlength = mychain->Size();
      cout.precision(6);
      for(int i = 1; i<chainlength; i++){
            cout << "element: " << i << " NLL: " << mychain->NLL(i) << mychain->Weight(i) << endl;
            const RooArgSet* myset = mychain->Get();
            myset->writeToStream(cout, false);
            double mass_scaled = TMath::Power(1.01,myset->getRealValue("beta_mass_dielectron2011"))*450.;
            cout << "mass_scaled:" << mass_scaled << endl;
      }
      cout << " ... end: dumping Markov Chain" << endl;
   }
   return;
}




