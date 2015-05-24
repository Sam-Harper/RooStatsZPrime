#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "TSystem.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "ModelConfiguratorZprime.hh"
#include "Resultator.hh"
#include "PoiRangeEstimator.hh"
#include "DataPruner.hh"
#include "Pixie.hh"
//#include <libconfig.hh>

#include "RooRealVar.h"
#include "RooNumIntConfig.h"

#include "CmdLineInt.hh"
#include <boost/algorithm/string.hpp>

void getWorkSpaces(const std::string& filename,const std::string& chanNamesStr,std::map<std::string,RooWorkspace*>& workSpaces);
void setupChanParams(const std::map<std::string,RooWorkspace*>& workSpaces,ModelConfiguratorZprime* modelConfig);

using namespace RooFit;


int main(int argc, char* argv[]) {
  gSystem->AddIncludePath("-I$ROOFITSYS/include");
   //without this settings root is having troubles to correctly normalize the model
   RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
   RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","61Points") ;
   RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",1000) ;
   RooAbsPdf::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
   RooAbsPdf::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","61Points") ;
   RooAbsPdf::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",1000) ;

   //alternative without changing the default algorithmF for integration
   //RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1.E-14);
   //RooAbsReal::defaultIntegratorConfig()->setEpsRel(1.E-14);

   //read options

   char inputFilename[256];
   char outFilename[256];
   int nrJobs,jobNr;
   char chanNames[256];
   int nrExpts,nrMCMCIt;
   double resMass;
   bool pruneData;
   char mode[256];

   CmdLineInt cmdLineInt(argv[0]);
   cmdLineInt.addNonOption(inputFilename,true," ","input filename");
   cmdLineInt.addNonOption(outFilename,true," ","output filename");
   cmdLineInt.addOption("n",&nrJobs,1,"Nr of Jobs to split input into");
   cmdLineInt.addOption("j",&jobNr,1,"Job Nr");
   cmdLineInt.addOption("chans",chanNames,"electron_eb_phys14","channels (seperated by a \":\") to run");
   cmdLineInt.addOption("nrExpts",&nrExpts,100,"nr of pseudo experiments to run");
   cmdLineInt.addOption("nrMCMCIt",&nrMCMCIt,10000,"nr of itterations used for the MCMC Integration (10000 is a good number)");
   cmdLineInt.addOption("resMass",&resMass,1000," new resonance mass");
   cmdLineInt.addOption("pruneData",&pruneData,false," trim the data to a window");
   cmdLineInt.addOption("mode",mode,"pseudo","mode = data or pseudo");
   if(!cmdLineInt.processCmdLine(argc,argv)) return 0; //exit if we havnt managed to get required parameters
   
 

   std::cout << "calculate some limits ..." <<std::endl;

   std::map<std::string, RooWorkspace*> workSpaces; // maps channel names to workspaces
   getWorkSpaces(inputFilename,chanNames,workSpaces);

   std::cout <<"configuring model"<<std::endl;
   ModelConfiguratorZprime * myConfigurator = new ModelConfiguratorZprime(workSpaces); 
   std::cout <<"configured model"<<std::endl;

   //define list with shared Vars
   std::string sharedVarsString = "peak,mass,ratio,width_p0,width_p1,width"; //FIXME:check
   myConfigurator->setSharedVars(sharedVarsString);

   // CONVENTION: poi needs to have the same name in all single channel workspaces
   std::string poiString = "ratio";
   myConfigurator->setPoiString(poiString);
   
   setupChanParams(workSpaces,myConfigurator); //sorts out nusiance and obs params
   
   //setup combined RooWorkSpace and ModelConfig
   myConfigurator->Setup();
   
   // FIX MASS HYPOTHESIS
   myConfigurator->setMassHypothesis(resMass);
   
   //safe combined WS

   //myConfigurator->WriteCombinedWS("CombinedWS.root");
   
   // ADJUST WORKSPACE FOR 7/8 TeV Combination
   
   Pixie * myPixie = new Pixie();
   myPixie->SetResolution_ZPSI(myConfigurator->getCombinedWS(), myConfigurator->getChannelNames());

   
//    //just removing the uncertainties from the list of nuisance parameters does not switch off their variation in the Markov Chain (probaly because the prior term is included as a part of the likelihood)
//    if (!usemassscaleuncer){
//      cout << "mass scale uncertainties are not applied!" << endl;
//      if(run_channel_dielectron_ebeb){
//        myConfigurator->setVarRange("beta_mass_dielectron_ebeb", -0.001, 0.001);
//        myConfigurator->setVar("beta_mass_dielectron_ebeb", 0.0);
//    }
//      if(run_channel_dielectron_ebee){
//        myConfigurator->setVarRange("beta_mass_dielectron_ebee", -0.001, 0.001);
//       myConfigurator->setVar("beta_mass_dielectron_ebee", 0.0);
//      }
//    }
   
   
   
   cout << ".. calculated some limits" << endl;
   
   
   //Setup DataPruner
   std::map<std::string , double> rangeMap;
   if(pruneData){
     for(auto& chan : workSpaces) rangeMap.insert(std::pair<std::string,double>(chan.first,400));
   }
   DataPruner* myDataPruner = new DataPruner(rangeMap,myConfigurator);
   
   
   double mass_low = 600;
   double mass_high = 2000;
   double mass_step = 10;
   double fit_range_low = 400;
   double fit_range_high = 2000;
   double width_factor = 1;
   
   // RUN Significance results
   Resultator * myResultator = new Resultator(myConfigurator,myDataPruner);
   myResultator->calculateRatioSignificance( mode, mass_low, mass_high, mass_step, nrExpts, outFilename, fit_range_low, fit_range_high, width_factor, false) ;
   
   cout << ".. did some testing" << endl;
   
   
   
}


void getWorkSpaces(const std::string& filename,const std::string& chanNamesStr,std::map<std::string,RooWorkspace*>& workSpaces)
{
  TFile file(filename.c_str(),"READ");

  std::vector<std::string> chanNames;
  boost::split(chanNames,chanNamesStr,boost::is_any_of(":"));
  for(auto& chanName : chanNames) {
    auto ws = static_cast<RooWorkspace*>(file.Get(chanName.c_str()));
    if(ws){
      workSpaces.insert(std::pair<std::string,RooWorkspace*>(chanName,ws));
      std::cout <<"added channel "<<chanName<<std::endl;
    }else{
      std::cout <<"channel "<<chanName<<" not found, exiting "<<std::endl;
      exit(0);
    }
  }

}

void setupChanParams(const std::map<std::string,RooWorkspace*>& workSpaces,ModelConfiguratorZprime* modelConfig)
{
  const std::string nuisanceParams("beta_nsig,beta_nbkg");
  const std::string globalObsParams("glob_nsig,glob_nbkg,glob_mass");
  const std::string obsParams("mass");

  for(auto& chan : workSpaces){
    modelConfig->setNuisanceParams(chan.first,nuisanceParams);
    modelConfig->setGlobalObs(chan.first,globalObsParams);
    modelConfig->setObservables(chan.first,obsParams);
  }

}
