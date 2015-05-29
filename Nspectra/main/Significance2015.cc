#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "TSystem.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "RooStatsZPrime/ModelConfiguratorZprime.hh"
#include "RooStatsZPrime/Resultator.hh"
#include "RooStatsZPrime/PoiRangeEstimator.hh"
#include "RooStatsZPrime/DataPruner.hh"
#include "RooStatsZPrime/Pixie.hh"
//#include <libconfig.hh>

#include "RooRealVar.h"
#include "RooNumIntConfig.h"

#include "Utility/CmdLineInt.hh"
#include <boost/algorithm/string.hpp>

void getWorkSpaces(const std::string& filename,const std::string& chanNamesStr,std::map<std::string,RooWorkspace*>& workSpaces);
void setupChanParams(const std::map<std::string,RooWorkspace*>& workSpaces,ModelConfiguratorZprime* modelConfig);
void rmMassScaleUncertainty(const std::map<std::string,RooWorkspace*>& workSpaces,ModelConfiguratorZprime* modelConfig);

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
   double minMass,maxMass,massStep;
   double fitRangeLow,fitRangeHigh;

   CmdLineInt cmdLineInt(argv[0]);
   cmdLineInt.addNonOption(inputFilename,true," ","input filename");
   cmdLineInt.addNonOption(outFilename,true," ","output filename");
   cmdLineInt.addOption("n",&nrJobs,1,"Nr of Jobs to split input into");
   cmdLineInt.addOption("j",&jobNr,1,"Job Nr");
   cmdLineInt.addOption("chans",chanNames,"electron_eb_phys14","channels (seperated by a \":\") to run");
   cmdLineInt.addOption("nrExpts",&nrExpts,100,"nr of pseudo experiments to run");
   cmdLineInt.addOption("pruneData",&pruneData,false," trim the data to a window");
   cmdLineInt.addOption("mode",mode,"pseudo","mode = data or pseudo");
   cmdLineInt.addOption("nrMCMCIt",&nrMCMCIt,10000,"nr of itterations used for the MCMC Integration for limit");
   cmdLineInt.addOption("resMass",&resMass,1000," new resonance mass for limit");
   cmdLineInt.addOption("minMass",&minMass,600,"min mass for sig scan");
   cmdLineInt.addOption("maxMass",&maxMass,2000,"max mass for sig scan");
   cmdLineInt.addOption("massStep",&massStep,10,"mass step for sig scan");
   cmdLineInt.addOption("fitRangeLow",&fitRangeLow,200,"min mass for fit range for sig scan");
   cmdLineInt.addOption("fitRangeHigh",&fitRangeHigh,4000,"max mass for fit range for sig scan");
   if(!cmdLineInt.processCmdLine(argc,argv)) return 0; //exit if we havnt managed to get required parameters
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   std::map<std::string, RooWorkspace*> workSpaces; // maps channel names to workspaces
   getWorkSpaces(inputFilename,chanNames,workSpaces);

   ModelConfiguratorZprime * myConfigurator = new ModelConfiguratorZprime(workSpaces); 


   //define list with shared Vars
   std::string sharedVarsString = "peak,mass,ratio,width_p0,width_p1,width"; //FIXME:check
   myConfigurator->setSharedVars(sharedVarsString);
   
   // CONVENTION: poi needs to have the same name in all single channel workspaces
   std::string poiString = "ratio";
   myConfigurator->setPoiString(poiString);
     
   setupChanParams(workSpaces,myConfigurator); //sorts out nusiance and obs params
   myConfigurator->Setup();
   rmMassScaleUncertainty(workSpaces,myConfigurator); 
   
 
   //Setup DataPruner
   std::map<std::string , double> rangeMap;
   if(pruneData){
     for(auto& chan : workSpaces) rangeMap.insert(std::pair<std::string,double>(chan.first,400));
   }
   DataPruner* myDataPruner = new DataPruner(rangeMap,myConfigurator);
   
   // RUN Significance results
   Resultator * myResultator = new Resultator(myConfigurator,myDataPruner);
   myResultator->calculateRatioSignificance( mode, minMass, maxMass, massStep, nrExpts, outFilename, fitRangeLow, fitRangeHigh, 1, false) ;
   
   return 0;
   
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

void rmMassScaleUncertainty(const std::map<std::string,RooWorkspace*>& workSpaces,ModelConfiguratorZprime* modelConfig)
{

  for(auto& chan : workSpaces){
    if(modelConfig->getCombinedWS()->var(("beta_mass_"+chan.first).c_str())){
      modelConfig->setVarRange("beta_mass_"+chan.first,-0.001,0.001);
      modelConfig->setVar("beta_mass_"+chan.first,0.0);
    }
  }

}
