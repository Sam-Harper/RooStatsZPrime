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



using namespace RooFit;
using namespace std;


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

   bool run_channel_dielectron_ebeb = true;
   bool run_channel_dielectron_ebee = true;
   bool run_channel_dimuon = true;
   bool writeplots = true;
   bool datapruning = true;
   bool usemassscaleuncer = true;

   double Zprimemass = static_cast<double>(atof(argv[4]));
   unsigned int ntoys  = atoi(argv[5]);
   unsigned int MCMCiter  = atoi(argv[6]);

   string filesuffix = argv[7];
   string mode = argv[8];
   string workspacedir = argv[9];

   int help_bool = atoi(argv[10]);
   if(help_bool == 0){
      writeplots = false;
   }
   std::string plotfile = "";
   if(writeplots){
      plotfile = plotfile + "plots_" + filesuffix + ".root";
   }

   help_bool = atoi(argv[11]);
   if(help_bool == 0){
      datapruning = false;
   }

   double poiRangeFactor = static_cast<double>(atof(argv[12]));

   help_bool = atoi(argv[13]);
   if(help_bool == 0){
      usemassscaleuncer = false;
   }

   int help_run_channel = atoi(argv[1]);
   if(help_run_channel == 0){
      run_channel_dielectron_ebeb = false;
   }
   else{
      run_channel_dielectron_ebeb = true;
   }
   help_run_channel = atoi(argv[2]);
   if(help_run_channel == 0){
      run_channel_dielectron_ebee = false;
   }
   else{
      run_channel_dielectron_ebee = true;
   }
   help_run_channel = atoi(argv[3]);
   if(help_run_channel == 0){
      run_channel_dimuon = false;
   }
   else{
      run_channel_dimuon = true;
   }


cout << "calculate some limits ..." << endl;

std::map<string, RooWorkspace*> WSvec; // maps channel names to workspaces

// DEFINE WOKSPACES AND CHANNEL NAMES

// CHANNEL: dielectron EBEB

const std::string filename1 = workspacedir + "/ws_dielectron_ebeb_ratio_Moriond13_v2.root"; // the root file containing the workspace FIXME:update
const std::string ws_name1 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname1 = "dielectron_ebeb"; // name of the channel -> to be used in your daugther class of ModelConfigurator

if(run_channel_dielectron_ebeb){

cout << "prepare with : run_channel_dielectron_ebeb" << endl;

TFile file1(filename1.c_str(), "read"); // construct TFile object to load the workspace
RooWorkspace * ws1 = (RooWorkspace *)file1.Get( ws_name1.c_str() ); // load the workspace
cout << "workspacename: " << ws1->GetName();


WSvec.insert( pair<string, RooWorkspace *>(channelname1, ws1) ); // insert channel info

}

// CHANNEL: dielectron EBEE

const std::string filename2 = workspacedir + "/ws_dielectron_ebee_ratio_Moriond13_v2.root"; // the root file containing the workspace FIXME:update
const std::string ws_name2 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname2 = "dielectron_ebee"; // name of the channel -> to be used in your daugther class of ModelConfigurator

if(run_channel_dielectron_ebee){

cout << "prepare with : run_channel_dielectron_ebee" << endl;

TFile file2(filename2.c_str(), "read"); // construct TFile object to load the workspace

RooWorkspace * ws2 = (RooWorkspace *)file2.Get( ws_name2.c_str() ); // load the workspace
cout << "workspacename: " << ws2->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname2, ws2) ); // insert channel info

}

// CHANNEL: dimuon

const std::string filename3 = workspacedir + "/ws_dimuon_ratio_Moriond13_v1.root"; // the root file containing the workspace FIXME:update
const std::string ws_name3 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname3 = "dimuon"; // name of the channel -> to be used in your daugther class of ModelConfigurator

if(run_channel_dimuon){
TFile file3(filename3.c_str(), "read"); // construct TFile object to load the workspace

RooWorkspace * ws3 = (RooWorkspace *)file3.Get( ws_name3.c_str() ); // load the workspace
cout << "workspacename: " << ws3->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname3, ws3) ); // insert channel info

}


// CREATE MODELCONFIGURATOR

ModelConfiguratorZprime * myConfigurator = new ModelConfiguratorZprime(WSvec); // construct ModelConfigurator

// ADDITIONAL INFORMATION FOR COMBINING THE CHANNELS

//define list with shared Vars
std::string sharedVarsString = "peak,mass,ratio,width_p0,width_p1,width"; //FIXME:check
// set shared variables 
// COMMENT: this structure may not be sufficient
myConfigurator->setSharedVars(sharedVarsString);
// COMMENT: for combination of 2011 and 2012 data we will also need to handle correlated Variables -> think about best way to do this

// define the parameter of interest
// CONVENTION: poi needs to have the same name in all single channel workspaces
std::string poiString = "ratio";
myConfigurator->setPoiString(poiString);


// ASSOCIATE VARIABLES WITH STATISTICAL MEANING

// CHANNEL: dielectron_ebeb

if(run_channel_dielectron_ebeb){

string nuisanceParamsString_channel1 = "";
if (usemassscaleuncer){
nuisanceParamsString_channel1 = "beta_nsig,beta_nbkg,beta_mass"; //FIXME:check
}
else{
nuisanceParamsString_channel1 = "beta_nsig,beta_nbkg"; //FIXME:check
}
string GlobalObsParamsString_channel1 = "glob_nsig,glob_nbkg,glob_mass"; //FIXME:check
string ObservableParamsString_channel1 = "mass";

myConfigurator->setNuisanceParams(channelname1, nuisanceParamsString_channel1);
myConfigurator->setGlobalObs(channelname1, GlobalObsParamsString_channel1);
myConfigurator->setObservables(channelname1, ObservableParamsString_channel1);

}

// CHANNEL: dielectron_ebee

if(run_channel_dielectron_ebee){

string nuisanceParamsString_channel2 = "";
if (usemassscaleuncer){
nuisanceParamsString_channel2 = "beta_nsig,beta_nbkg,beta_mass"; //FIXME:check
}
else{
nuisanceParamsString_channel2 = "beta_nsig,beta_nbkg"; //FIXME:check
}
string GlobalObsParamsString_channel2 = "glob_nsig,glob_nbkg,glob_mass"; //FIXME:check
string ObservableParamsString_channel2 = "mass";

myConfigurator->setNuisanceParams(channelname2, nuisanceParamsString_channel2);
myConfigurator->setGlobalObs(channelname2, GlobalObsParamsString_channel2);
myConfigurator->setObservables(channelname2, ObservableParamsString_channel2);

}

// CHANNEL: dimuon

if(run_channel_dimuon){

string nuisanceParamsString_channel3 = "beta_nsig,beta_nbkg"; //FIXME:check
string GlobalObsParamsString_channel3 = "glob_nsig,glob_nbkg"; //FIXME:check
string ObservableParamsString_channel3 = "mass";

myConfigurator->setNuisanceParams(channelname3, nuisanceParamsString_channel3);
myConfigurator->setGlobalObs(channelname3, GlobalObsParamsString_channel3);
myConfigurator->setObservables(channelname3, ObservableParamsString_channel3);

}



//setup combined RooWorkSpace and ModelConfig
myConfigurator->Setup();

// FIX MASS HYPOTHESIS
myConfigurator->setMassHypothesis(Zprimemass);

//safe combined WS

//myConfigurator->WriteCombinedWS("CombinedWS.root");

// ADJUST WORKSPACE FOR 7/8 TeV Combination

Pixie * myPixie = new Pixie();


//Setup DataPruner for customized lower mass thresholds
// std::map<std::string , double> Rangemap;
// 
// if(datapruning){
//    if(run_channel1){ Rangemap.insert( pair<std::string, double>("dimuon2011", 200.) );}
//    if(run_channel2){ Rangemap.insert( pair<std::string, double>("dielectron2011", 200.) );}
//    if(run_channel3){ Rangemap.insert( pair<std::string, double>("dimuon2012", 200.) );}
//    if(run_channel4){ Rangemap.insert( pair<std::string, double>("dielectron2012", 200.) );}
// }
// 
// DataPruner * mydatapruner = new DataPruner(Rangemap);
// 
// if(datapruning){
//    //For testing: adjust some other parameters: nbkg_dielectron2011, nbkg_dimuon2011
//    //myConfigurator->setVar("peak",1000);
//    //myConfigurator->setVar("nbkg_est_dielectron2011",458);
//    //myConfigurator->setVar("nbkg_est_dimuon2011",570);
//    //myConfigurator->setVarRange("mass", 350., 2000.);
//    myConfigurator->setVarRange("mass", 200., 3000.);
// }

//Setup DataPruner based on signal width
DataPruner * mydatapruner;
if(datapruning){
//CONVEBTION: the mass thresholds for at least 400 events rae currently hardcoded -> move to the single channel cfg files
   double thresh_min_ee_ebeb = 500.;
   double thresh_min_ee_ebee = 500.;
   double thresh_min_mm = 550.;
   double thresh_min = 20000. ;
   if(run_channel_dielectron_ebeb && run_channel_dielectron_ebee && run_channel_dimuon ){ thresh_min = std::min(thresh_min_ee_ebeb, std::min(thresh_min_ee_ebee, thresh_min_mm) ) ;}
   if(run_channel_dielectron_ebeb && run_channel_dielectron_ebee && !(run_channel_dimuon) ){ thresh_min = std::min(thresh_min_ee_ebeb, thresh_min_ee_ebee) ;}
   if(run_channel_dielectron_ebeb && !(run_channel_dielectron_ebee) && run_channel_dimuon ){ thresh_min = std::min(thresh_min_ee_ebeb, thresh_min_mm) ;}
   if( !(run_channel_dielectron_ebeb) && run_channel_dielectron_ebee && run_channel_dimuon ){ thresh_min = std::min(thresh_min_ee_ebee, thresh_min_mm) ;}
   if( run_channel_dielectron_ebeb && !(run_channel_dielectron_ebee) && !(run_channel_dimuon) ){ thresh_min = thresh_min_ee_ebeb ;}
   if( !(run_channel_dielectron_ebeb) && run_channel_dielectron_ebee && !(run_channel_dimuon) ){ thresh_min = thresh_min_ee_ebee ;}
   if( !(run_channel_dielectron_ebeb) && !(run_channel_dielectron_ebee) && run_channel_dimuon){ thresh_min = thresh_min_mm ;}
   cout << "thresh_min for data pruning is : " << thresh_min << endl;
   mydatapruner = new DataPruner(6.0, myConfigurator, thresh_min, -1);
}
else{
   //default initialization that does not do any pruning
   std::map<std::string , double> Rangemap;
   mydatapruner = new DataPruner(Rangemap, myConfigurator);
}
//some settings for tests
// myConfigurator->setVarRange("beta_nsig_dielectron2012", -3.5, 3.5);
// myConfigurator->setVarRange("beta_nbkg_dielectron2012", -3.5, 3.5);
// myConfigurator->setVarRange("beta_mass_dielectron2012", -3.5, 3.5);
// myConfigurator->setVarRange("beta_nsig_dielectron2011", -3.5, 3.5);
// myConfigurator->setVarRange("beta_nbkg_dielectron2011", -3.5, 3.5);
// myConfigurator->setVarRange("beta_mass_dielectron2011", -0.67, -0.63);
// myConfigurator->setVar("beta_mass_dielectron2011",-0.65);
// myConfigurator->setVarRange("glob_mass_dielectron2011", -0.0001, 0.0001);
// myConfigurator->setVar("glob_mass_dielectron2011",0.0);
// myConfigurator->setVarRange("glob_mass_dielectron2012", -0.0001, 0.0001);
// myConfigurator->setVar("glob_mass_dielectron2012",0.0);

//adjust signal resolution for the ZPSI model

myPixie->SetResolution_ZPSI(myConfigurator->getCombinedWS(), myConfigurator->getChannelNames());


//just removing the uncertainties from the list of nuisance parameters does not switch off their variation in the Markov Chain (probaly because the prior term is included as a part of the likelihood)
if (!usemassscaleuncer){
   cout << "mass scale uncertainties are not applied!" << endl;
   if(run_channel_dielectron_ebeb){
      myConfigurator->setVarRange("beta_mass_dielectron_ebeb", -0.001, 0.001);
      myConfigurator->setVar("beta_mass_dielectron_ebeb", 0.0);
   }
   if(run_channel_dielectron_ebee){
      myConfigurator->setVarRange("beta_mass_dielectron_ebee", -0.001, 0.001);
      myConfigurator->setVar("beta_mass_dielectron_ebee", 0.0);
   }
}

// RUN Bayesian limits

Resultator * myResultator = new Resultator(myConfigurator, mydatapruner, plotfile);
myResultator->setNbinsPosterior(300);

//Estimate reasonable POI range
PoiRangeEstimator * myPoiRangeEstimator = new PoiRangeEstimator(myConfigurator, myResultator);
double poiUpperLimitGuess = myPoiRangeEstimator->GetPoiUpper("multi", Zprimemass);
cout << "estimate for reasonable upper limit of poi range: " << poiUpperLimitGuess << endl; 
delete myPoiRangeEstimator;
myResultator->SetPoiUpperLimitByHand(poiUpperLimitGuess*poiRangeFactor);
myResultator->setPoiRangeFactor( poiRangeFactor);
//myResultator->SetPoiUpperLimitByHand(0.09); //this line is just for testing
myResultator->calculateMCMClimit( MCMCiter, 100, ntoys, filesuffix, mode);

cout << ".. calculated some limits" << endl;

}
