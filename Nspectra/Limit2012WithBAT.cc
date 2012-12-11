#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "RooWorkspace.h"
#include "TFile.h"

#include "ModelConfiguratorZprime.hh"
#include "ResultatorWithBAT.hh"
#include "PoiRangeEstimator.hh"
#include "DataPruner.hh"
#include "Pixie.hh"
//#include <libconfig.hh>

#include "RooRealVar.h"
#include "RooNumIntConfig.h"

using namespace RooFit;
using namespace std;


int main(int argc, char* argv[]) {

   RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");

   //read options

   bool run_channel1 = true;
   bool run_channel2 = true;
   bool run_channel3 = true;
   bool run_channel4 = true;
   bool writeplots = true;
   bool datapruning = true;
   bool usemassscaleuncer = true;

   double Zprimemass = static_cast<double>(atof(argv[5]));
   unsigned int ntoys  = atoi(argv[6]);
   unsigned int MCMCiter  = atoi(argv[7]);

   string filesuffix = argv[8];
   string mode = argv[9];
   string workspacedir = argv[10];

   int help_bool = atoi(argv[11]);
   if(help_bool == 0){
      writeplots = false;
   }
   std::string plotfile = "";
   if(writeplots){
      plotfile = plotfile + "plots_" + filesuffix + ".root";
   }

   help_bool = atoi(argv[12]);
   if(help_bool == 0){
      datapruning = false;
   }

   double poiRangeFactor = static_cast<double>(atof(argv[13]));

   help_bool = atoi(argv[14]);
   if(help_bool == 0){
      usemassscaleuncer = false;
   }

   int help_run_channel = atoi(argv[1]);
   if(help_run_channel == 0){
      run_channel1 = false;
   }
   else{
      run_channel1 = true;
   }
   help_run_channel = atoi(argv[2]);
   if(help_run_channel == 0){
      run_channel2 = false;
   }
   else{
      run_channel2 = true;
   }
   help_run_channel = atoi(argv[3]);
   if(help_run_channel == 0){
      run_channel3 = false;
   }
   else{
      run_channel3 = true;
   }
   help_run_channel = atoi(argv[4]);
   if(help_run_channel == 0){
      run_channel4 = false;
   }
   else{
      run_channel4 = true;
   }



cout << "calculate some limits ..." << endl;

std::map<string, RooWorkspace*> WSvec; // maps channel names to workspaces

// DEFINE WOKSPACES AND CHANNEL NAMES

// CHANNEL: dimuon2011

const std::string filename1 = workspacedir + "/ws_dimuon_ratio_full2011v1.root"; // the root file containing the workspace
const std::string ws_name1 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname1 = "dimuon2011"; // name of the channel -> to be used in your daugther class of ModelConfigurator
TFile file1(filename1.c_str(), "read"); // construct TFile object to load the workspace

if(run_channel1){

RooWorkspace * ws1 = (RooWorkspace *)file1.Get( ws_name1.c_str() ); // load the workspace
cout << "workspacename: " << ws1->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname1, ws1) ); // insert channel info

}

// CHANNEL: dimuon2012

const std::string filename3 = workspacedir + "/ws_dimuon_ratio_prelim_2012.root"; // the root file containing the workspace
const std::string ws_name3 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname3 = "dimuon2012"; // name of the channel -> to be used in your daugther class of ModelConfigurator
TFile file3(filename3.c_str(), "read"); // construct TFile object to load the workspace

if(run_channel3){

RooWorkspace * ws3 = (RooWorkspace *)file3.Get( ws_name3.c_str() ); // load the workspace
cout << "workspacename: " << ws3->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname3, ws3) ); // insert channel info

}

// CHANNEL: dielectron2011

const std::string filename2 = workspacedir + "/ws_dielectron_ratio_full2011v1.root"; // the root file containing the workspace
const std::string ws_name2 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname2 = "dielectron2011"; // name of the channel -> to be used in your daugther class of ModelConfigurator
TFile file2(filename2.c_str(), "read"); // construct TFile object to load the workspace

if(run_channel2){

RooWorkspace * ws2 = (RooWorkspace *)file2.Get( ws_name2.c_str() ); // load the workspace
cout << "workspacename: " << ws2->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname2, ws2) ); // insert channel info

}

// CHANNEL: dielectron2012

const std::string filename4 = workspacedir + "/ws_dielectron_ratio_prelim_2012.root"; // the root file containing the workspace
const std::string ws_name4 = "myWS"; // the name of the workspace TObject to be used
const std::string channelname4 = "dielectron2012"; // name of the channel -> to be used in your daugther class of ModelConfigurator
TFile file4(filename4.c_str(), "read"); // construct TFile object to load the workspace

if(run_channel4){

RooWorkspace * ws4 = (RooWorkspace *)file4.Get( ws_name4.c_str() ); // load the workspace
cout << "workspacename: " << ws4->GetName();

WSvec.insert( pair<string, RooWorkspace *>(channelname4, ws4) ); // insert channel info

}

// CREATE MODELCONFIGURATOR

ModelConfiguratorZprime * myConfigurator = new ModelConfiguratorZprime(WSvec); // construct ModelConfigurator

// ADDITIONAL INFORMATION FOR COMBINING THE CHANNELS

//define list with shared Vars
std::string sharedVarsString = "peak,mass,ratio,width_p0,width_p1,width";
// set shared variables 
// COMMENT: this structure may not be sufficient
myConfigurator->setSharedVars(sharedVarsString);
// COMMENT: for combination of 2011 and 2012 data we will also need to handle correlated Variables -> think about best way to do this

// define the parameter of interest
// CONVENTION: poi needs to have the same name in all single channel workspaces
std::string poiString = "ratio";
myConfigurator->setPoiString(poiString);


// ASSOCIATE VARIABLES WITH STATISTICAL MEANING

// CHANNEL: dimuon2011

if(run_channel1){

//see twobody.C line 1790 ff
string nuisanceParamsString_channel1 = "beta_nsig,beta_nbkg";
string GlobalObsParamsString_channel1 = "glob_nsig,glob_nbkg";
string ObservableParamsString_channel1 = "mass";

myConfigurator->setNuisanceParams(channelname1, nuisanceParamsString_channel1);
myConfigurator->setGlobalObs(channelname1, GlobalObsParamsString_channel1);
myConfigurator->setObservables(channelname1, ObservableParamsString_channel1);

}

// CHANNEL: dimuon2012

if(run_channel3){

string nuisanceParamsString_channel3 = "beta_nsig,beta_nbkg";
string GlobalObsParamsString_channel3 = "glob_nsig,glob_nbkg";
string ObservableParamsString_channel3 = "mass";

myConfigurator->setNuisanceParams(channelname3, nuisanceParamsString_channel3);
myConfigurator->setGlobalObs(channelname3, GlobalObsParamsString_channel3);
myConfigurator->setObservables(channelname3, ObservableParamsString_channel3);

}

// CHANNEL: dielectron2011

if(run_channel2){

string nuisanceParamsString_channel2 = "";
if (usemassscaleuncer){
nuisanceParamsString_channel2 = "beta_nsig,beta_nbkg,beta_mass";
}
else{
nuisanceParamsString_channel2 = "beta_nsig,beta_nbkg";
}
string GlobalObsParamsString_channel2 = "glob_nsig,glob_nbkg,glob_mass";
string ObservableParamsString_channel2 = "mass";

myConfigurator->setNuisanceParams(channelname2, nuisanceParamsString_channel2);
myConfigurator->setGlobalObs(channelname2, GlobalObsParamsString_channel2);
myConfigurator->setObservables(channelname2, ObservableParamsString_channel2);

}

// CHANNEL: dielectron2012

if(run_channel4){

string nuisanceParamsString_channel4 = "";
if (usemassscaleuncer){
nuisanceParamsString_channel4 = "beta_nsig,beta_nbkg,beta_mass";
}
else{
nuisanceParamsString_channel4 = "beta_nsig,beta_nbkg";
}
string GlobalObsParamsString_channel4 = "glob_nsig,glob_nbkg,glob_mass";
string ObservableParamsString_channel4 = "mass";

myConfigurator->setNuisanceParams(channelname4, nuisanceParamsString_channel4);
myConfigurator->setGlobalObs(channelname4, GlobalObsParamsString_channel4);
myConfigurator->setObservables(channelname4, ObservableParamsString_channel4);

}

//setup combined RooWorkSpace and ModelConfig
myConfigurator->Setup();

// FIX MASS HYPOTHESIS
myConfigurator->setMassHypothesis(Zprimemass);

//safe combined WS

//myConfigurator->WriteCombinedWS("CombinedWS.root");

// ADJUST WORKSPACE FOR 7/8 TeV Combination

Pixie * myPixie = new Pixie();

if( (run_channel3 || run_channel4 ) && (run_channel1 || run_channel2) ){
   if (run_channel1){
      myPixie->SetupFor7and8Combination(myConfigurator->getCombinedWS(),"dimuon2011");
   }
   if (run_channel2){
      myPixie->SetupFor7and8Combination(myConfigurator->getCombinedWS(),"dielectron2011");
   }
}

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
   mydatapruner = new DataPruner(4.0, myConfigurator);
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



//just removing the uncertainties from the list of nuisance parameters does not switch off their variation in the Markov Chain (probaly because the prior term is included as a part of the likelihood)
if (!usemassscaleuncer){
   cout << "mass scale uncertainties are not applied!" << endl;
   if(run_channel2){
      myConfigurator->setVarRange("beta_mass_dielectron2011", -0.001, 0.001);
      myConfigurator->setVar("beta_mass_dielectron2011", 0.0);
   }
   if(run_channel4){
      myConfigurator->setVarRange("beta_mass_dielectron2012", -0.001, 0.001);
      myConfigurator->setVar("beta_mass_dielectron2012", 0.0);
   }
}

// RUN Bayesian limits

ResultatorWithBAT * myResultator = new ResultatorWithBAT(myConfigurator, mydatapruner, plotfile);
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
