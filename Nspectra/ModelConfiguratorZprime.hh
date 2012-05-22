//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef MODELCONFIGURATORZPRIME_hh
#define MODELCONFIGURATORZPRIME_hh


#include<vector>
#include<string>
#include<map>

#include "RooWorkspace.h"
#include "ModelConfigurator.hh"

using namespace RooFit;
//using namespace RooStats;
using namespace std;

class ModelConfiguratorZprime : public ModelConfigurator {

   public: 
      // constructor
      ModelConfiguratorZprime(std::map<string, RooWorkspace*>  WSmap);
      // destructor
      ~ModelConfiguratorZprime();

      // creates combined workspace and ModelConfig for the analysis
      void Setup();

      void Setup_channel(std::string channelname);
      void Setup_ModelConfig();

      void setMassHypothesis(double peak); // set the Z' mass hypothesis
      double getMassHypothesis(); // get the Z' mass hypothesis
 
   private:

      const std::string sublegend;

      double _peak; // the Z' mass hypothesis



};
#endif