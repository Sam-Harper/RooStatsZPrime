//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#ifndef PIXIE_hh
#define PIXIE_hh


#include<vector>
#include<string>
#include<map>
#include "RooWorkspace.h"


#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

class Pixie {

   public: 
      // constructor
      Pixie();
      // destructor
      ~Pixie();

      std::vector<std::string> commaSepStringToStringVec(std::string commaSepString);
      std::vector<std::string> commaSepStringToStringForRooArgSet(std::string channelname, std::string commaSepString, std::vector<std::string> sharedVars);
      void SetupFor7and8Combination(RooWorkspace * ws, std::string channelname);
      void SetupFor7and8Combination_RS(RooWorkspace * ws, std::string channelname);
      void SetResolution_RS(RooWorkspace * ws, double kdivMplred, std::vector<std::string> channelnames);
      void SetResolution_ZPSI(RooWorkspace * ws, std::vector<std::string> channelnames);
      void CheckLikelihood(RooWorkspace * ws, ModelConfig * modelconfig);
 
   private:

      const std::string legend;


};
#endif
