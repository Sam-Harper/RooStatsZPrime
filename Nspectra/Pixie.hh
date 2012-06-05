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


using namespace std;

class Pixie {

   public: 
      // constructor
      Pixie();
      // destructor
      ~Pixie();

      std::vector<std::string> commaSepStringToStringVec(std::string commaSepString);
      std::vector<std::string> commaSepStringToStringForRooArgSet(std::string channelname, std::string commaSepString, std::vector<std::string> sharedVars);

 
   private:

      const std::string legend;


};
#endif
