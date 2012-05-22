//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "Pixie.hh"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

using namespace std;


//--------------------Constructor
Pixie::Pixie() : legend("Pixie"){
;
}

Pixie::~Pixie() {
;
}

std::vector<std::string> Pixie::commaSepStringToStringVec(std::string commaSepString) {
   std::vector<std::string> stringVec;
   boost::split(stringVec, commaSepString, boost::is_any_of(","));
   return stringVec;
}

std::vector<std::string> Pixie::commaSepStringToStringForRooArgSet(std::string channelname, std::string commaSepString, std::vector<std::string> sharedVars) {
   std::string funclegend = "commaSepStringToStringVecPlusChannelName() - ";
   funclegend += channelname;

   std::vector<std::string> stringVec;
   stringVec = commaSepStringToStringVec(commaSepString);

      //check if shared variables have been defined
      if( sharedVars.empty() ){
         cout << legend << funclegend << " WARNING: shared variables have not yet been set (only ok if there are none)" << endl;
      }

   for(std::vector<std::string>::iterator it= stringVec.begin(); it != stringVec.end(); it++ ){
      //check if variable is contained in the sharedVars vector
      if(std::find(sharedVars.begin(), sharedVars.end(), *(it)) == sharedVars.end()){
         (*it)=(*it)+"_"+channelname;
      }
   }

   return stringVec;
}
