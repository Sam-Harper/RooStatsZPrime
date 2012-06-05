//
// authors: Gena Kukartsev, Stefan A. Schmitz
//
// project started: April 2012
// 

#include "PoiRangeEstimator.hh"
#include "TStopwatch.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>

#include <numeric>
#include <iterator>

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "TH1F.h"
#include "RooRealVar.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;


//--------------------Constructor
PoiRangeEstimator::PoiRangeEstimator(ModelConfiguratorZprime* configurator, Resultator * myResultator) : legend("PoiRangeEstimator - ") {
   cout << legend << " constructed" << endl;

   _configurator = configurator;
   _myResultator = myResultator;
   _myWS = _configurator->getCombinedWS() ;
   _myModelConfigCopy =  new ModelConfig(*(_configurator->getCombinedModelConfig())); //COMMENT: not sure why we need a copy for this, keep it for now as in twobody.C

   mpPlrInt = 0;
   data = _myResultator->getDataBox()->createObservedData();
}

PoiRangeEstimator::PoiRangeEstimator(ModelConfiguratorZprime* configurator, RooDataSet * inputdata) : legend("PoiRangeEstimator - ") {
   cout << legend << " constructed" << endl;

   _configurator = configurator;
   _myWS = _configurator->getCombinedWS() ;
   _myModelConfigCopy =  new ModelConfig(*(_configurator->getCombinedModelConfig())); //COMMENT: not sure why we need a copy for this, keep it for now as in twobody.C

   mpPlrInt = 0;
   _myResultator = 0; 
   data = inputdata;
}

PoiRangeEstimator::~PoiRangeEstimator() {
  delete mpPlrInt;
  delete _myModelConfigCopy;
}

Double_t PoiRangeEstimator::GetPoiUpper(std::string channel, Double_t peak){
  //
  // Estimate a good value for the upper boundary of the range of POI
  // this is known to work well for up to two channels -> should still be a reasonable estimate for 4 channels

   std::string funclegend = "GetPoiUpper(std::string channel, Double_t peak) - ";
   cout << legend << funclegend << endl;

  double _range = 1.0;

   //FIXME: remove hardcoded channel names!!!

  if (channel.find("dimuon2011")!=std::string::npos || 
      channel.find("dielectron2011")!=std::string::npos){
    _range = GetPoiUpperSimple(channel, peak);
  }
  else if (channel.find("multi")!=std::string::npos){
   std::vector<double> rangevec;
   std::vector<std::string> channelnamevec = _configurator->getChannelNames();
   for (std::vector<std::string>::iterator vecIt = channelnamevec.begin(); vecIt != channelnamevec.end(); vecIt++){
      rangevec.push_back(GetPoiUpperSimple(*vecIt, peak));
   }
    _range = *std::max_element( rangevec.begin(), rangevec.end() ); 

  }
  else{
    std::cout << legend << "Unknown channel "
         << channel << " , exiting"
         << std::endl;
    std::exit(-1);
  }

  std::cout << legend << "crude estimate for poi range: "
       << 3.0*_range << std::endl;
  std::cout << legend 
       << "this will be used if the profile likelihood ratio estimate fails"
       << std::endl;

  Double_t result = 0.1;

  // estimate limit with profile likelihood ratio and
  // set the range to 3 times the limit
  GetPlrInterval(0.95);
  
  // query intervals
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooRealVar * _poi = (RooRealVar *)_myModelConfigCopy->GetParametersOfInterest()->first();
  double upper_limit = mpPlrInt->UpperLimit( *_poi );
  RooMsgService::instance().setGlobalKillBelow(msglevel);
  
  // safety in case upper limit == 0
  if (upper_limit<std::numeric_limits<double>::min()){
    upper_limit = _range;
  }

  result = 3.0*upper_limit;

  return result;
}

double PoiRangeEstimator::GetPoiUpperSimple(std::string channel, Double_t peak){
  //
  // Estimate a good value for the upper boundary of the range of POI
  // using just data in a window corresponding to the signal region
  //

   std::string funclegend = "GetPoiUpperSimple(std::string channel, Double_t peak) - ";
   cout << legend << funclegend << endl;

  char buf[128];

  double _range = 1.0;

  // get data yield under the signal peak
  _myWS->var("peak")->setVal(peak);
  sprintf(buf,"width_%s",channel.c_str());
  double _width = _myWS->function(buf)->getVal();
  sprintf(buf,"sigma_%s",channel.c_str());
  double _sigma = _myWS->function(buf)->getVal();
  double c_min = peak - sqrt(_width*_width + _sigma*_sigma);
  double c_max = peak + sqrt(_width*_width + _sigma*_sigma);
  TH1F h_count("h_count", "h_count", 100, c_min, c_max);
  // mData[channel]->fillHistogram(&h_count, *_myModelConfigCopy->GetObservables()); //COMMENT: no idea what this line is good for, so I commented it out
  double n_count = h_count.Integral(1,100);

  // ad-hoc fix when there are no events in window
  if (n_count < 0.3) n_count = 0.3;
  std::cout << legend << "event yield in data in the window: "
       << n_count << std::endl;

  // compute the corresponding POI range
  sprintf(buf,"nsig_scale_%s",channel.c_str());
  double _nsig_scale = _myWS->var(buf)->getVal();
  sprintf(buf,"nz_%s",channel.c_str());
  double _nz = _myWS->var(buf)->getVal();
  sprintf(buf,"eff_%s",channel.c_str());
  double _eff = _myWS->function(buf)->getVal();
  double n_excess = 3.0 * sqrt(n_count)/0.68; // let signal excess be ~ uncertainty on BG
  _range = n_excess / _nsig_scale / _nz / _eff;

  return _range;
}

LikelihoodInterval * PoiRangeEstimator::GetPlrInterval( double conf_level ){
  //
  // Profile likelihood ratio interval calculations
  //

  delete mpPlrInt;
  
  //_myModelConfigCopy->Print(); // seems like unnecessary output at this point

  ProfileLikelihoodCalculator plc(*data, *_myModelConfigCopy);
  plc.SetConfidenceLevel(conf_level);
  mpPlrInt = plc.GetInterval();

  return mpPlrInt;
}