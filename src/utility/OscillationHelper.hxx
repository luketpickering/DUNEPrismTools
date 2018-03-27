#ifndef DP_OSCILLATIONHELPER_HXX_SEEN
#define DP_OSCILLATIONHELPER_HXX_SEEN

#include "BargerPropagator.h"

struct OscillationHelper {
  enum NuTypes {
    kNuebarType = -1,
    kNumubarType = -2,
    kNutaubarType = -3,
    kNueType = 1,
    kNumuType = 2,
    kNutauType = 3,
  };

  double DipAngle_degrees;      // = 5.8;
  double OscParams[6];  // = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
  double LengthParam;   // = 0xdeadbeef;

  bool IsSetUp;

  BargerPropagator bp;
  NuTypes FromType, ToType;

  NuTypes GetNuType(int pdg);

  void Setup(std::string const &FileWithConfTree);
  void Setup(double OscParams[6], double DipAngle_degrees=5.8);

  OscillationHelper() : IsSetUp(false) {};

  void SetOscillationChannel(int PDGFrom, int PDGTo);
  double GetWeight(double ENu_GeV);
};

#endif
