#ifndef DP_OSCILLATIONHELPER_HXX_SEEN
#define DP_OSCILLATIONHELPER_HXX_SEEN

#include "BargerPropagator.h"

#include "TDirectory.h"
#include "TH1D.h"

#ifdef USE_FHICL
#include "fhiclcpp/ParameterSet.h"
#endif

#include <memory>

struct OscillationHelper {
  enum NuTypes {
    kNuebarType = -1,
    kNumubarType = -2,
    kNutaubarType = -3,
    kNueType = 1,
    kNumuType = 2,
    kNutauType = 3,
  };

  double DipAngle_degrees; // = 5.8;
  double OscParams[6];     // = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
  double LengthParam;      // = 0xdeadbeef;

  bool IsSetUp;

  BargerPropagator bp;
  Int_t FromPDG, ToPDG;
  NuTypes FromType, ToType;

  NuTypes GetNuType(int pdg);

  void Setup(std::string const &FileWithConfTree);
  void Setup(double OscParams[6], double DipAngle_degrees = 5.8);
#ifdef USE_FHICL
  void Setup(fhicl::ParameterSet const &);
#endif

  OscillationHelper() : IsSetUp(false){};
  OscillationHelper(OscillationHelper const &other) {
    DipAngle_degrees = other.DipAngle_degrees;
    std::copy_n(other.OscParams, 6, OscParams);
    LengthParam = other.LengthParam;
    IsSetUp = other.IsSetUp;
    FromPDG = other.FromPDG;
    ToPDG = other.ToPDG;
    FromType = other.FromType;
    ToType = other.ToType;
  }

  void SetOscillationChannel(int PDGFrom, int PDGTo);
  double GetWeight(double ENu_GeV);
  void OscillateHistogram(std::unique_ptr<TH1D> &h);

  void WriteConfigTree(TDirectory *f);
};

#endif
