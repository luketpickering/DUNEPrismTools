#ifndef EFFCORRECTOR_HXX_SEEN
#define EFFCORRECTOR_HXX_SEEN

#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

#include "TH2D.h"

#include <string>
#include <vector>

struct EffCorrector {

  enum ModeEnum { kEHadrAbsPos = 1, kEHadrDetPos, kENonNeutronHadrDetPos,
    kEHadrVisDetPos };

  std::vector<BoundingBox> StopActiveRegions;

  TH2 *MuonKinematics_seleff;
  TH2 *ShowerKinematics_seleff;

  ModeEnum mode;

  EffCorrector(ModeEnum mode, std::string const &EffCorrectorFile,
   StopConfig &sc);

  double GetMuonKinematicsEffWeight(DepositsSummary const & edr);
  double GetHadronKinematicsEffWeight(DepositsSummary const & edr);
};

#endif
