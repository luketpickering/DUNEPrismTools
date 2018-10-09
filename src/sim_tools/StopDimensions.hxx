#ifndef STOPDIMENSIONS_HXX_SEEN
#define STOPDIMENSIONS_HXX_SEEN

#include "TH3D.h"

struct StopDimensions {
  int NXSteps;
  double DetMin[3];
  double DetMax[3];
  double VetoGap[3];
  int NYSteps;
  int NZSteps;
  StopDimensions()
      : NXSteps(400),
        DetMin{0, 0, 0},
        DetMax{0, 0, 0},
        VetoGap{0, 0, 0},
        NYSteps(3),
        NZSteps(3) {}
  StopDimensions(std::string &fname);
  TH3D* BuildDetectorMap();
};

#endif
