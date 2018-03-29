#include "StopDimensions.hxx"

#include "SimConfigTreeReader.hxx"

#include <vector>
#include <iostream>
#include <algorithm>

// #define DEBUG

StopDimensions::StopDimensions(std::string &fname) : NYSteps(3),
NZSteps(3) {

  SimConfig sc("SimConfigTree", fname);

  NXSteps = sc.NXSteps;
  std::copy_n(sc.DetMin, 3, DetMin);
  std::copy_n(sc.DetMax, 3, DetMax);
  std::copy_n(sc.VetoGap, 3, VetoGap);

}

TH3D* StopDimensions::BuildDetectorMap() {
  std::vector<double> YBins = {DetMin[1], DetMin[1] + VetoGap[1],
                               DetMax[1] - VetoGap[1], DetMax[1]};
  std::vector<double> ZBins = {DetMin[2], DetMin[2] + VetoGap[2],
                               DetMax[2] - VetoGap[2], DetMax[2]};

  std::vector<double> XBins;
  double step =
      ((DetMax[0] - VetoGap[0]) - (DetMin[0] + VetoGap[0])) / double(NXSteps - 2);
  XBins.push_back(DetMin[0]);
  XBins.push_back(DetMin[0] + VetoGap[0]);

#ifdef DEBUG
  std::cout << "Building detector map: X { " << XBins[0] << ", " << XBins.back()
            << ", " << std::flush;
#endif

  for (int i = 0; i < (NXSteps - 2); ++i) {
    XBins.push_back(XBins.back() + step);
#ifdef DEBUG
    std::cout << XBins.back() << ((i != (NXSteps - 1)) ? ", " : "")
              << std::flush;
#endif
  }
  XBins.push_back(XBins.back() + VetoGap[0]);
#ifdef DEBUG

  std::cout << XBins.back() << "}" << std::endl;
  std::cout << "\t\t: Y {" << YBins[0] << ", " << YBins[1] << ", " << YBins[2]
            << ", " << YBins[3] << " }" << std::endl;
  std::cout << "\t\t: Z {" << ZBins[0] << ", " << ZBins[1] << ", " << ZBins[2]
            << ", " << ZBins[3] << " }" << std::endl;
#endif

  if (NYSteps == 0) {
    NYSteps = 3;
  }
  if (NZSteps == 0) {
    NZSteps = 3;
  }

  if (NYSteps != 3) {
    YBins.clear();
    double step_y = (DetMax[1] - DetMin[1]) / double(NYSteps);
    YBins.push_back(DetMin[1]);
    for (int i = 0; i < NYSteps; ++i) {
      YBins.push_back(YBins.back() + step_y);
    }
  }

  if (NZSteps != 3) {
    ZBins.clear();
    double step_z = (DetMax[2] - DetMin[2]) / double(NZSteps);
    ZBins.push_back(DetMin[2]);
    for (int i = 0; i < NYSteps; ++i) {
      ZBins.push_back(ZBins.back() + step_z);
    }
  }

  TH3D* Dm =
      new TH3D("dm", "", (XBins.size() - 1), XBins.data(), (YBins.size() - 1),
               YBins.data(), (ZBins.size() - 1), ZBins.data());
  Dm->SetDirectory(nullptr);
  Dm->SetName("dm_c");

  return Dm;
}
