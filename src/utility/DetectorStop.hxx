#ifndef DP_DETECTORSTOP_HXX_SEEN
#define DP_DETECTORSTOP_HXX_SEEN

#include "TXMLEngine.h"

#include <string>
#include <vector>

namespace {
static std::string const rptagname = "RunPlan";
static std::string const dstagname = "Detector";
static std::string const sstagname = "Stops";
static std::string const stagname = "Stop";
}

struct DetectorStop {
  double ActiveExent[3];
  double CenterPosition[3];
  double POTExposure;

  DetectorStop();

  bool ConfigureDetector(XMLNodePointer_t node);
  DetectorStop CloneDetectorConfig();
  bool ConfigureStop(XMLNodePointer_t node);
};

std::vector<DetectorStop> ReadDetectorStopConfig(
    std::string const &fname, std::string const &RPName = "");

#endif
