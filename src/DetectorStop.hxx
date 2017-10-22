#ifndef DP_DETECTORSTOP_HXX_SEEN
#define DP_DETECTORSTOP_HXX_SEEN

#include "Utils.hxx"

#include "TH1D.h"
#include "TXMLEngine.h"

#include <array>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
static std::string const rptagname = "RunPlan";
static std::string const dstagname = "Detector";
static std::string const sstagname = "Stops";
static std::string const stagname = "Stop";
}

struct DetectorStop {
  double MeasurementRegionWidth;

  double DetectorFiducialWidth;
  double DetectorFiducialHeight;
  double DetectorFiducialDepth;
  double FiducialVolumeDensity;

  double LateralOffset;
  double POTExposure;

  std::map<int, std::vector<TH1D *> > Fluxes;
  std::map<int, std::map<std::string, std::vector<TH1D *> > > EventRates;

  DetectorStop()
      : MeasurementRegionWidth(0xdeadbeef),
        DetectorFiducialWidth(0xdeadbeef),
        DetectorFiducialHeight(0xdeadbeef),
        DetectorFiducialDepth(0xdeadbeef),
        FiducialVolumeDensity(0xdeadbeef),
        LateralOffset(0xdeadbeef),
        POTExposure(0xdeadbeef),
        Fluxes(),
        EventRates() {}

  bool ConfigureDetector(XMLNodePointer_t node) {
    std::array<bool, 5> found;
    MeasurementRegionWidth = GetXMLAttributeValue<double>(
        node, "MeasurementRegionWidth_m", found[0]);
    DetectorFiducialWidth =
        GetXMLAttributeValue<double>(node, "DetectorFiducialWidth_m", found[1]);
    DetectorFiducialHeight = GetXMLAttributeValue<double>(
        node, "DetectorFiducialHeight_m", found[2]);
    DetectorFiducialDepth =
        GetXMLAttributeValue<double>(node, "DetectorFiducialDepth_m", found[3]);
    FiducialVolumeDensity = GetXMLAttributeValue<double>(
        node, "FiducialVolumeDensity_kgm3", found[4]);

    if (std::count(found.begin(), found.end(), true) != 5) {
      std::cerr << "[ERROR]: When reading " << dstagname
                << " node, could not find all expected attributes: {"
                << "MeasurementRegionWidth_m: "
                << (found[0] ? "found" : "not found")
                << ", DetectorFiducialWidth_m: "
                << (found[1] ? "found" : "not found")
                << ", DetectorFiducialHeight_m: "
                << (found[2] ? "found" : "not found")
                << ", DetectorFiducialDepth_m: "
                << (found[3] ? "found" : "not found")
                << ", FiducialVolumeDensity_kgm3: "
                << (found[4] ? "found" : "not found") << " }" << std::endl;

      return false;
    }
    return true;
  }

  DetectorStop CloneDetectorConfig() {
    DetectorStop ds;

    ds.MeasurementRegionWidth = MeasurementRegionWidth;
    ds.DetectorFiducialWidth = DetectorFiducialWidth;
    ds.DetectorFiducialHeight = DetectorFiducialHeight;
    ds.DetectorFiducialDepth = DetectorFiducialDepth;
    ds.FiducialVolumeDensity = FiducialVolumeDensity;

    return ds;
  }
  bool ConfigureStop(XMLNodePointer_t node) {
    std::array<bool, 2> found;
    LateralOffset =
        GetXMLAttributeValue<double>(node, "LateralOffset_m", found[0]);
    POTExposure = GetXMLAttributeValue<double>(node, "POTExposure", found[1]);

    if (std::count(found.begin(), found.end(), true) != 2) {
      std::cerr << "[ERROR]: When reading " << stagname
                << " node, could not find all expected attributes: {"
                << "LateralOffset_m: " << (found[0] ? "found" : "not found")
                << ", POTExposure: " << (found[1] ? "found" : "not found")
                << "}" << std::endl;

      return false;
    }
    return true;
  }

  size_t GetNMeasurementSlices() {
    return floor(MeasurementRegionWidth / (DetectorFiducialWidth / 2.0)) * 2;
  }

  double GetAbsoluteOffsetOfSlice(size_t i) {
    if (i >= GetNMeasurementSlices()) {
      std::cout << "[WARN]: Requested position of slice: " << i
                << ", but current detector configuration only specifies: "
                << GetNMeasurementSlices() << std::endl;
      return 0xdeadbeef;
    }
    return (LateralOffset - ((DetectorFiducialWidth / 2.0) *
                             MeasurementRegionWidth * (double(i) + 0.5)));
  }

  void AddSliceFlux(size_t i, int species, TH1D const *flux) {
    if (i >= GetNMeasurementSlices()) {
      std::cout << "[WARN]: Requested position of slice: " << i
                << ", but current detector configuration only specifies: "
                << GetNMeasurementSlices() << std::endl;
      return;
    }

    if (!Fluxes.count(species)) {
      Fluxes[species] = std::vector<TH1D *>();
    }

    if (Fluxes[species].size() < i) {
      Fluxes[species].resize(i);
    }
    Fluxes[species][i] = static_cast<TH1D *>(flux->Clone());
    Fluxes[species][i]->SetDirectory(nullptr);
  }

  void PredictEventRates(std::vector<TH1D *> XSecComponents) { throw; }

  void Write() {
    TDirectory *ogDir = gDirectory;

    if (!ogDir) {
      std::cout
          << "[WARN]; No open TFile, cannot write detector stop information."
          << std::endl;
      return;
    }
    std::stringstream ss("");

    ss << "stop_" << LateralOffset << "_m";

    gDirectory->mkdir(ss.str().c_str());
    gDirectory->cd(ss.str().c_str());

    for (int species : {-14, -12, 12, 14}) {
      for (size_t fl_it = 0; fl_it < Fluxes[species].size(); ++fl_it) {
        TH1D *fl = Fluxes[species][fl_it];
        ss.str("");
        ss << GetSpeciesName(species) << "_flux_"
           << GetAbsoluteOffsetOfSlice(fl_it) << "_m";

        fl->Write(ss.str().c_str(), TObject::kOverwrite);
        fl->SetDirectory(nullptr);
      }

      if (!EventRates.count(species) || !EventRates[species].size()) {
        continue;
      }
      for (size_t evr_it = 0; evr_it < EventRates[species]["total"].size();
           ++evr_it) {
        TH1D *evr = EventRates[species]["total"][evr_it];
        ss.str("");
        ss << GetSpeciesName(species) << "_evrate_"
           << GetAbsoluteOffsetOfSlice(evr_it) << "_m";

        evr->Write(ss.str().c_str(), TObject::kOverwrite);
        evr->SetDirectory(nullptr);
      }
    }

    if (EventRates.size()) {
      gDirectory->mkdir("event_rates");
      gDirectory->cd("event_rates");

      for (int species : {-14, -12, 12, 14}) {
        if (!EventRates.count(species)) {
          continue;
        }
        for (std::map<std::string, std::vector<TH1D *> >::iterator evr_it =
                 EventRates[species].begin();
             evr_it != EventRates[species].end(); ++evr_it) {
          if (evr_it->first == "total") {
            continue;
          }

          for (size_t ev_it = 0; ev_it < evr_it->second.size(); ++ev_it) {
            TH1D *evr = evr_it->second[ev_it];
            ss.str("");
            ss << GetSpeciesName(species) << "_evrate_"
               << GetAbsoluteOffsetOfSlice(ev_it) << "_m_" << evr_it->first;

            evr->Write(ss.str().c_str(), TObject::kOverwrite);
            evr->SetDirectory(nullptr);
          }
        }
      }
    }
    if (ogDir) {
      ogDir->cd();
    }
  }

  ~DetectorStop() {
    for (auto species : Fluxes) {
      for (auto fl : species.second) {
        delete fl;
      }
    }
    for (auto species : EventRates) {
      for (auto evr_vect : species.second) {
        for (auto evr : evr_vect.second) {
          delete evr;
        }
      }
    }
  }
};

inline std::vector<DetectorStop> ReadDetectorStopConfig(
    std::string const &fname) {
  TXMLEngine xE;
  xE.SetSkipComments(true);
  XMLDocPointer_t doc = xE.ParseFile(fname.c_str());
  if (!doc) {
    std::cout << "[ERROR]: Attempted to parse XML file: " << fname
              << ", but failed." << std::endl;
    exit(1);
  }

  std::vector<DetectorStop> stops;

  XMLNodePointer_t rootNode = xE.DocGetRootElement(doc);

  XMLNodePointer_t root_child = xE.GetChild(rootNode);
  while (root_child) {  // Look for run plan node
    if (rptagname == xE.GetNodeName(root_child)) {
      XMLNodePointer_t rp_child = xE.GetChild(root_child);
      DetectorStop detDefinition;
      bool foundDetDefinition = false;
      while (rp_child) {
        if (dstagname == xE.GetNodeName(rp_child)) {
          foundDetDefinition = detDefinition.ConfigureDetector(rp_child);
          std::cout << "[INFO]: Read detector definition with "
                    << detDefinition.GetNMeasurementSlices()
                    << " measurement slices. { MeasurementRegionWidth: "
                    << detDefinition.MeasurementRegionWidth
                    << ", DetectorFiducialWidth: "
                    << detDefinition.DetectorFiducialWidth
                    << ", DetectorFiducialHeight: "
                    << detDefinition.DetectorFiducialHeight
                    << ", DetectorFiducialDepth: "
                    << detDefinition.DetectorFiducialDepth
                    << ", FiducialVolumeDensity: "
                    << detDefinition.FiducialVolumeDensity << " }."
                    << std::endl;

          rp_child = xE.GetNext(rp_child);
          continue;
        }
        if (sstagname == xE.GetNodeName(rp_child)) {
          if (!foundDetDefinition) {
            std::cout << "[WARN]: Ignoring " << sstagname << " node because "
                      << dstagname << " has not been encountered yet."
                      << std::endl;
            rp_child = xE.GetNext(rp_child);
            continue;
          }

          XMLNodePointer_t stops_child = xE.GetChild(rp_child);
          while (stops_child) {
            if (stagname == xE.GetNodeName(stops_child)) {
              DetectorStop ds = detDefinition.CloneDetectorConfig();
              if (ds.ConfigureStop(stops_child)) {
                std::cout << "[INFO]: Read stop at offset: " << ds.LateralOffset
                          << std::endl;
                stops.push_back(ds);
              } else {
                std::cout << "[WARN]: Failed to parse stop definition."
                          << std::endl;
              }
            }
            stops_child = xE.GetNext(stops_child);
          }
        }
        rp_child = xE.GetNext(rp_child);
      }
      break;
    }
    root_child = xE.GetNext(root_child);
  }

  xE.FreeDoc(doc);

  return stops;
}

#endif
