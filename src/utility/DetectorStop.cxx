#include "DetectorStop.hxx"

#include "XMLUtility.hxx"

#include <algorithm>
#include <array>
#include <iostream>

DetectorStop::DetectorStop() {
  std::fill_n(ActiveExent,3,0);
  std::fill_n(CenterPosition,3,0);
  POTExposure=0xdeadbeef;
}

bool DetectorStop::ConfigureDetector(XMLNodePointer_t node) {
  std::array<bool, 3> found;
  ActiveExent[0] =
      GetXMLAttributeValue<double>(node, "DetectorActiveWidth_m", found[0]);
  ActiveExent[1] = GetXMLAttributeValue<double>(
      node, "DetectorActiveHeight_m", found[1]);
  ActiveExent[2] =
      GetXMLAttributeValue<double>(node, "DetectorActiveDepth_m", found[2]);

  if (std::count(found.begin(), found.end(), true) != 3) {
    std::cerr << "[ERROR]: When reading " << dstagname
              << " node, could not find all expected attributes: {"
              << "DetectorWidth_m: "
              << (found[0] ? "found" : "not found")
              << ", DetectorHeight_m: "
              << (found[1] ? "found" : "not found")
              << ", DetectorDepth_m: "
              << (found[2] ? "found" : "not found")
              << " }" << std::endl;

    return false;
  }
  return true;
}

DetectorStop DetectorStop::CloneDetectorConfig() {
  DetectorStop ds;

  std::copy_n(ActiveExent,3,ds.ActiveExent);
  std::copy_n(CenterPosition,3,ds.CenterPosition);
  ds.POTExposure = POTExposure;
  return ds;
}

bool DetectorStop::ConfigureStop(XMLNodePointer_t node) {
  std::array<bool, 4> found;
  CenterPosition[0] =
      GetXMLAttributeValue<double>(node, "LateralOffset_m", found[0]);
  CenterPosition[1] = GetXMLAttributeValue<double>(node, "VerticalOffset_m", found[1]);
  CenterPosition[2] = 0;
  POTExposure = GetXMLAttributeValue<double>(node, "POTExposure", found[3]);

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


std::vector<DetectorStop> ReadDetectorStopConfig(
    std::string const &fname, std::string const &RPName) {
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
      bool found;
      std::string name =
          GetXMLAttributeValue<std::string>(root_child, "Name", found);
      if (RPName.size()) {
        if (!found) {
          std::cout << "[INFO]: Ignoring unnamed run plan." << std::endl;
          root_child = xE.GetNext(root_child);
          continue;
        }
        if (RPName != name) {
          std::cout << "[INFO]: Ignoring run plan named: " << name << std::endl;
          root_child = xE.GetNext(root_child);
          continue;
        }
      }
      XMLNodePointer_t rp_child = xE.GetChild(root_child);
      DetectorStop detDefinition;
      bool foundDetDefinition = false;
      while (rp_child) {
        if (dstagname == xE.GetNodeName(rp_child)) {
          foundDetDefinition = detDefinition.ConfigureDetector(rp_child);
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
                std::cout << "[INFO]: Read stop at offset: " << ds.CenterPosition[0]
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
