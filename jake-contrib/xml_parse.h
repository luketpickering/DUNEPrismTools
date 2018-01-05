#include "TXMLEngine.h"
#include <array>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>


struct DetectorStop {
  double detectorSizeX;
  double detectorSizeY;
  double detectorSizeZ;

  double fiducialGapX;
  double fiducialGapY;
  double fiducialGapZ;

  double shift;
};

template <typename T>
inline T GetXMLAttributeValue(XMLNodePointer_t node,
                              std::string const &attrName, bool &found) {
  TXMLEngine xE;
  xE.SetSkipComments(true);

  XMLAttrPointer_t attr = xE.GetFirstAttr(node);
  while (attr) {
    // Find match
    if (attrName == xE.GetAttrName(attr)) {
      std::istringstream ss(xE.GetAttrValue(attr));
      T rtn;
      ss >> rtn;

      if (ss.fail()) {
        break;
      }

      found = true;
      return rtn;
    }

    // Get Next Attribute
    attr = xE.GetNextAttr(attr);
  }

  found = false;
  return T();
}

static std::string const rptagname = "RunPlan";
static std::string const dstagname = "Detector";
static std::string const sstagname = "Stops";
static std::string const stagname = "Stop";
static std::string const fullname = "Full";

inline std::vector<DetectorStop> ReadDetectorStopConfig(std::string const &fname, std::string const &RPName = ""){
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

  while (root_child){
    if (rptagname == xE.GetNodeName(root_child) ){
      bool found_rp;
      std::string name = GetXMLAttributeValue<std::string>(root_child, "Name", found_rp);
      
      if (RPName.size()) {
        if(!found_rp){
          std::cout<<"Ignoring unnamed run plan." << std::endl;
          root_child = xE.GetNext(root_child);
        }
        if(RPName != name) {
          std::cout << "Ignoring run plan named: " << name << std::endl;
          root_child = xE.GetNext(root_child);
          continue;
        }
      }
    

      XMLNodePointer_t rp_child = xE.GetChild(root_child);
      DetectorStop detDefinition;
//      bool foundDetDefinition = false;
     
      std::array<bool, 7> found;

      while(rp_child){
        if ( dstagname == xE.GetNodeName(rp_child) ){
          detDefinition.detectorSizeX = GetXMLAttributeValue<double>(rp_child,"detectorSizeX",found[0]);
          detDefinition.detectorSizeY = GetXMLAttributeValue<double>(rp_child,"detectorSizeY",found[1]);
          detDefinition.detectorSizeZ = GetXMLAttributeValue<double>(rp_child,"detectorSizeZ",found[2]);

          detDefinition.fiducialGapX = GetXMLAttributeValue<double>(rp_child,"fiducialGapX",found[3]); 
          detDefinition.fiducialGapY = GetXMLAttributeValue<double>(rp_child,"fiducialGapY",found[4]); 
          detDefinition.fiducialGapZ = GetXMLAttributeValue<double>(rp_child,"fiducialGapZ",found[5]); 

          detDefinition.shift = GetXMLAttributeValue<double>(rp_child,"shift",found[6]); 

          std::cout << detDefinition.detectorSizeX<<std::endl<< 
          detDefinition.detectorSizeY<<std::endl<<
          detDefinition.detectorSizeZ<<std::endl<<

          detDefinition.fiducialGapX <<std::endl<<
          detDefinition.fiducialGapY <<std::endl<<
          detDefinition.fiducialGapZ <<std::endl<<

          detDefinition.shift<<std::endl;

          stops.push_back(detDefinition);

          rp_child = xE.GetNext(rp_child);
          continue;
        }
        rp_child = xE.GetNext(rp_child);     
      }
      break;
    }
    root_child = xE.GetNext(root_child);
  }

  return stops;
}

inline DetectorStop GetFullDetectorConfig(std::string const &fname, std::string const &RPName = ""){
  TXMLEngine xE;
  xE.SetSkipComments(true);
  XMLDocPointer_t doc = xE.ParseFile(fname.c_str());

  if (!doc) {
    std::cout << "[ERROR]: Attempted to parse XML file: " << fname
                  << ", but failed." << std::endl;
                  exit(1);
  }


  XMLNodePointer_t rootNode = xE.DocGetRootElement(doc);

  XMLNodePointer_t root_child = xE.GetChild(rootNode);
  DetectorStop detDefinition;

  while (root_child){
    if (fullname == xE.GetNodeName(root_child) ){
      bool found_rp;
      std::string name = GetXMLAttributeValue<std::string>(root_child, "Name", found_rp);
      
      if (RPName.size()) {
        if(!found_rp){
          std::cout<<"Ignoring unnamed run plan." << std::endl;
          root_child = xE.GetNext(root_child);
        }
        if(RPName != name) {
          std::cout << "Ignoring run plan named: " << name << std::endl;
          root_child = xE.GetNext(root_child);
          continue;
        }
      }
    

      XMLNodePointer_t rp_child = xE.GetChild(root_child);
//      bool foundDetDefinition = false;
     
      std::array<bool, 7> found;

      while(rp_child){
        if ( dstagname == xE.GetNodeName(rp_child) ){
          detDefinition.detectorSizeX = GetXMLAttributeValue<double>(rp_child,"detectorSizeX",found[0]);
          detDefinition.detectorSizeY = GetXMLAttributeValue<double>(rp_child,"detectorSizeY",found[1]);
          detDefinition.detectorSizeZ = GetXMLAttributeValue<double>(rp_child,"detectorSizeZ",found[2]);

          detDefinition.fiducialGapX = GetXMLAttributeValue<double>(rp_child,"fiducialGapX",found[3]); 
          detDefinition.fiducialGapY = GetXMLAttributeValue<double>(rp_child,"fiducialGapY",found[4]); 
          detDefinition.fiducialGapZ = GetXMLAttributeValue<double>(rp_child,"fiducialGapZ",found[5]); 

          detDefinition.shift = GetXMLAttributeValue<double>(rp_child,"shift",found[6]); 

          std::cout << "FOUND FULL DETECTOR\n" <<
          detDefinition.detectorSizeX<<std::endl<< 
          detDefinition.detectorSizeY<<std::endl<<
          detDefinition.detectorSizeZ<<std::endl<<

          detDefinition.fiducialGapX <<std::endl<<
          detDefinition.fiducialGapY <<std::endl<<
          detDefinition.fiducialGapZ <<std::endl<<

          detDefinition.shift<<std::endl;


          rp_child = xE.GetNext(rp_child);
          continue;
        }
        rp_child = xE.GetNext(rp_child);     
      }
      break;
    }
    root_child = xE.GetNext(root_child);
  }


  return detDefinition;
}
