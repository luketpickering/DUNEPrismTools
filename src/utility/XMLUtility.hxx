#ifndef XMLUTILITY_HXX_SEEN
#define XMLUTILITY_HXX_SEEN

#include "StringParserUtility.hxx"

#include "TXMLEngine.h"

#include <sstream>
#include <algorithm>
#include <cstdlib>

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
template <>
inline bool GetXMLAttributeValue<bool>(XMLNodePointer_t node,
                                       std::string const &attrName,
                                       bool &found) {
  TXMLEngine xE;
  xE.SetSkipComments(true);

  XMLAttrPointer_t attr = xE.GetFirstAttr(node);
  while (attr) {
    // Find match
    if (attrName == xE.GetAttrName(attr)) {
      std::string val = xE.GetAttrValue(attr);
      std::transform(val.begin(), val.end(), val.begin(), ::tolower);

      if ((val == "true")) {
        found = true;
        return true;
      }

      if ((val == "false")) {
        found = true;
        return false;
      }

      std::istringstream ss(val);
      bool rtn;
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
  return bool();
}
template <typename T>
inline T GetXMLAttributeValue(XMLNodePointer_t node,
                              std::string const &attrName) {
  bool dummy;
  return GetXMLAttributeValue<T>(node, attrName, dummy);
}

template <>
inline bool GetXMLAttributeValue<bool>(XMLNodePointer_t node,
                                       std::string const &attrName) {
  bool dummy;
  return GetXMLAttributeValue<bool>(node, attrName, dummy);
}

template <typename T>
std::vector<T> GetXMLAttributeList(XMLNodePointer_t node,
                                   std::string const &attrName, bool &found) {
  TXMLEngine xE;
  xE.SetSkipComments(true);

  XMLAttrPointer_t attr = xE.GetFirstAttr(node);
  while (attr) {
    // Find match
    if (attrName == xE.GetAttrName(attr)) {
      std::vector<T> rtn = ParseToVect<T>(xE.GetAttrValue(attr), ",");

      if (!rtn.size()) {
        break;
      }

      found = true;
      return rtn;
    }

    // Get Next Attribute
    attr = xE.GetNextAttr(attr);
  }

  found = false;
  return std::vector<T>();
}

template <typename T>
std::vector<T> GetXMLAttributeList(XMLNodePointer_t node,
                                   std::string const &attrName) {
  bool dummy;
  return GetXMLAttributeList<T>(node, attrName, dummy);
}

#endif
