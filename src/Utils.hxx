#ifndef DP_UTILS_HXX_SEEN
#define DP_UTILS_HXX_SEEN

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TXMLEngine.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template <typename T>
inline T str2T(std::string const &str) {
  std::istringstream stream(str);
  T d;
  stream >> d;

  if (stream.fail()) {
    std::cerr << "[WARN]: Failed to parse string: " << str
              << " as requested type." << std::endl;
    return T();
  }

  return d;
}

template <typename T>
inline std::string to_str(T const &inp) {
  std::stringstream stream("");
  stream << inp;
  return stream.str();
}

template <>
inline bool str2T<bool>(std::string const &str) {
  if ((str == "true") || (str == "True") || (str == "TRUE")) {
    return true;
  }

  if ((str == "false") || (str == "False") || (str == "FALSE")) {
    return false;
  }

  std::istringstream stream(str);
  bool d;
  stream >> d;

  if (stream.fail()) {
    std::cerr << "[WARN]: Failed to parse string: " << str
              << " as requested type." << std::endl;
    return false;
  }

  return d;
}

template <typename T>
inline std::vector<T> ParseToVect(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<T> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) {
      continue;
    }
    vals.push_back(str2T<T>(temp_string));
  }
  return vals;
}

template <typename T>
inline void AppendVect(std::vector<T> &target, std::vector<T> const &toApp) {
  for (size_t i = 0; i < toApp.size(); ++i) {
    target.push_back(toApp[i]);
  }
}

// Converts "5_10:1" into a vector containing: 5,6,7,8,9,10
inline std::vector<double> BuildDoubleList(std::string const &str) {
  std::vector<std::string> steps = ParseToVect<std::string>(str, ":");
  if (steps.size() != 2) {
    return ParseToVect<double>(str, ",");
  }
  double step = str2T<double>(steps[1]);

  std::vector<double> range = ParseToVect<double>(steps[0], "_");
  if (steps.size() != 2) {
    std::cout
        << "[ERROR]: When attempting to parse bin range descriptor: \" " << str
        << "\", couldn't determine range. Expect form: <bin1low>_<binXUp>:step"
        << std::endl;
    exit(1);
  }

  int nsteps = (range[1] - range[0]) / step;

  std::vector<double> rtn;
  for (int step_it = 0; step_it <= nsteps; ++step_it) {
    rtn.push_back(range[0] + step * step_it);
  }
  return rtn;
}

inline TFile *CheckOpenFile(std::string const &fname, char const *opts="") {
  TFile *inpF = new TFile(fname.c_str(),opts);
  if (!inpF || !inpF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open input file: " << fname << std::endl;
    exit(1);
  }
  return inpF;
}

template <class TH>
inline TH *GetHistogram(TFile *f, std::string const &fhname) {
  TH *inpH = dynamic_cast<TH *>(f->Get(fhname.c_str()));

  if (!inpH) {
    std::cout << "[ERROR]: Couldn't get TH: " << fhname
              << " from input file: " << f->GetName() << std::endl;
    exit(1);
  }

  inpH = static_cast<TH *>(inpH->Clone());
  inpH->SetDirectory(nullptr);
  return inpH;
}

template <class TH>
inline TH *GetHistogram(std::string const &fname, std::string const &hname) {
  TDirectory *ogDir = gDirectory;

  TFile *inpF = CheckOpenFile(fname);

  TH *h = GetHistogram<TH>(inpF, hname);

  inpF->Close();
  delete inpF;

  if (ogDir) {
    ogDir->cd();
  }

  return h;
}

template <class TH>
inline std::vector<TH *> GetHistograms(std::string const &fname,
                                         std::string const &hnamePattern) {
  std::vector<TH *> histos;

  TFile *inpF = CheckOpenFile(fname);

  TRegexp matchExp(hnamePattern.c_str(), true);

  TIter next(inpF->GetListOfKeys());
  TKey *k;
  while ((k = dynamic_cast<TKey *>(next()))) {
    TObject *obj = k->ReadObj();
    TH *obj_hist = dynamic_cast<TH *>(obj);
    if (!obj_hist) {
      continue;
    }

    Ssiz_t len = 0;
    if (matchExp.Index(obj->GetName(), &len) != Ssiz_t(-1)) {
      histos.push_back(GetHistogram<TH>(inpF, obj->GetName()));
    }
  }

  inpF->Close();
  delete inpF;

  return histos;
}

inline void SumHistograms(TH1D *summed, double *coeffs,
                          std::vector<TH1D *> const &histos) {
  size_t nf = histos.size();
  summed->Reset();
  for (size_t f_it = 0; f_it < nf; ++f_it) {
    summed->Add(histos[f_it], coeffs[f_it]);

    for (Int_t bi_it = 0; bi_it < summed->GetXaxis()->GetNbins() + 1; ++bi_it) {
      double bc = summed->GetBinContent(bi_it);
      double be = summed->GetBinError(bi_it);

      if ((bc != bc) || (be != be)) {
        std::cout << "Produced bad histo: bin " << bi_it
                  << " after adding coeff " << f_it << " = " << coeffs[f_it]
                  << std::endl;
        throw;
      }
    }
  }
}

inline double FindHistogramPeak(TH1D *hist, double resolution,
                                std::string const &WriteName) {
  double min = hist->GetXaxis()->GetBinLowEdge(1);
  double max = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());

  Int_t NSteps = ceil((max - min) / resolution);
  max = min + NSteps * resolution;

  TGraph interp(1);
  interp.Set(hist->GetXaxis()->GetNbins());

  for (Int_t bi_it = 1; bi_it < hist->GetXaxis()->GetNbins(); ++bi_it) {
    interp.SetPoint(bi_it - 1, hist->GetXaxis()->GetBinCenter(bi_it),
                    hist->GetBinContent(bi_it));
  }

  double peak = std::numeric_limits<double>::min();
  double peak_E = 0;

  for (Int_t i = 0; i < NSteps; ++i) {
    if (interp.Eval(min + i * resolution) > peak) {
      peak = interp.Eval(min + i * resolution);
      peak_E = min + i * resolution;
    }
  }

  if (WriteName.length()) {
    interp.Write(WriteName.c_str(), TObject::kOverwrite);
  }

  return peak_E;
}

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
std::vector<T> ParseStringList(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<T> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) {
      continue;
    }
    std::istringstream ss(temp_string);
    vals.push_back(T());
    ss >> vals.back();
  }

  return vals;
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
      std::vector<T> rtn = ParseStringList<T>(xE.GetAttrValue(attr), ",");

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

std::string GetSpeciesName(int pdg) {
  switch (pdg) {
    case -14:
      return "numubar";
    case -12:
      return "nuebar";
    case 12:
      return "nue";
    case 14:
      return "numu";
    default:
      return "unknown";
  }
}

#endif
