#pragma once

#include "TH2Jagged.h"

#include "ROOTUtility.hxx"

std::vector<Double_t> ParseJaggedUniformAxis(fhicl::ParameterSet const &ps);

/// Builds a TH2Jagged from a supplied fhicl::ParameterSet
///
/// Format for the Parameter set:
// {
//   xattr: "MYX"
//   yattr: "MYY"
//   xisuniform: true # default
//   binning: [
//     {
//       MYX: "0_1:0.5,1_10:1"
//       MYY: "0_5:0.25"
//     },
//     {
//       MYX: "10_20:1.5"
//       MYY: "0_5:0.5"
//     }, ...
//   ]
//   Title: "MyJagged;xtitle;ytitle" #optional
//   Name: "MyJagged" #optional
// }
template <class TH2T>
inline TH2Jagged<TH2T> *BuildJagged(fhicl::ParameterSet const &ps) {
  bool xisunform = ps.get<bool>("xisuniform", true);
  std::string const &uattr = ps.get<std::string>(xisunform ? "xattr" : "yattr");
  std::string const &nuattr =
      ps.get<std::string>(xisunform ? "yattr" : "xattr");

  std::vector<fhicl::ParameterSet> const &binning =
      ps.get<std::vector<fhicl::ParameterSet>>("binning");

  std::vector<std::vector<Double_t>> UBinnings;
  std::vector<std::vector<Double_t>> NUBinnings;

  for (fhicl::ParameterSet const &NUSlice : binning) {
    UBinnings.push_back(BuildBinEdges(NUSlice.get<std::string>(uattr)));
    NUBinnings.push_back(BuildBinEdges(NUSlice.get<std::string>(nuattr)));
  }

  std::vector<Double_t> UBinEdges = ParseJaggedUniformAxis(ps);

  // Stack NU binnings
  std::vector<Int_t> NUNBins;
  std::vector<std::vector<Double_t>> NUBinEdges;
  for (size_t sit = 0; sit < UBinnings.size(); ++sit) {
    for (size_t uit = 0; uit < (UBinnings.at(sit).size() - 1); ++uit) {
      NUNBins.push_back(NUBinnings.at(sit).size() - 1);
      NUBinEdges.push_back(NUBinnings.at(sit));
    }
  }

  std::vector<Double_t const *> NUBinEdges_2D;
  for (std::vector<Double_t> const &nube : NUBinEdges) {
    NUBinEdges_2D.push_back(nube.data());
  }

  std::string name = ps.get<std::string>("Name", "MyJagged");
  std::string title = ps.get<std::string>("Title", "");

  return new TH2Jagged<TH2T>(name.c_str(), title.c_str(), UBinEdges.size() - 1,
                             UBinEdges.data(), NUNBins.data(),
                             NUBinEdges_2D.data(), xisunform);
}

template <typename ST> struct hist_types {};

template <> struct hist_types<double> {
  using TH2T = TH2D;
  using TH1T = TH1D;
};

template <> struct hist_types<float> {
  using TH2T = TH2F;
  using TH1T = TH1F;
};

template <typename ST>
inline TH1 *BuildHistogram(fhicl::ParameterSet const &ps,
                           std::string const &Name, std::string const &Title,
                           std::string const &xattr = "energy",
                           std::string const &yattr = "off_axis",
                           TDirectory *dir = nullptr) {
  if (ps.has_key("xattr")) {

    TH1 *j = BuildJagged<ST>(ps);
    j->SetName(Name.c_str());
    j->SetTitle(Title.c_str());
    return j;

  } else if (ps.has_key(xattr) && ps.has_key(yattr)) {

    std::vector<double> XBinning = BuildBinEdges(ps.get<std::string>(xattr));
    std::vector<double> YBinning = BuildBinEdges(ps.get<std::string>(yattr));

    if (XBinning.size() < 2) {
      std::cout << "[ERROR]: Failed to find any X bins in xbinning descriptor: "
                << ps.get<std::string>(xattr) << std::endl;
      throw;
    }

    if (YBinning.size() < 2) {
      TH1 *h = new typename hist_types<ST>::TH1T(
          Name.c_str(), Title.c_str(), (XBinning.size() - 1), XBinning.data());
      h->SetDirectory(dir);
      return h;
    }

    TH1 *h = new typename hist_types<ST>::TH2T(
        Name.c_str(), Title.c_str(), (XBinning.size() - 1), XBinning.data(),
        (YBinning.size() - 1), YBinning.data());

    h->SetDirectory(dir);
    return h;

  } else if (ps.has_key(xattr)) {

    std::vector<double> XBinning = BuildBinEdges(ps.get<std::string>(xattr));
    TH1 *h = new typename hist_types<ST>::TH1T(
        Name.c_str(), Title.c_str(), (XBinning.size() - 1), XBinning.data());
    h->SetDirectory(dir);
    return h;
  }
  std::cout << "[ERROR]: Failed to build histogram as couldn't find expected "
               "keys (\"xattr\", \""
            << xattr << "\", \"" << yattr << "\"): " << ps.to_string()
            << std::endl;
  throw;
}

std::vector<TAxis> DetermineOptimalNonUniformBinning(
    TAxis const &UniformAxis,
    std::vector<std::tuple<float, float, double>> const &data,
    size_t NPerBin = 10000 /*1% poisson error*/, double minNU = 0xd34db33f,
    double maxNU = 0xd34db33f, double minNuBinWidth = 0);
