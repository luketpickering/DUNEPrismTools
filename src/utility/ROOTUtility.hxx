#ifndef ROOTUTILITY_HXX_SEEN
#define ROOTUTILITY_HXX_SEEN

#include "StringParserUtility.hxx"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TTree.h"

#include <string>
#include <iostream>
#include <vector>
#include <limits>

TFile *CheckOpenFile(std::string const &fname, char const *opts = "");

template <class T>
inline TChain *OpenTChainWithFileList(std::string const &tname,
                                      std::string const &flist, T &NFiles) {
  TChain *chain = new TChain(tname.c_str());
  std::vector<std::string> fs = ParseToVect<std::string>(flist, ",");
  NFiles = 0;
  for (std::string const &fdesc : fs) {
    std::cout << "[INFO]: Reading TChain(" << tname
              << ") files from descriptor: " << fdesc << std::endl;
    NFiles += chain->Add(fdesc.c_str());
  }
  if (!NFiles || !chain->GetEntries()) {
    std::cout << "[WARN]: Failed to add any files or entries to TChain ("
              << tname << ") from descriptor = \"" << flist << "\"."
              << std::endl;
    delete chain;
    chain = nullptr;
  }
  return chain;
}

TChain *OpenTChainWithFileList(std::string const &tname,
                                      std::string const &flist);

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

std::vector<std::pair<double, TH1D *> > SplitTH2D(
    TH2D *t2, bool AlongY, double min = -std::numeric_limits<double>::max(),
    double max = std::numeric_limits<double>::max());

std::vector<std::pair<double, TH1D *> > InterpolateSplitTH2D(
    TH2D *t2, bool AlongY, std::vector<double> Vals);

std::pair<Int_t, Int_t> GetProjectionBinRange(
    std::pair<double, double> ValRange, TAxis *axis);

std::vector<TH1D *> MergeSplitTH2D(
    TH2D *t2, bool AlongY, std::vector<std::pair<double, double> > Vals);

bool CheckTTreeHasBranch(TTree *tree, std::string const &BranchName);

void SumHistograms(TH1D *summed, double *coeffs,
                          std::vector<TH1D *> const &histos);

double FindHistogramPeak(TH1D *hist, double resolution,
                                std::string const &WriteName);

#endif
