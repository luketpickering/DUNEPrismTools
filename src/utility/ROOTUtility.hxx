#ifndef ROOTUTILITY_HXX_SEEN
#define ROOTUTILITY_HXX_SEEN

#include "StringParserUtility.hxx"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TMatrixD.h"
#include "TRegexp.h"
#include "TTree.h"

#ifdef USE_EIGEN
#include "Eigen/Core"
#endif

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

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

template <class TH>
inline std::unique_ptr<TH> GetHistogram_uptr(TFile *f,
                                             std::string const &fhname) {
  TH *fh = dynamic_cast<TH *>(f->Get(fhname.c_str()));

  if (!fh) {
    std::cout << "[ERROR]: Couldn't get TH: " << fhname
              << " from input file: " << f->GetName() << std::endl;
    exit(1);
  }

  std::unique_ptr<TH> inpH =
      std::unique_ptr<TH>(static_cast<TH *>(fh->Clone()));

  inpH->SetDirectory(nullptr);
  return inpH;
}

template <class TH>
inline std::unique_ptr<TH> GetHistogram_uptr(std::string const &fname,
                                             std::string const &hname) {
  TDirectory *ogDir = gDirectory;

  TFile *inpF = CheckOpenFile(fname);

  std::unique_ptr<TH> h = GetHistogram_uptr<TH>(inpF, hname);

  inpF->Close();
  delete inpF;

  if (ogDir) {
    ogDir->cd();
  }

  return h;
}

template <class TH>
inline std::vector<std::unique_ptr<TH>>
GetHistograms_uptr(std::string const &fname, std::string const &hnamePattern) {
  std::vector<std::unique_ptr<TH>> histos;

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
      histos.push_back(GetHistogram_uptr<TH>(inpF, obj->GetName()));
    }
  }

  inpF->Close();
  delete inpF;

  return histos;
}

std::vector<std::pair<std::pair<double, double>, TH1D *>>
SplitTH2D(TH2D const *t2, bool AlongY,
          double min = -std::numeric_limits<double>::max(),
          double max = std::numeric_limits<double>::max());

std::vector<std::pair<double, TH1D *>>
InterpolateSplitTH2D(TH2D *t2, bool AlongY, std::vector<double> Vals);

std::pair<Int_t, Int_t>
GetProjectionBinRange(std::pair<double, double> ValRange, TAxis *axis);

std::vector<TH1D *> MergeSplitTH2D(TH2D *t2, bool AlongY,
                                   std::vector<std::pair<double, double>> Vals);

std::vector<std::unique_ptr<TH1D>>
MergeSplitTH2D(std::unique_ptr<TH2D> &t2, bool AlongY,
               std::vector<std::pair<double, double>> Vals);

TH2D *SliceNormTH2D(TH2D *t2, bool AlongY);

bool CheckTTreeHasBranch(TTree *tree, std::string const &BranchName);

template <class TH>
void SumHistograms(TH1D *summed, double const *coeffs,
                   std::vector<TH *> const &histos) {
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

template <class TH>
void SumHistograms(TH1D *summed, double const *coeffs,
                   std::vector<std::unique_ptr<TH>> const &histos) {
  size_t nf = histos.size();
  summed->Reset();
  for (size_t f_it = 0; f_it < nf; ++f_it) {
    summed->Add(histos[f_it].get(), coeffs[f_it]);

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

void FindTH1Peaks(TH1D const *flux, int &left, int &right, int n);

double FindHistogramPeak(TH1D *hist, double resolution,
                         std::string const &WriteName);

std::vector<double> Getstdvector(TH2 const *);
std::vector<double> Getstdvector(TH1 const *);

template <typename T>
void Mergestdvector(std::vector<T> &a, std::vector<T> const &b) {
  for (T const &b_i : b) {
    a.push_back(b_i);
  }
}

template <typename THN>
std::vector<double> Getstdvector(std::vector<std::unique_ptr<THN>> const &rhv) {
  std::vector<double> ev;

  for (std::unique_ptr<THN> const &rh : rhv) {
    Mergestdvector(ev, Getstdvector(rh.get()));
  }

  return ev;
}

#ifdef USE_EIGEN
Eigen::MatrixXd GetEigenMatrix(TMatrixD const *);
Eigen::VectorXd GetEigenFlatVector(std::vector<double> const &);
TMatrixD GetTMatrixD(Eigen::MatrixXd const &);
#endif

#endif
