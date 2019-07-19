#ifndef ROOTUTILITY_HXX_SEEN
#define ROOTUTILITY_HXX_SEEN

#include "StringParserUtility.hxx"

#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TMatrixD.h"
#include "TRegexp.h"
#include "TTree.h"

#ifdef USE_EIGEN
#include "Eigen/Core"
#endif

#include "fhiclcpp/ParameterSet.h"

#include <cstdio>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
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
inline TH *GetHistogram(TFile *f, std::string const &fhname,
                        bool exit_on_fail = true) {
  TH *inpH = dynamic_cast<TH *>(f->Get(fhname.c_str()));

  if (!inpH) {
    if (exit_on_fail) {
      std::cout << "[ERROR]: Couldn't get TH: " << fhname
                << " from input file: " << f->GetName() << std::endl;
      exit(1);
    } else {
      return nullptr;
    }
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
inline std::unique_ptr<TH>
GetHistogram_uptr(TFile *f, std::string const &fhname, bool no_throw = false) {
  TH *fh = dynamic_cast<TH *>(f->Get(fhname.c_str()));

  if (!fh) {
    if (no_throw) {
      return std::unique_ptr<TH>(nullptr);
    }
    std::cout << "[ERROR]: Couldn't get TH: " << fhname
              << " from input file: " << f->GetName() << std::endl;
    throw;
  }

  std::unique_ptr<TH> inpH =
      std::unique_ptr<TH>(static_cast<TH *>(fh->Clone()));

  inpH->SetDirectory(nullptr);
  return inpH;
}

template <class TH>
inline std::unique_ptr<TH> GetHistogram_uptr(std::string const &fname,
                                             std::string const &hname,
                                             bool no_throw = false) {
  TDirectory *ogDir = gDirectory;

  TFile *inpF = CheckOpenFile(fname);

  std::unique_ptr<TH> h = GetHistogram_uptr<TH>(inpF, hname, no_throw);

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

template <class RType>
inline std::unique_ptr<RType>
GetROOTObject_uptr(TFile *f, std::string const &fhname, bool no_throw = false) {
  RType *fh = dynamic_cast<RType *>(f->Get(fhname.c_str()));

  if (!fh) {
    if (no_throw) {
      return std::unique_ptr<RType>(nullptr);
    }
    std::cout << "[ERROR]: Couldn't get RType: " << fhname
              << " from input file: " << f->GetName() << std::endl;
    exit(1);
  }

  std::unique_ptr<RType> inpH =
      std::unique_ptr<RType>(static_cast<RType *>(fh->Clone()));

  return inpH;
}

template <class RType>
inline std::unique_ptr<RType> GetROOTObject_uptr(std::string const &fname,
                                                 std::string const &hname,
                                                 bool no_throw = false) {
  TDirectory *ogDir = gDirectory;

  TFile *inpF = CheckOpenFile(fname);

  std::unique_ptr<RType> h = GetROOTObject_uptr<RType>(inpF, hname, no_throw);

  inpF->Close();
  delete inpF;

  if (ogDir) {
    ogDir->cd();
  }

  return h;
}

std::vector<std::pair<std::pair<double, double>, TH1D *>>
SplitTH2D(TH2D const *t2, bool AlongY,
          double min = -std::numeric_limits<double>::max(),
          double max = std::numeric_limits<double>::max());

std::vector<std::pair<double, TH1D *>>
InterpolateSplitTH2D(TH2D *t2, bool AlongY, std::vector<double> Vals);

std::pair<Int_t, Int_t>
GetProjectionBinRange(std::pair<double, double> ValRange, TAxis *axis);

std::vector<TH1D *>
MergeSplitTH2D(TH2D *t2, bool AlongY,
               std::vector<std::pair<double, double>> const &Vals);

std::vector<std::unique_ptr<TH1D>>
MergeSplitTH2D(std::unique_ptr<TH2D> &t2, bool AlongY,
               std::vector<std::pair<double, double>> const &Vals);

std::unique_ptr<TH2D>
ReMergeSplitTH2D(std::vector<std::unique_ptr<TH1D>> const &,
                 std::vector<std::pair<double, double>> const &Vals,
                 std::string const &name, std::string const &title = "");

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

std::vector<double> Getstdvectorerror(TH2 const *);
std::vector<double> Getstdvectorerror(TH1 const *);

template <typename T>
inline void Mergestdvector(std::vector<T> &a, std::vector<T> const &b) {
  for (T const &b_i : b) {
    a.push_back(b_i);
  }
}

template <typename T>
inline void Mergestdvector(std::vector<T> &a, std::vector<T> &&b) {
  for (T &b_i : b) {
    a.push_back(std::move(b_i));
  }
}

template <typename THN>
inline std::vector<double>
Getstdvector(std::vector<std::unique_ptr<THN>> const &rhv) {
  std::vector<double> ev;

  for (std::unique_ptr<THN> const &rh : rhv) {
    Mergestdvector(ev, Getstdvector(rh.get()));
  }

  return ev;
}

template <typename THN>
inline std::vector<double>
Getstdvectorerror(std::vector<std::unique_ptr<THN>> const &rhv) {
  std::vector<double> ev;

  for (std::unique_ptr<THN> const &rh : rhv) {
    Mergestdvector(ev, Getstdvectorerror(rh.get()));
  }

  return ev;
}

size_t FillHistFromstdvector(TH2 *, std::vector<double> const &,
                             size_t offset = 0,
                             std::vector<double> const &error = {});
size_t FillHistFromstdvector(TH1 *, std::vector<double> const &,
                             size_t offset = 0,
                             std::vector<double> const &error = {});
template <typename THN>
inline size_t FillHistFromstdvector(std::vector<std::unique_ptr<THN>> &rhv,
                                    std::vector<double> const &vals,
                                    size_t offset = 0,
                                    std::vector<double> const &error = {}) {

  for (std::unique_ptr<THN> const &rh : rhv) {
    offset = FillHistFromstdvector(rh.get(), vals, offset, error);
  }
  return offset;
}

template <typename THN>
inline std::vector<std::unique_ptr<THN>>
CloneHistVector(std::vector<std::unique_ptr<THN>> const &rhv,
                std::string const &post_fix = "") {
  std::vector<std::unique_ptr<THN>> clonev;
  for (std::unique_ptr<THN> const &rh : rhv) {
    clonev.emplace_back(static_cast<THN *>(
        rh->Clone((std::string(rh->GetName()) + post_fix).c_str())));
    clonev.back()->SetDirectory(nullptr);
  }
  return clonev;
}

#ifdef USE_EIGEN

size_t
FillHistFromEigenVector(TH2 *, Eigen::VectorXd const &, size_t offset = 0,
                        Eigen::VectorXd const &error = Eigen::VectorXd());
size_t
FillHistFromEigenVector(TH1 *, Eigen::VectorXd const &, size_t offset = 0,
                        Eigen::VectorXd const &error = Eigen::VectorXd());
template <typename THN>
inline size_t
FillHistFromEigenVector(std::vector<std::unique_ptr<THN>> &rhv,
                        Eigen::VectorXd const &vals, size_t offset = 0,
                        Eigen::VectorXd const &error = Eigen::VectorXd()) {

  for (std::unique_ptr<THN> const &rh : rhv) {
    offset = FillHistFromEigenVector(rh.get(), vals, offset, error);
  }
  return offset;
}

Eigen::MatrixXd GetEigenMatrix(TMatrixD const *);
inline Eigen::MatrixXd GetEigenMatrix(TMatrixD const &m) {
  return GetEigenMatrix(&m);
}
Eigen::MatrixXd GetEigenMatrix(TH2 const *);
Eigen::VectorXd GetEigenFlatVector(std::vector<double> const &);
std::unique_ptr<TMatrixD> GetTMatrixD(Eigen::MatrixXd const &);
#endif

template <size_t n>
inline std::array<double, n + 1>
GetPolyFitCoeffs(std::vector<double> const &xvals,
                 std::vector<double> const &yvals, bool LimitNDOF = true) {

  size_t nd = std::min(xvals.size(), yvals.size());
  TGraph g(nd);
  for (size_t i = 0; i < nd; ++i) {
    g.SetPoint(i, xvals[i], yvals[i]);
  }

  size_t dof =
      std::min(size_t(9), LimitNDOF ? std::min(xvals.size() - 2, n) : n);

  static char polstr[10];
  sprintf(polstr, "pol%d", int(dof));

  TFitResultPtr r = g.Fit(polstr, "Q0S");

  std::array<double, n + 1> rtn;
  std::fill_n(rtn.begin(), n + 1, 0);
  if (!r->IsValid()) {
    std::cout << "[FIT FAILURE]: " << polstr
              << " fit failed, returning empty result." << std::endl;
    return rtn;
  }
  for (size_t i = 0; i < (dof + 1); ++i) {
    rtn[i] = r->Parameter(i);
  }

  return rtn;
}

template <size_t n>
inline std::array<double, n + 1>
GetPolyFitCoeffs(std::vector<std::pair<double, double>> const &xyvals,
                 bool LimitNDOF = true) {

  size_t nd = xyvals.size();
  TGraph g(nd);
  for (size_t i = 0; i < nd; ++i) {
    g.SetPoint(i, xyvals[i].first, xyvals[i].second);
  }

  size_t dof =
      std::min(size_t(9), LimitNDOF ? std::min(xyvals.size() - 2, n) : n);

  static char polstr[10];
  sprintf(polstr, "pol%d", int(dof));

  TFitResultPtr r = g.Fit(polstr, "Q0S");

  std::array<double, n + 1> rtn;
  std::fill_n(rtn.begin(), n + 1, 0);
  if (!r->IsValid()) {
    std::cout << "[FIT FAILURE]: " << polstr
              << " fit failed, returning empty result." << std::endl;
    return rtn;
  }
  for (size_t i = 0; i < (dof + 1); ++i) {
    rtn[i] = r->Parameter(i);
  }

  return rtn;
}

template <size_t n>
inline std::array<double, n + 1>
GetPolyFitCoeffs(std::vector<std::tuple<double, double, double>> const &xyevals,
                 bool LimitNDOF = true) {

  size_t nd = xyevals.size();
  TGraphErrors g(nd);
  for (size_t i = 0; i < nd; ++i) {
    g.SetPoint(i, std::get<0>(xyevals[i]), std::get<1>(xyevals[i]));
    g.SetPointError(i, 0, std::get<2>(xyevals[i]));
  }

  size_t dof =
      std::min(size_t(9), LimitNDOF ? std::min(xyevals.size() - 2, n) : n);

  static char polstr[10];
  sprintf(polstr, "pol%d", int(dof));

  TFitResultPtr r = g.Fit(polstr, "Q0S");

  std::array<double, n + 1> rtn;
  std::fill_n(rtn.begin(), n + 1, 0);

  if (!r->IsValid()) {
    std::cout << "[FIT FAILURE]: " << polstr
              << " fit failed, returning empty result." << std::endl;
    return rtn;
  }

  for (size_t i = 0; i < (dof + 1); ++i) {
    rtn[i] = r->Parameter(i);
  }

  return rtn;
}

std::unique_ptr<TH1> THToF(std::unique_ptr<TH1> &hin);

inline double GetTH1RatioError(TH1 *num, TH1 *denom, Int_t bi_it) {
  double num_v = num->GetBinContent(bi_it);
  double num_e = num->GetBinError(bi_it);

  double denom_v = denom->GetBinContent(bi_it);
  double denom_e = denom->GetBinError(bi_it);

  return ((num_v / denom_v) * sqrt(pow(num_e, 2) / pow(num_v, 2) +
                                   pow(denom_e, 2) / pow(denom_v, 2)));
}

#endif
