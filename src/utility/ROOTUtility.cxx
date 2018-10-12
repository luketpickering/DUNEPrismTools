#include "ROOTUtility.hxx"

#include "TGraph.h"
#include "TObjArray.h"

#include <numeric>

TFile *CheckOpenFile(std::string const &fname, char const *opts) {
  TFile *inpF = new TFile(fname.c_str(), opts);
  if (!inpF || !inpF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open input file: " << fname << std::endl;
    exit(1);
  }
  return inpF;
}

TChain *OpenTChainWithFileList(std::string const &tname,
                               std::string const &flist) {
  Int_t dummy = 0;
  return OpenTChainWithFileList(tname, flist, dummy);
}

std::vector<std::pair<std::pair<double, double>, TH1D *>>
SplitTH2D(TH2D const *t2, bool AlongY, double min, double max) {
  std::vector<std::pair<std::pair<double, double>, TH1D *>> split;

  for (Int_t bi_it = 1;
       bi_it < (AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetNbins() + 1;
       ++bi_it) {
    if (((AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetBinCenter(bi_it) <
         min) ||
        (AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetBinCenter(bi_it) > max) {
      continue;
    }

    split.push_back(std::make_pair(
        std::make_pair(
            (AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetBinLowEdge(bi_it),
            (AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetBinUpEdge(bi_it)),
        (AlongY ? t2->ProjectionX(to_str(bi_it).c_str(), bi_it, bi_it)
                : t2->ProjectionY(to_str(bi_it).c_str(), bi_it, bi_it))));
    split.back().second->SetDirectory(NULL);
  }

  return split;
}

std::vector<std::pair<double, TH1D *>>
InterpolateSplitTH2D(TH2D *t2, bool AlongY, std::vector<double> Vals) {
  std::vector<std::pair<double, TH1D *>> split;

  for (double v : Vals) {
    TH1D *dummyProj = (AlongY ? t2->ProjectionX(to_str(v).c_str(), 1, 1)
                              : t2->ProjectionY(to_str(v).c_str(), 1, 1));
    dummyProj->Reset();

    for (Int_t bi_it = 1;
         bi_it < (AlongY ? t2->GetXaxis() : t2->GetYaxis())->GetNbins() + 1;
         ++bi_it) {
      dummyProj->SetBinContent(
          bi_it,
          t2->Interpolate(
              (AlongY ? dummyProj->GetXaxis()->GetBinCenter(bi_it) : v),
              (AlongY ? v : dummyProj->GetXaxis()->GetBinCenter(bi_it))));
      Int_t OBin = (AlongY ? t2->GetYaxis() : t2->GetXaxis())->FindFixBin(v);

      if (AlongY) {
        dummyProj->SetBinError(bi_it, t2->GetBinError(bi_it, OBin));
      } else {
        dummyProj->SetBinError(bi_it, t2->GetBinError(OBin, bi_it));
      }
    }

    split.push_back(std::make_pair(v, dummyProj));
    split.back().second->SetDirectory(NULL);
  }

  return split;
}

std::pair<Int_t, Int_t>
GetProjectionBinRange(std::pair<double, double> ValRange, TAxis *axis) {
  Int_t low_bin = axis->FindFixBin(ValRange.first);
  if (fabs(axis->GetBinUpEdge(low_bin) - ValRange.first) < 1E-5) {
    // std::cout << "[INFO]: Axis-selected bin = " << low_bin
    //           << ", with high edge = " << axis->GetBinUpEdge(low_bin)
    //           << ", starting range was requested as " << ValRange.first
    //           << ". Up-shifting low bin." << std::endl;
    low_bin += 1;
  }
  Int_t high_bin = axis->FindFixBin(ValRange.second);
  if (fabs(axis->GetBinLowEdge(high_bin) - ValRange.second) < 1E-5) {
    // std::cout << "[INFO]: Axis-selected bin = " << high_bin
    //           << ", with low edge = " << axis->GetBinLowEdge(high_bin)
    //           << ", ending range was requested as " << ValRange.second
    //           << ". Down-shifting high bin." << std::endl;
    high_bin -= 1;
  }

  if (fabs(axis->GetBinLowEdge(low_bin) - ValRange.first) > 1E-5) {
    std::cout << "[ERROR]: Chose first projection bin = " << low_bin
              << ", with low edge = " << axis->GetBinLowEdge(low_bin)
              << ", but the starting range was requested as " << ValRange.first
              << std::endl;
    exit(1);
  }
  if (fabs(axis->GetBinUpEdge(high_bin) - ValRange.second) > 1E-5) {
    std::cout << "[ERROR]: Chose last projection bin = " << high_bin
              << ", with up edge = " << axis->GetBinLowEdge(high_bin)
              << ", but the ending range was requested as " << ValRange.second
              << std::endl;
    exit(1);
  }

  if (low_bin == 0) {
    std::cout << "[ERROR]: Low bin is underflow bin, 2D flux is not adequate "
                 "for this splitting scheme"
              << std::endl;
    exit(1);
  }
  if (high_bin == (axis->GetNbins() + 1)) {
    std::cout << "[ERROR]: High bin is overflow bin, 2D flux is not adequate "
                 "for this splitting scheme"
              << std::endl;
    exit(1);
  }
  return std::make_pair(low_bin, high_bin);
}

std::vector<TH1D *>
MergeSplitTH2D(TH2D *t2, bool AlongY,
               std::vector<std::pair<double, double>> Vals) {
  std::vector<TH1D *> split;

  for (std::pair<double, double> v : Vals) {
    std::pair<Int_t, Int_t> binr =
        GetProjectionBinRange(v, (AlongY ? t2->GetYaxis() : t2->GetXaxis()));

    split.push_back(
        (AlongY ? t2->ProjectionX(
                      (to_str(v.first) + "_to_" + to_str(v.second)).c_str(),
                      binr.first, binr.second)
                : t2->ProjectionY(
                      (to_str(v.first) + "_to_" + to_str(v.second)).c_str(),
                      binr.first, binr.second)));
    split.back()->Scale(1.0 / double(binr.second - binr.first + 1));
    split.back()->SetDirectory(NULL);
  }

  return split;
}

std::vector<std::unique_ptr<TH1D>>
MergeSplitTH2D(std::unique_ptr<TH2D> &t2, bool AlongY,
               std::vector<std::pair<double, double>> Vals) {
  std::vector<TH1D *> split = MergeSplitTH2D(t2.get(), AlongY, Vals);

  std::vector<std::unique_ptr<TH1D>> split_uptr;

  for (TH1D *h : split) {
    split_uptr.emplace_back(h);
  }

  return split_uptr;
}

TH2D *SliceNormTH2D(TH2D *t2, bool AlongY) {
  TH2D *cpy = static_cast<TH2D *>(t2->Clone(
      (std::string(t2->GetName()) + (AlongY ? "_normYSlice" : "_normXSlice"))
          .c_str()));

  for (Int_t sl_bin_it = 1;
       sl_bin_it < (AlongY ? t2->GetYaxis() : t2->GetXaxis())->GetNbins();
       ++sl_bin_it) {

    double total = 0;

    for (Int_t bin_it = 1;
         bin_it < (AlongY ? t2->GetXaxis() : t2->GetYaxis())->GetNbins();
         ++bin_it) {
      total += (AlongY ? t2->GetBinContent(bin_it + 1, sl_bin_it + 1)
                       : t2->GetBinContent(sl_bin_it + 1, bin_it + 1));
    }

    for (Int_t bin_it = 1;
         bin_it < (AlongY ? t2->GetXaxis() : t2->GetYaxis())->GetNbins();
         ++bin_it) {
      (AlongY ? cpy->SetBinContent(
                    bin_it + 1, sl_bin_it + 1,
                    t2->GetBinContent(bin_it + 1, sl_bin_it + 1) / total)
              : cpy->SetBinContent(
                    sl_bin_it + 1, bin_it + 1,
                    t2->GetBinContent(sl_bin_it + 1, bin_it + 1) / total));
    }
  }
  return cpy;
}

bool CheckTTreeHasBranch(TTree *tree, std::string const &BranchName) {
  TObjArray *branches = tree->GetListOfBranches();

  for (Int_t b_it = 0; b_it < branches->GetEntries(); ++b_it) {
    if (BranchName == branches->At(b_it)->GetName()) {
      return true;
    }
  }
  return false;
}

double FindHistogramPeak(TH1D *hist, double resolution,
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

void FindTH1Peaks(TH1D const *flux, int &left, int &right, int n) {
  // std::cout << "[INFO] Looking for peaks..." << std::endl;

  std::unique_ptr<TH1D> temp = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(flux->Clone("peakfindingtemp")));
  temp->SetDirectory(nullptr);
  temp->Smooth(10);

  double threshold = (temp->Integral()) / (5 * (temp->GetNbinsX()));

  // std::cout << "[INFO] Peak threshold " << threshold << std::endl;

  int nfound = 0;
  double content[3] = {0};

  for (int bin_ind = temp->GetNbinsX(); bin_ind > 0 && nfound < n; bin_ind--) {
    content[2] = temp->GetBinContent(bin_ind - 1);
    if ((content[0] < content[1]) && (content[1] > content[2]) &&
        (content[1] > threshold)) {
      if (nfound == 0) {
        right = bin_ind;
      }
      if (nfound == (n - 1)) {
        left = bin_ind;
      }
      nfound++;
      // std::cout << "[INFO] found a peak of height " << content[1] << " at bin
      // "
      // << bin_ind << std::endl;
    }
    content[0] = content[1];
    content[1] = content[2];
  }
}

std::vector<double> Getstdvector(TH2 const *rh) {
  std::vector<double> ev;

  for (Int_t y_it = 0; y_it < rh->GetYaxis()->GetNbins(); ++y_it) {
    for (Int_t x_it = 0; x_it < rh->GetXaxis()->GetNbins(); ++x_it) {
#ifdef DEBUG_GETFLATVECTOR
      if (rh->GetBinContent(x_it + 1, y_it + 1) &&
          !std::isnormal(rh->GetBinContent(x_it + 1, y_it + 1))) {
        std::cout << "[ERROR]: Got bad bin entry at (" << x_it << ", " << y_it
                  << ") = " << rh->GetBinContent(x_it + 1, y_it + 1)
                  << std::endl;
        throw;
      }
#endif
      ev.push_back(rh->GetBinContent(x_it + 1, y_it + 1));
    }
  }

  return ev;
}

std::vector<double> Getstdvector(TH1 const *rh) {
  Int_t dim = rh->GetDimension();
  if (dim == 1) {
    std::vector<double> ev;
    for (Int_t x_it = 0; x_it < rh->GetXaxis()->GetNbins(); ++x_it) {
#ifdef DEBUG_GETFLATVECTOR
      if (rh->GetBinContent(x_it + 1) &&
          !std::isnormal(rh->GetBinContent(x_it + 1))) {
        std::cout << "[ERROR]: Got bad bin entry at (" << x_it
                  << ") = " << rh->GetBinContent(x_it + 1) << std::endl;
        throw;
      }
#endif
      ev.push_back(rh->GetBinContent(x_it + 1));
    }
    return ev;
  } else if (dim == 2) {
    return Getstdvector(static_cast<TH2 const *>(rh));
  }
  std::cout << "[ERROR]: GetEigenFlatVector cannot handle THND where N = "
            << dim << std::endl;
  throw;
}

#ifdef USE_EIGEN

Eigen::MatrixXd GetEigenMatrix(TMatrixD const *rm) {
  Eigen::MatrixXd em(rm->GetNrows(), rm->GetNcols());

  for (Int_t ir = 0; ir < rm->GetNrows(); ++ir) {
    for (Int_t ic = 0; ic < rm->GetNcols(); ++ic) {
      em(ir, ic) = (*rm)[ir][ic];
    }
  }

  return em;
}

Eigen::VectorXd GetEigenFlatVector(std::vector<double> const &v) {
  Eigen::VectorXd ev(v.size());
  size_t idx = 0;
  size_t N = v.size();
  for (size_t i = 0; i < N; ++i) {
    ev(idx++) = v[i];
  }
  return ev;
}

std::unique_ptr<TMatrixD> GetTMatrixD(Eigen::MatrixXd const &em) {
  std::unique_ptr<TMatrixD> rm(new TMatrixD(em.rows(), em.cols()));

  for (Int_t ir = 0; ir < rm->GetNrows(); ++ir) {
    for (Int_t ic = 0; ic < rm->GetNcols(); ++ic) {
      (*rm)[ir][ic] = em(ir, ic);
    }
  }

  return rm;
}

#endif
