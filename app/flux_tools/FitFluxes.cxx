#include "FluxFitResultsTreeReader.hxx"
#include "OscillationParametersTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitter.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMinuit.h"

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

std::string InputTargetFluxFile, InputTargetFluxName;
std::string InputFluxFile, InputFluxName;
std::string OutputFile, OutputDirectory="";
std::string InpCoeffFile="", InpCoeffDir="";

bool UPDATEOutputFile = false;

bool IsGauss = false;
double GaussC, GaussW;
double TargetPeakNorm;

double CoeffLimit = 30;
double RegFactor = 0xdeadbeef;

int MaxMINUITCalls = 50000;

bool FitBetweenFoundPeaks = false;

double FitBetween_low = 0xdeadbeef, FitBetween_high = 0xdeadbeef;
int FitBinLow, FitBinHigh;

std::vector<std::pair<double, double> > XRanges;
std::vector<bool> ApplyReg;

std::vector<TH1D *> Fluxes;
TH1D *SummedFlux;
TH1D *TargetFlux;
TH1D *TargetFlux_orignorm;
TH1D *OscFlux;
TH1D *OscFlux_orignorm;

TF1 *TargetGauss;

enum OutOfRangeModeEnum {
  kIgnore = 0,
  kZero,
  kExponentialDecay,
  kGaussianDecay
};
/// If using a FitBetween mode:
/// 0: Ignore all bins outside range
/// 1: Try to force bins to 0
/// 2: Exponential decay from target flux at closest kept bin.
int OutOfRangeMode = kZero;
enum OutOfRangeSideEnum { kBoth = 0, kLeft, kRight };

/// If using an OOR mode:
/// 0: Include both out of ranges
/// 1: Only include out of range to the left of the range
/// 2: Only include out of range to the right of the range
int OutOfRangeSide = kBoth;
double ExpDecayRate = 3;
double OORFactor = 1;

bool UseNuPrismChi2 = false;

int MergeOAABins = 0, MergeENuBins = 0;

void FindPeaks(TH1D *OscFlux, int &left, int &right, int n = 3) {
  std::cout << "[INFO] Looking for peaks..." << std::endl;

  TH1D *temp = new TH1D();

  OscFlux->Copy(*temp);
  temp->Smooth(10);

  double threshold = (temp->Integral()) / (5 * (temp->GetNbinsX()));

  // double threshold = 1.e-16;

  std::cout << "[INFO] Peak threshold " << threshold << std::endl;

  int nfound = 0;
  double content[3] = {0};

  for (int bin_ind = temp->GetNbinsX(); bin_ind > 0 && nfound < n; bin_ind--) {
    content[2] = temp->GetBinContent(bin_ind - 1);
    if (content[0] < content[1] && content[1] > content[2] &&
        content[1] > threshold) {
      if (nfound == 0) right = bin_ind;
      if (nfound == n - 1) left = bin_ind;
      nfound++;
      std::cout << "[INFO] found a peak of height " << content[1] << " at bin "
                << bin_ind << std::endl;
    }
    content[0] = content[1];
    content[1] = content[2];
  }
}

void BuildTargetFlux(TH1D *OscFlux) {
  TargetFlux = static_cast<TH1D *>(OscFlux->Clone());
  TargetFlux->SetDirectory(NULL);

  for (Int_t bi_it = 1; bi_it < FitBinLow; ++bi_it) {
    if ((OutOfRangeMode != kIgnore) &&
        (OutOfRangeSide == kBoth || OutOfRangeSide == kLeft)) {
      double target = 0;
      if (OutOfRangeMode == kExponentialDecay) {
        double enu_first_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(FitBinLow);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_first_counted_bin = TargetFlux->GetBinContent(FitBinLow);
        double enu_bottom_bin = TargetFlux->GetXaxis()->GetBinCenter(1);
        double sigma5_range = enu_first_counted_bin - enu_bottom_bin;
        target =
            content_first_counted_bin *
            exp(-ExpDecayRate * (enu_first_counted_bin - enu) / sigma5_range);
      }
      if (OutOfRangeMode == kGaussianDecay) {
        double enu_first_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(FitBinLow);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_first_counted_bin = TargetFlux->GetBinContent(FitBinLow);
        double enu_bottom_bin = TargetFlux->GetXaxis()->GetBinCenter(1);
        double sigma5_range = enu_first_counted_bin - enu_bottom_bin;
        target =
            content_first_counted_bin *
            exp(-ExpDecayRate * (enu_first_counted_bin - enu) *
                (enu_first_counted_bin - enu) / (sigma5_range * sigma5_range));
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }

  for (Int_t bi_it = FitBinLow; bi_it < (FitBinHigh + 1); ++bi_it) {
    TargetFlux->SetBinContent(bi_it, OscFlux->GetBinContent(bi_it));
    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }

  for (Int_t bi_it = (FitBinHigh + 1);
       bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
    if ((OutOfRangeMode != kIgnore) &&
        (OutOfRangeSide == kBoth || OutOfRangeSide == kRight)) {
      double target = 0;
      if (OutOfRangeMode == kExponentialDecay) {
        double enu_last_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(FitBinHigh);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_last_counted_bin = TargetFlux->GetBinContent(FitBinHigh);
        double enu_top_bin = TargetFlux->GetXaxis()->GetBinCenter(
            TargetFlux->GetXaxis()->GetNbins());
        double sigma5_range = enu_top_bin - enu_last_counted_bin;
        target =
            content_last_counted_bin *
            exp(-ExpDecayRate * (enu - enu_last_counted_bin) / sigma5_range);
      }
      if (OutOfRangeMode == kGaussianDecay) {
        double enu_last_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(FitBinHigh);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_last_counted_bin = TargetFlux->GetBinContent(FitBinHigh);
        double enu_top_bin = TargetFlux->GetXaxis()->GetBinCenter(
            TargetFlux->GetXaxis()->GetNbins());
        double sigma5_range = enu_top_bin - enu_last_counted_bin;
        target =
            content_last_counted_bin *
            exp(-ExpDecayRate * (enu - enu_last_counted_bin) *
                (enu - enu_last_counted_bin) / (sigma5_range * sigma5_range));
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }
}

size_t NCalls = 0;
double Chi2_last, OOR_last, reg_last;

void TargetSumChi2(int &nDim, double *gout, double &result, double coeffs[],
                   int flg) {
  SumHistograms(SummedFlux, coeffs, Fluxes);
  NCalls++;

  double sumdiff = 0;
  double OOR = 0;

  if ((OutOfRangeMode != kIgnore)) {
    if ((OutOfRangeSide == kBoth || OutOfRangeSide == kLeft)) {
      for (Int_t bi_it = 1; bi_it < FitBinLow; ++bi_it) {
        double diff, err;
        if (!UseNuPrismChi2) {
          err =
              (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
               SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it));

          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);

        } else {
          err = SummedFlux->GetBinContent(bi_it)
                    ? (0.0001 * SummedFlux->GetBinContent(bi_it) *
                       SummedFlux->GetBinContent(bi_it))
                    : (0.0001 * TargetPeakNorm * TargetPeakNorm);
          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);
        }

        if (!err) {
          continue;
        }

        if (diff && !std::isnormal(diff)) {
          std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                    << " [ORR LowE], targ: " << TargetFlux->GetBinContent(bi_it)
                    << ":" << TargetFlux->GetBinError(bi_it)
                    << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it) << ", err = " << err
                    << ", diff = " << diff << std::endl;
          throw;
        }

        sumdiff += diff;
        OOR += diff;
      }
    }

    if ((OutOfRangeSide == kBoth || OutOfRangeSide == kRight)) {
      for (Int_t bi_it = (FitBinHigh + 1);
           bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
        double diff, err;
        if (!UseNuPrismChi2) {
          err = SummedFlux->GetBinContent(bi_it)
                    ? (0.0001 * SummedFlux->GetBinContent(bi_it) *
                       SummedFlux->GetBinContent(bi_it))
                    : (0.0001 * TargetPeakNorm * TargetPeakNorm);

          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);

        } else {
          err = (0.0001 * SummedFlux->GetBinContent(bi_it) *
                 SummedFlux->GetBinContent(bi_it));
          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);
        }

        if (!err) {
          continue;
        }

        if (diff && !std::isnormal(diff)) {
          std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                    << " [ORR HighE], targ: "
                    << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it)
                    << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it) << ", err = " << err
                    << ", diff = " << diff << std::endl;
          throw;
        }

        sumdiff += diff;
        OOR += diff;
      }
    }
  }

  for (Int_t bi_it = FitBinLow; bi_it < (FitBinHigh + 1); ++bi_it) {
    double diff, err;
    if (!UseNuPrismChi2) {
      err = (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
             SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it));

      diff = ((TargetFlux->GetBinContent(bi_it) -
               SummedFlux->GetBinContent(bi_it)) *
              (TargetFlux->GetBinContent(bi_it) -
               SummedFlux->GetBinContent(bi_it)) /
              err);
    } else {
      err = SummedFlux->GetBinContent(bi_it)
                ? (0.0001 * SummedFlux->GetBinContent(bi_it) *
                   SummedFlux->GetBinContent(bi_it))
                : (0.0001 * TargetPeakNorm * TargetPeakNorm);

      diff = ((TargetFlux->GetBinContent(bi_it) -
               SummedFlux->GetBinContent(bi_it)) *
              (TargetFlux->GetBinContent(bi_it) -
               SummedFlux->GetBinContent(bi_it)) /
              err);
    }

    if (diff && !std::isnormal(diff)) {
      std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it)
                << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it) << ", err = " << err
                << ", diff = " << diff << std::endl;
      throw;
    }

    sumdiff += diff;

    if (sumdiff && !std::isnormal(sumdiff)) {
      std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it)
                << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it) << std::endl;
      throw;
    }
  }

  double reg = 0;
  if (RegFactor != 0xdeadbeef) {
    for (size_t i = 0; i < Fluxes.size(); i++) {
      if (ApplyReg.size() && !ApplyReg[i]) {
        continue;
      }
      reg += pow((coeffs[i] - coeffs[i + 1]) / RegFactor, 2);
    }
  }

  result = sumdiff + reg;

  Chi2_last = result;
  OOR_last = OOR;
  reg_last = reg;

  std::cout << "[INFO]: Target flux chi2: " << result << " (reg = " << reg
            << ", OOR = " << OOR << " )." << std::endl;
}

void TargetSumGauss(int &nDim, double *gout, double &result, double coeffs[],
                    int flg) {
  NCalls++;
  SumHistograms(SummedFlux, coeffs, Fluxes);

  double sumdiff = 0;
  for (Int_t bi_it = FitBinLow; bi_it < FitBinHigh + 1; ++bi_it) {
    double bi_c_E = SummedFlux->GetXaxis()->GetBinCenter(bi_it);

    double GaussEval = TargetGauss->Eval(bi_c_E);
    double SummedBinContent = SummedFlux->GetBinContent(bi_it);

    if (!UseNuPrismChi2) {
      double Uncert = pow(GaussEval / 20.0, 2) + pow(TargetPeakNorm / 30.0, 2);
      sumdiff += (pow(GaussEval - SummedBinContent, 2) / Uncert);
    } else {
      double GUncert =
          pow(0.0001 * GaussEval, 2) + pow(0.00005 * TargetPeakNorm, 2);
      sumdiff += (pow(GaussEval - SummedBinContent, 2) / GUncert);
    }

    if (sumdiff && !std::isnormal(sumdiff)) {
      std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                << ", gauss: " << GaussEval << ":" << TargetPeakNorm
                << ", sum: " << SummedFlux->GetBinContent(bi_it) << ":"
                << SummedFlux->GetBinError(bi_it) << std::endl;
      throw;
    }
  }

  double reg = 0;
  if (RegFactor != 0xdeadbeef) {
    for (size_t i = 0; i < Fluxes.size(); i++) {
      if (ApplyReg.size() && !ApplyReg[i]) {
        continue;
      }
      reg += pow((coeffs[i] - coeffs[i + 1]) / RegFactor, 2);
    }
  }
  sumdiff += reg;

  Chi2_last = result;
  OOR_last = 0;
  reg_last = reg;

  result = sumdiff;
  std::cout << "[INFO]: Gauss chi2: " << sumdiff << " (reg = " << reg << " )."
            << std::endl;
}

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "  Input options:                                               "
         "           \n"
         "\n"
         "\t-f <ROOT file,FluxHist2DName>      : Input 2D flux histogram,"
         " Y bins \n"
         "\t                                     correspond to different "
         "fluxes.\n"
         "\n"
         "\t-MX  <nbins to merge>              : Merge neutrino energy bins "
         "before splitting into\n"
         "\t                                     fluxes."
         "\n\n"
         "\t-M  <OA1>:<OA_W>,<OA2>_<OAN>:<OA_W>,<OAN+1>:<OA_W>,...\n "
         "\t                                   : Merge bins in off axis flux "
         "positions. Each position \n"
         "\t                                     or position range specifies a "
         "slice width. The \n"
         "\t                                     corresponding absolute slice "
         "ranges must match up to\n"
         "\t                                     merge-able Y bin edges from "
         "histogram passed to -f.\n"
         "\t                                     You will be notified if they "
         "don't.\n"
         "\n"
         "\t-A <FitOutput.root[,dirname]>      : Start a new fit from the "
         "results of an \n"
         "\t                                     old fit. Optional dirname "
         "corresponds to \n"
         "\t                                     -d option. (Tip: set -n 0 "
         "to apply previous \n"
         "\t                                     results to new inputs without "
         "running a fit.)\n"
         "\n"
         "  Output options:                                              "
         "           \n"
         "\t-[o|a] <ROOT file>                 : The output root file. "
         "Using -o will\n"
         "\t                                     overwrite a file of the "
         "same name, \n"
         "\t                                     -a will append the fit "
         "result to   \n"
         "\t                                     the file.               "
         "           \n"
         "\n"
         "\t-d <directory name>                : If passed, fit result "
         "will be put  \n"
         "\t                                     into a subdirectory of "
         "the root    \n"
         "\t                                     file.                   "
         "           \n"
         "\n"
         "  Target options:                                              "
         "           \n"
         "\t-t <ROOT file,hist name>           : The histogram of the "
         "target flux to\n"
         "\t                                     fit to. This file should "
         "contain an \n"
         "\t                                     oscillation parameter config "
         "tree generated \n"
         "\t                                     by dp_OscillateFlux.\n"
         "\n"
         "\t-g <mean,width>                    : Use a gaussian target "
         "distribution \n"
         "\t                                     instead of a target "
         "flux shape.    \n"
         "\n"
         "  Fitter options:                                            "
         "           \n"
         "\t-n <MaxCalls=50000>                : The maximum number of "
         "MINUIT       \n"
         "\t                                     evaluations before "
         "giving up the   \n"
         "\t                                     fit."
         "\n"
         "\n"
         "\t-c <CoeffLimit=30>                 : Parameter limits of "
         "flux component \n"
         "\t                                     coefficients.           "
         "           \n"
         "\n"
         "  Figure of merit options:                                     "
         "           \n"
         "\t-C                                 : Use NuPrism tools Chi2.\n\n"
         "\t-rg <regularisation factor>        : Adds neighbouring       "
         "coefficient\n"
         "\t                                     regularisation.         "
         "           \n"
         "\n"
         "\t-l <min val>,<max val>             : Fit between min and max."
         " Outside \n"
         "\t                                     of this range, -m"
         "determines\n"
         "\t                                     behavior.             \n"
         "\n"
         "\t-p                                 : Fit between the first"
         "and third oscillation \n"
         "\t                                     peaks (Only useful with -t) \n"
         "\n"
         "\t-m <0,1,2,3>                       : Out of range behavior.  "
         "          \n"
         "\t                                     0: Ignore out of range "
         "bins.      \n"
         "\t                                     1: Force out of range "
         "bins to 0.  \n"
         "\t                                     2: Exponential decay "
         "outside fit  \n"
         "\t                                        region. Decay rate "
         "is          \n"
         "\t                                        determined by -ed."
         "\n"
         "\t                                     3: Gaussian decay "
         "outside fit  \n"
         "\t                                        region. Decay width "
         "is          \n"
         "\t                                        determined by -ed."
         "\n"
         "\n"
         "\t-ed <decay rate>                   : For -m [2,3], controls "
         "decay rate.    \n"
         "\t                                     Default = 3, larger is "
         "faster     \n"
         "\t                                     decay.\n"
         "\n"
         "\t-of <out of range factor>          : Allow out of range to "
         "contribute less to the FOM\n "
         "\t                                     by this factor.\n"
         "\n"
         "\t-ms <out of range side>            : 0 = Include both low and high "
         "out of range E, \n"
         "\t                                     1 = include low E, 2 = "
         "include "
         "high E.\n"
         "\n"
      << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-g") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for g, expected 2." << std::endl;
        throw;
      }
      GaussC = params[0];
      GaussW = params[1];
      if ((GaussC <= 0) || (GaussW <= 0)) {
        std::cout << "[ERROR]: Recieved " << argv[opt]
                  << " argument for -g, expected \"<GaussMean>,<GaussWidth>\", "
                     "where both values are > 0."
                  << std::endl;
        throw;
      }
      IsGauss = true;
      TargetGauss = new TF1("tFunc", "gaus", 0, 20);
    } else if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for l, expected 2." << std::endl;
        throw;
      }
      FitBetween_low = params[0];
      FitBetween_high = params[1];
    } else if (std::string(argv[opt]) == "-p") {
      FitBetweenFoundPeaks = true;
    } else if (std::string(argv[opt]) == "-t") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -t, expected 2." << std::endl;
        throw;
      }
      InputTargetFluxFile = params[0];
      InputTargetFluxName = params[1];
    } else if (std::string(argv[opt]) == "-A") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      InpCoeffFile = params[0];
      if (params.size() > 1) {
        InpCoeffDir = params[1];
      } else {
        InpCoeffDir = "";
      }
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
      UPDATEOutputFile = false;
    } else if (std::string(argv[opt]) == "-a") {
      OutputFile = argv[++opt];
      UPDATEOutputFile = true;
    } else if (std::string(argv[opt]) == "-d") {
      OutputDirectory = argv[++opt];
    } else if (std::string(argv[opt]) == "-f") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -f, expected 2." << std::endl;
        throw;
      }
      InputFluxFile = params[0];
      InputFluxName = params[1];
    } else if (std::string(argv[opt]) == "-n") {
      MaxMINUITCalls = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-c") {
      CoeffLimit = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-m") {
      OutOfRangeMode = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-rg") {
      RegFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ed") {
      ExpDecayRate = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-of") {
      OORFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ms") {
      OutOfRangeSide = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-C") {
      UseNuPrismChi2 = true;
    } else if (std::string(argv[opt]) == "-MX") {
      MergeENuBins = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-M") {
      XRanges = BuildRangesList(argv[++opt]);
    } else if ((std::string(argv[opt]) == "-?") ||
               std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      throw;
    }
    opt++;
  }
}
int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!OutputFile.length()) {
    std::cout << "[ERROR]: No output file specified." << std::endl;
    throw;
  }

  if((!InputFluxFile.length()) || (!InputFluxName.length())){
    std::cout << "[ERROR]: No input flux file or histogram name specified."
              << std::endl;
    throw;
  }

  if (!IsGauss) { // Set up target flux
    if((!InputTargetFluxFile.length()) || (!InputTargetFluxName.length())){
      std::cout << "[ERROR]: No target flux file or histogram name specified."
                << std::endl;
      throw;
    }

    OscFlux = GetHistogram<TH1D>(InputTargetFluxFile, InputTargetFluxName);

    if (FitBetweenFoundPeaks) {
      FindPeaks(OscFlux, FitBinLow, FitBinHigh);
    } else {
      if (FitBetween_low == 0xdeadbeef) {
        FitBinLow = 1;
        FitBinHigh = OscFlux->GetXaxis()->GetNbins();
      } else {
        FitBinLow = OscFlux->GetXaxis()->FindFixBin(FitBetween_low);
        FitBinHigh = OscFlux->GetXaxis()->FindFixBin(FitBetween_high);
      }
    }

    std::cout << "[INFO]: Fitting between bins: " << FitBinLow << "--"
      << FitBinHigh << std::endl;

    BuildTargetFlux(OscFlux);
  }

    TH2D *Flux2D = GetHistogram<TH2D>(InputFluxFile, InputFluxName);
    if (!Flux2D) {
      std::cout << "[ERROR]: Found no input flux with name: \""
                << InputFluxName << "\" in file: \"" << InputFluxFile << "\"."
                << std::endl;
      throw;
    }

  if (MergeENuBins) {
    Flux2D->RebinX(MergeENuBins);
    Flux2D->Scale(1.0 / double(MergeENuBins));
  }

  if (XRanges.size()) {
    Fluxes = MergeSplitTH2D(Flux2D, true, XRanges);
  } else {
    std::vector<std::pair<std::pair<double,double>, TH1D *> > SplitFlux2D =
      SplitTH2D(Flux2D,true);
    for (size_t i = 0; i < SplitFlux2D.size(); ++i) {
      XRanges.push_back(SplitFlux2D[i].first);
      Fluxes.push_back(SplitFlux2D[i].second);
    }
  }
  for (size_t i = 0; i < XRanges.size(); ++i) {
    std::cout << "[INFO]: Built flux for slice " << i << " between "
              << XRanges[i].first << ", and "
              << XRanges[i].second << " m."
              << std::endl;
    if (i != 0) {
      // Don't regularise across gaps
      if (fabs(XRanges[i].first -
               XRanges[i - 1].second) > 1E-5) {
        ApplyReg.back() = false;
      }
    }
    ApplyReg.push_back(true);
  }

  if (!Fluxes.size()) {
    std::cout << "[ERROR]: Found no input fluxes." << std::endl;
    throw;
  }

  SummedFlux = static_cast<TH1D *>(Fluxes[0]->Clone());
  SummedFlux->SetDirectory(NULL);
  SummedFlux->Reset();

  if (IsGauss) {
    if (FitBetween_low == 0xdeadbeef) {
      FitBinLow = 1;
      FitBinHigh = SummedFlux->GetXaxis()->GetNbins();
    } else {
      FitBinLow = SummedFlux->GetXaxis()->FindFixBin(FitBetween_low);
      FitBinHigh = SummedFlux->GetXaxis()->FindFixBin(FitBetween_high);
    }
  }

  double *coeffs = new double[Fluxes.size()];
  if (InpCoeffFile.size()) {
    if (InpCoeffDir.size() && (InpCoeffDir.back() != '/')) {
      InpCoeffDir += "/";
    }

    std::string treeName = InpCoeffDir + "CoeffTree";

    std::cout << "[INFO]: Reading input previous fit results." << std::endl;
    SliceConfig sc(treeName, InpCoeffFile);

    std::vector< std::pair<double,double> > inpXRanges = sc.GetXRanges();
    std::vector<double> inpCoeffs = sc.GetCoeffs();

    if(inpXRanges.size() != XRanges.size()){
      std::cout << "[ERROR]: Input fit result had " << inpXRanges.size()
        << " off-axis slices, but here we have found " << XRanges.size()
        << ". Input results are incompatible." << std::endl;
      throw;
    }

    for(size_t f_it = 0; f_it < inpXRanges.size(); ++f_it){
      std::cout << "\tFlux window[" << f_it << "] = {"
        << inpXRanges[f_it].first << " -- " << inpXRanges[f_it].second
        << "}, Coeff = " << inpCoeffs[f_it] << std::endl;
      if( (fabs(inpXRanges[f_it].first-XRanges[f_it].first) > 1E-5 ) ||
          (fabs(inpXRanges[f_it].second-XRanges[f_it].second) > 1E-5)   ){
        std::cout << "[ERROR]: Here we found Flux window[" << f_it << "] = {"
          << XRanges[f_it].first << " -- " << XRanges[f_it].second
          << "}. Input results are incompatible." << std::endl;
          throw;
      }
      coeffs[f_it] = inpCoeffs[f_it];
    }

  } else {
    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      coeffs[flux_it] = 0;
    }
  }

  if (!IsGauss) {
    // Rescale the target to a similar size to the fluxes.
    double OnAxisPeak = Fluxes[0]->GetMaximum();
    double TargetMax = TargetFlux->GetMaximum();
    double NDOverFDFitScaleFactor = OnAxisPeak / TargetMax;
    TargetFlux_orignorm = static_cast<TH1D *>(
      TargetFlux->Clone("TargetFlux_orignorm"));
    TargetFlux_orignorm->SetDirectory(nullptr);
    TargetFlux->Scale(NDOverFDFitScaleFactor);

    OscFlux_orignorm = static_cast<TH1D *>(OscFlux->Clone("OscFlux_orignorm"));
    OscFlux_orignorm->SetDirectory(nullptr);
    OscFlux->Scale(NDOverFDFitScaleFactor);
  }

  TFitter *minimizer = nullptr;
  int fitstatus = -1;

  if (MaxMINUITCalls != 0) {
    minimizer = new TFitter(Fluxes.size());

    double targetenu = IsGauss ? GaussC : TargetFlux->GetXaxis()->GetBinCenter(
                                              TargetFlux->GetMaximumBin());

    size_t flux_with_closest_peak_it = -1;
    double flux_peak_dist = std::numeric_limits<double>::max();
    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      double maxenu = Fluxes[flux_it]->GetXaxis()->GetBinCenter(
          Fluxes[flux_it]->GetMaximumBin());
      double peak_difference = fabs(maxenu - targetenu);
      if (peak_difference < flux_peak_dist) {
        flux_peak_dist = peak_difference;
        flux_with_closest_peak_it = flux_it;
      }
    }

    std::cout << "[INFO]: flux_with_closest_peak_it: "
              << flux_with_closest_peak_it
              << ", flux_peak_dist = " << flux_peak_dist << std::endl;

    TargetPeakNorm = Fluxes[flux_with_closest_peak_it]->GetMaximum();

    if (!InpCoeffFile.size()) {
      for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
        if (abs(flux_it - flux_with_closest_peak_it) ==
            1) {  // Neighbours are negative
          coeffs[flux_it] = -0.35;
        } else if (abs(flux_it - flux_with_closest_peak_it) == 0) {
          coeffs[flux_it] = 1;
        } else {  // Others start free
          coeffs[flux_it] = 0;
        }
      }
    }

    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      std::cout << "[INFO]: fluxpar" << flux_it
                << " uses flux: " << Fluxes[flux_it]->GetName() << std::endl;

      std::stringstream ss("");
      ss << "fluxpar" << flux_it;

      minimizer->SetParameter(flux_it, ss.str().c_str(), coeffs[flux_it], 0.1,
                              -CoeffLimit, CoeffLimit);
    }

    if (IsGauss) {
      TargetGauss->SetParameter(0, TargetPeakNorm);
      TargetGauss->SetParameter(1, GaussC);
      TargetGauss->SetParameter(2, GaussW);
      std::cout << "[INFO]: Peak: " << TargetPeakNorm
                << ", Central value: " << GaussC << ", width: " << GaussW
                << std::endl;
    }

    minimizer->SetMaxIterations(MaxMINUITCalls);
    minimizer->GetMinuit()->SetPrintLevel(0);

    minimizer->SetFCN(IsGauss ? TargetSumGauss : TargetSumChi2);

    double alist[2] = {double(MaxMINUITCalls), 1E-5 * TargetPeakNorm};

    fitstatus = minimizer->ExecuteCommand("MIGRAD", alist, 2);
    if (fitstatus) {
      std::cout << "[WARN]: Failed to find minimum (STATUS: " << fitstatus
                << ")." << std::endl;
    } else {
      minimizer->PrintResults(1, 0);
    }

    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      coeffs[flux_it] = minimizer->GetParameter(flux_it);
    }
  } // end fit

  SumHistograms(SummedFlux, coeffs, Fluxes);

  TFile *oupF =
      CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE");
  TDirectory *oupD = oupF;

  if (OutputDirectory.length()) {
    oupD = oupF->mkdir(OutputDirectory.c_str());
  }

  TDirectory *wD = oupD->mkdir("weighted_fluxes");
  wD->cd();

  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
    cl->Scale(coeffs[flux_it]);
    cl->SetDirectory(wD);
  }

  if (minimizer) {
    std::pair< std::vector<double>, std::vector<double> > XRangeBinCoeffs =
      SliceConfig::BuildXRangeBinsCoeffs(XRanges, coeffs);
    TH1D *coeffsH = new TH1D("Coeffs", "Coeffs;Off-axis position;Weight",
                             (XRangeBinCoeffs.first.size()-1),
                             XRangeBinCoeffs.first.data());

    for (size_t coeff_it = 0; coeff_it < (XRangeBinCoeffs.first.size()-1); coeff_it++) {
      coeffsH->SetBinContent(coeff_it + 1,
                             XRangeBinCoeffs.second[coeff_it]);
      coeffsH->SetBinError(coeff_it + 1, 0);
    }
  }

  oupD->cd();

  SummedFlux->SetDirectory(oupD);
  SummedFlux->SetName("BestFit");
  if (IsGauss) {
    TGraph *tGauss = new TGraph(1);
    tGauss->Set(1E4 - 1);

    double min = SummedFlux->GetXaxis()->GetBinLowEdge(1);
    double step = (SummedFlux->GetXaxis()->GetBinUpEdge(
                       SummedFlux->GetXaxis()->GetNbins()) -
                   SummedFlux->GetXaxis()->GetBinLowEdge(1)) /
                  double(1E4);
    for (size_t i = 1; i < 1E4; ++i) {
      double enu = min + i * step;
      double val = TargetGauss->Eval(enu);
      if (val != val) {
        continue;
      }
      tGauss->SetPoint(i - 1, enu, val);
    }

    tGauss->Write("Target_gauss");

  } else {
    TargetFlux->SetDirectory(oupD);
    TargetFlux->SetName("Target");

    TargetFlux_orignorm->SetDirectory(oupF);
    TargetFlux_orignorm->SetName("Target_input_normalisation");

    OscFlux->SetDirectory(oupD);
    OscFlux->SetName("InputFlux");

    OscFlux_orignorm->SetDirectory(oupD);
    OscFlux_orignorm->SetName("InputFlux_input_normalisation");
  }

  std::cout << "Used " << Fluxes.size() << " fluxes to fit "
            << (FitBinHigh - FitBinLow) << " bins. (Chi2 = " << Chi2_last
            << ", OOR = " << OOR_last << ", Reg = " << reg_last << ")."
            << std::endl;

  if (fitstatus > 0) {
    std::cout << "[WARN]: Failed to find minimum (STATUS: " << fitstatus << ")."
              << std::endl;
  }

  TTree *CoeffTree = new TTree("CoeffTree", "");

  double XRange[2];
  double Coeff;

  CoeffTree->Branch("XRange", &XRange, "XRange[2]/D");
  CoeffTree->Branch("Coeff", &Coeff, "Coeff/D");

  for (size_t i = 0; i < Fluxes.size(); ++i) {
    XRange[0] = XRanges[i].first * 100.0;
    XRange[1] = XRanges[i].second * 100.0;
    Coeff = coeffs[i];

    CoeffTree->Fill();
  }

  // Write out configuration trees.

  TTree *SliceConfigTree = new TTree("SliceConfigTree", "");
  SliceConfig *sc = SliceConfig::MakeTreeWriter(SliceConfigTree);
  for(size_t i = 0; i < Fluxes.size(); ++i){
    sc->XRange[0] = XRanges[i].first * 100.0; // Assumes tree written in off-axis position
    sc->XRange[1] = XRanges[i].second * 100.0; // Assumes tree written in off-axis position
    sc->Coeff = coeffs[i];

    SliceConfigTree->Fill();
  }

  TTree *FluxFitResultsTree = new TTree("FluxFitResultsTree", "");
  FluxFitResultsTreeReader *fr = FluxFitResultsTreeReader::MakeTreeWriter(FluxFitResultsTree, IsGauss);

  fr->NFluxes = Fluxes.size();
  fr->NIterations = NCalls;
  fr->Chi2 = Chi2_last;
  fr->RegularisationPenalty = reg_last;
  fr->FitRange[0] = FitBetween_low;
  fr->FitRange[1] = FitBetween_high;

  if (!IsGauss) {
    fr->OutOfRangePenalty = OOR_last;
    fr->NDOverFDFitScaleFactor = (SummedFlux->GetMaximum() / OscFlux_orignorm->GetMaximum());
  } else {
    fr->GaussCenter_GeV = GaussC;
    fr->GaussWidth_GeV = GaussW;
  }
  FluxFitResultsTree->Fill();

  if (!IsGauss) {
    OscillationParameters op_in("OscConfigTree", InputTargetFluxFile);
    TTree *OscConfigTree = new TTree("OscConfigTree", "");
    OscillationParameters *op_copy = OscillationParameters::MakeTreeWriter(OscConfigTree);
    op_copy->Copy(op_in);
    OscConfigTree->Fill();
  }

  oupF->Write();
  oupF->Close();
}
