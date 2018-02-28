#include "DetectorStop.hxx"
#include "GENIESplineReader.hxx"
#include "Utils.hxx"

#include "TF1.h"
#include "TFile.h"
#include "TFitter.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TKey.h"
#include "TList.h"
#include "TMinuit.h"

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

std::string FluxesFile, inpFluxHistsPattern;
std::string FluxHist2DName;
std::string inpFile, inpHistName;
std::string oupFile;
std::string oupDir;
std::vector<DetectorStop> detStops;
std::vector<std::pair<size_t, size_t> > UsedDetStopSlices;
std::string runPlanCfg;
std::string InpCoeffFile, InpCoeffDir;

std::string InpGConfigFile;
std::string GXMLFile;

bool UPDATEOutputFile = false;

std::vector<std::pair<std::string, std::string> > XSecComponentInputs;
std::map<std::string, TH1D *> XSecComponents;
TH1D *TotalXSec;

bool IsGauss = false;
double GaussC, GaussW;
double TargetPeakNorm;

double CoeffLimit = 30;
double RegFactor = 0xdeadbeef;

int MaxMINUITCalls = 50000;

double FitBetween_low = 0xdeadbeef, FitBetween_high = 0xdeadbeef;
int binLow, binHigh;
bool MultiplyChi2ContribByBinWidth = false;

double detYZ_m = 0xdeadbeef;
double detDensity_kgm3 = 0xdeadbeef;
double POT = 0xdeadbeef;

std::vector<std::pair<double, double> > IncludedOffAxisRange_2D_inputs;
std::vector<double> FluxOffaxisPositions_2D_inputs;
std::vector<bool> ApplyReg;
std::vector<double> SliceXWidth_m;

double referencePOT = std::numeric_limits<double>::min();
std::vector<double> MeasurementFactor;
std::vector<TH1D *> Fluxes;

TH1D *SummedFlux;

TH1D *TargetFlux;
TH1D *OscFlux;

TF1 *TargetGauss;

std::vector<std::pair<std::string, TGraph> > GENIEXSecs;

enum OutOfRangeModeEnum { kIgnore = 0, kZero, kExponentialDecay };
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

std::vector<double> InterpolatedOAAValues;

void BuildTargetFlux(TH1D *OscFlux) {
  TargetFlux = static_cast<TH1D *>(OscFlux->Clone());
  TargetFlux->SetDirectory(NULL);

  for (Int_t bi_it = 1; bi_it < binLow; ++bi_it) {
    if ((OutOfRangeMode != kIgnore) &&
        (OutOfRangeSide == kBoth || OutOfRangeSide == kLeft)) {
      double target = 0;
      if (OutOfRangeMode == kExponentialDecay) {
        double enu_first_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(binLow);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_first_counted_bin = TargetFlux->GetBinContent(binLow);
        double enu_bottom_bin = TargetFlux->GetXaxis()->GetBinCenter(1);
        double sigma5_range = enu_first_counted_bin - enu_bottom_bin;
        target =
            content_first_counted_bin *
            exp(-ExpDecayRate * (enu_first_counted_bin - enu) / sigma5_range);
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }

  for (Int_t bi_it = binLow; bi_it < (binHigh + 1); ++bi_it) {
    TargetFlux->SetBinContent(bi_it, OscFlux->GetBinContent(bi_it));
    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }

  for (Int_t bi_it = (binHigh + 1);
       bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
    if ((OutOfRangeMode != kIgnore) &&
        (OutOfRangeSide == kBoth || OutOfRangeSide == kRight)) {
      double target = 0;
      if (OutOfRangeMode == kExponentialDecay) {
        double enu_last_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(binHigh);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_last_counted_bin = TargetFlux->GetBinContent(binHigh);
        double enu_top_bin = TargetFlux->GetXaxis()->GetBinCenter(
            TargetFlux->GetXaxis()->GetNbins());
        double sigma5_range = enu_top_bin - enu_last_counted_bin;
        target =
            content_last_counted_bin *
            exp(-ExpDecayRate * (enu - enu_last_counted_bin) / sigma5_range);
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, 0.01 * TargetFlux->GetBinContent(bi_it));
  }
}

void TargetSumChi2(int &nDim, double *gout, double &result, double coeffs[],
                   int flg) {
  SumHistograms(SummedFlux, coeffs, Fluxes);

  double sumdiff = 0;
  double OOR = 0;

  if ((OutOfRangeMode != kIgnore)) {
    if ((OutOfRangeSide == kBoth || OutOfRangeSide == kLeft)) {
      for (Int_t bi_it = 1; bi_it < binLow; ++bi_it) {
        double diff;
        if (!UseNuPrismChi2) {
          double err =
              (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
               SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it));

          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);

        } else {
          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              (0.0001 * SummedFlux->GetBinContent(bi_it) *
                               SummedFlux->GetBinContent(bi_it)));
        }

        if (diff && !std::isnormal(diff)) {
          std::cout << "[INFO]: Found invalid diff, bin " << bi_it
                    << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it)
                    << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it) << std::endl;
          throw;  // exit(1);
        }

        sumdiff += diff;
        OOR += diff;
      }
    }

    if ((OutOfRangeSide == kBoth || OutOfRangeSide == kRight)) {
      for (Int_t bi_it = (binHigh + 1);
           bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
        double diff;
        if (!UseNuPrismChi2) {
          double err =
              (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
               SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it));

          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              err);

        } else {
          diff = OORFactor * ((TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) *
                              (TargetFlux->GetBinContent(bi_it) -
                               SummedFlux->GetBinContent(bi_it)) /
                              (0.0001 * SummedFlux->GetBinContent(bi_it) *
                               SummedFlux->GetBinContent(bi_it)));
        }

        if (diff && !std::isnormal(diff)) {
          std::cout << "[INFO]: Found invalid diff, bin " << bi_it
                    << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it)
                    << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                    << TargetFlux->GetBinError(bi_it) << std::endl;
          throw;  // exit(1);
        }

        sumdiff += diff;
        OOR += diff;
      }
    }
  }

  for (Int_t bi_it = binLow; bi_it < (binHigh + 1); ++bi_it) {
    if (!UseNuPrismChi2) {
      double err =
          (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
           SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it));

      sumdiff += ((TargetFlux->GetBinContent(bi_it) -
                   SummedFlux->GetBinContent(bi_it)) *
                  (TargetFlux->GetBinContent(bi_it) -
                   SummedFlux->GetBinContent(bi_it)) /
                  err);
    } else {
      sumdiff += ((TargetFlux->GetBinContent(bi_it) -
                   SummedFlux->GetBinContent(bi_it)) *
                  (TargetFlux->GetBinContent(bi_it) -
                   SummedFlux->GetBinContent(bi_it)) /
                  (0.0001 * SummedFlux->GetBinContent(bi_it) *
                   SummedFlux->GetBinContent(bi_it)));
    }

    if (sumdiff && !std::isnormal(sumdiff)) {
      std::cout << "[INFO]: Found invalid diff, bin " << bi_it
                << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it)
                << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it) << std::endl;
      throw;  // exit(1);
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

  std::cout << "[INFO]: Target flux chi2: " << result << " (reg = " << reg
            << ", OOR = " << OOR << " )." << std::endl;
}

void TargetSumGauss(int &nDim, double *gout, double &result, double coeffs[],
                    int flg) {
  SumHistograms(SummedFlux, coeffs, Fluxes);

  double sumdiff = 0;
  for (Int_t bi_it = binLow; bi_it < binHigh + 1; ++bi_it) {
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
      std::cout << "[INFO]: Found invalid diff, bin " << bi_it
                << ", gauss: " << GaussEval << ":" << TargetPeakNorm
                << ", sum: " << SummedFlux->GetBinContent(bi_it) << ":"
                << SummedFlux->GetBinError(bi_it) << std::endl;
      throw;  // exit(1);
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
         "\t-f <ROOT file[,hist name pattern]> : The input flux "
         "histograms used in  \n"
         "\t                                     the fit. If a pattern "
         "is not       \n"
         "\t                                     passed, then all TH1Ds "
         "in the input\n"
         "\t                                     are used.               "
         "           \n"
         "\n"
         "\t-h <HistName>                      : Input 2D flux histogram,"
         " Y bins \n"
         "\t                                     correspond to different "
         "fluxes.\n"
         "\n"
         "\t-r <RunPlan.XML>                   : An XML file specifying "
         "a run plan  \n"
         "\t                                     to make event rate "
         "predictions and \n"
         "\t                                     measured spectra for. "
         "See          \n"
         "\t                                     documentation for XML "
         "structure.   \n"
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
         "  Target options:                                              "
         "           \n"
         "\t-t <ROOT file,hist name>           : The histogram of the "
         "target flux to\n"
         "\t                                     fit.                    "
         "           \n"
         "\n"
         "\t-g <mean,width>                    : Use a gaussian target "
         "distribution \n"
         "\t                                     instead of a target "
         "flux shape.    \n"
         "\n"
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
         "\t-x <ROOT file, hist name>          : Add xsec component for "
         "making event\n"
         "\t                                     rate predictions.       "
         "           \n"
         "\n"
         "\t-rg <regularisation factor>        : Adds neighbouring       "
         "coefficient\n"
         "\t                                     regularisation.         "
         "           \n"
         "\n"
         "\t-l <min val>,<max val>              : Fit between min and max."
         "Outside \n"
         "\t                                     of this range, -m       "
         "determines\n"
         "\t                                     behavior.             \n"
         "\n"
         "\t-m <0,1,2>                         : Out of range behavior.  "
         "          \n"
         "\t                                     0: Ignore out of range "
         "bins.      \n"
         "\t                                     1: Force out of range "
         "bins to 0.  \n"
         "\t                                     2: exponential decay "
         "outside fit  \n"
         "\t                                        region. decay rate "
         "is          \n"
         "\t                                        determined by -ed."
         "\n"
         "\n"
         "\t-ed <decay rate>                   : For -m 2, controls "
         "decay rate.    \n"
         "\t                                     Default = 3, larger is "
         "faster     \n"
         "\t                                     decay.\n"
         "\n"
         "\t-of <out of range factor>          : Allow out of range to "
         "contribute less to the chi2\n "
         "\t                                     by this factor.\n"
         "\n"
         "\t-ms <out of range side>            : 0 = Include both low and high "
         "out of range E, \n"
         "\t                                     1 = include low E, 2 = "
         "include "
         "high E.\n"
         "\n"
         "\t-O <OAP_min,OAP_max[,OAP_min2,OAP_max2]> : Specify regions of off "
         "axis angle to include \n"
         "\t                                           in fit. Will only be "
         "applied to a `-h` style input.\n"
         "\n"
         "\t-gc <geniexsecconfig.xml>          : A xml config file containing "
         "xsec categories built from GENIE\n"
         "\t                                     spline names. An example "
         "should be found in \n"
         "\t                                     "
         "${DUNEPRISMTOOLSROOT}/configs/gxsec.conf.xml\n"
         "\n"
         "\t-gx <gxml.xml>                     : The GENIE XML file to read the"
         " splines from.\n"
         "\n"
         "\t-dYZ                               : YZ cross-sectional area of "
         "the detector (for event rate predictions).\n"
         "\n"
         "\t-dD                                : Detector target density (in "
         "kg_m3).\n"
         "\n"
         "\t-P                                 : POT exposure (for event rate "
         "predictions).\n"
         "\n"
         "\t-C                                 : Use NuPrism tools Chi2.\n\n"
         "\t-MY  <nbins to merge>              : Merge off-axis angle bins "
         "before splitting into\n\n"
         "\t                                     fluxes. Works with -h inputs.\n"
         "\n"
         "\t-MX  <nbins to merge>              : Merge neutrino energy bins "
         "before splitting into\n"
         "\t                                     fluxes. Works with -h inputs."
         "\n\n"
         "\t-I  <OA1>,<OA2>_<OAN>:<oa_step>,<OAN+1>,...\n "
         "\t                                   : Interpolate off axis flux "
         "positions from a -h \n"
         "\t                                     input. e.g. 1,2_6:2,7 would "
         "fit fluxes at \n"
         "\t                                     1, 2, 4, 6, and 7 m off axis."
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
        throw;  // exit(1);
      }
      GaussC = params[0];
      GaussW = params[1];
      if ((GaussC <= 0) || (GaussW <= 0)) {
        std::cout << "[ERROR]: Recieved " << argv[opt]
                  << " argument for -g, expected \"<GaussMean>,<GaussWidth>\", "
                     "where both values are > 0."
                  << std::endl;
        throw;  // exit(1);
      }
      IsGauss = true;
      TargetGauss = new TF1("tFunc", "gaus", 0, 20);
    } else if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for l, expected 2." << std::endl;
        throw;  // exit(1);
      }
      FitBetween_low = params[0];
      FitBetween_high = params[1];
    } else if (std::string(argv[opt]) == "-t") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -t, expected 2." << std::endl;
        throw;  // exit(1);
      }
      inpFile = params[0];
      inpHistName = params[1];
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
      oupFile = argv[++opt];
      UPDATEOutputFile = false;
    } else if (std::string(argv[opt]) == "-a") {
      oupFile = argv[++opt];
      UPDATEOutputFile = true;
    } else if (std::string(argv[opt]) == "-d") {
      oupDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-f") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() < 2) {
        FluxesFile = params[0];
        inpFluxHistsPattern = "*";
      } else {
        FluxesFile = params[0];
        inpFluxHistsPattern = params[1];
      }
    } else if (std::string(argv[opt]) == "-n") {
      MaxMINUITCalls = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-c") {
      CoeffLimit = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-m") {
      OutOfRangeMode = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-r") {
      runPlanCfg = argv[++opt];
    } else if (std::string(argv[opt]) == "-x") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() < 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -x, expected at least 2." << std::endl;
        throw;  // exit(1);
      }
      for (size_t xs_it = 1; xs_it < params.size(); ++xs_it) {
        XSecComponentInputs.push_back(std::make_pair(params[0], params[xs_it]));
      }
    } else if (std::string(argv[opt]) == "-h") {
      FluxHist2DName = argv[++opt];
    } else if (std::string(argv[opt]) == "-rg") {
      RegFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ed") {
      ExpDecayRate = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-B") {
      MultiplyChi2ContribByBinWidth = true;
    } else if (std::string(argv[opt]) == "-O") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");

      if (params.size() & 1) {
        std::cout << "[ERROR]: -O recieved: \"" << argv[opt]
                  << "\", which contains an odd number of inputs, expected an "
                     "even number."
                  << std::endl;
        throw;  // exit(1);
      }

      for (size_t i = 0; i < params.size(); i += 2) {
        IncludedOffAxisRange_2D_inputs.push_back(
            std::make_pair(params[i], params[i + 1]));
      }
    } else if (std::string(argv[opt]) == "-of") {
      OORFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ms") {
      OutOfRangeSide = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-gc") {
      InpGConfigFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-gx") {
      GXMLFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-dYZ") {
      detYZ_m = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-dD") {
      detDensity_kgm3 = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-P") {
      POT = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-C") {
      UseNuPrismChi2 = true;
    } else if (std::string(argv[opt]) == "-MY") {
      MergeOAABins = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-MX") {
      MergeENuBins = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-I") {
      std::vector<std::string> interpDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");
      InterpolatedOAAValues.clear();
      for (size_t vbd_it = 0; vbd_it < interpDescriptors.size(); ++vbd_it) {
        AppendVect(InterpolatedOAAValues,
                   BuildDoubleList(interpDescriptors[vbd_it]));
      }

      for (size_t bin_it = 1; bin_it < InterpolatedOAAValues.size(); ++bin_it) {
        if (InterpolatedOAAValues[bin_it] ==
            InterpolatedOAAValues[bin_it - 1]) {
          std::cout << "[INFO]: Removing duplciate interpolated oaa values ["
                    << bin_it << "] = " << InterpolatedOAAValues[bin_it]
                    << std::endl;
          InterpolatedOAAValues.erase(InterpolatedOAAValues.begin() + bin_it);
        }
      }

      for (size_t bin_it = 1; bin_it < InterpolatedOAAValues.size(); ++bin_it) {
        if (InterpolatedOAAValues[bin_it] < InterpolatedOAAValues[bin_it - 1]) {
          std::cout << "[ERROR]: Interpolated oaa value #" << bin_it << " = "
                    << InterpolatedOAAValues[bin_it] << ". however, #"
                    << (bin_it - 1) << " = "
                    << InterpolatedOAAValues[bin_it - 1] << std::endl;
          exit(1);
        }
      }
    } else if ((std::string(argv[opt]) == "-?") ||
               std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      throw;  // exit(1);
    }
    opt++;
  }
}
int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!IsGauss && (!inpFile.length() || !inpHistName.length())) {
    std::cout << "[ERROR]: No input file or input histogram name specified."
              << std::endl;
    throw;  // exit(1);
  }

  if (!oupFile.length()) {
    std::cout << "[ERROR]: No output file specified." << std::endl;
    throw;  // exit(1);
  }

  if (!IncludedOffAxisRange_2D_inputs.size() && inpHistName.length()) {
    IncludedOffAxisRange_2D_inputs.push_back(
        std::make_pair(-std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max()));
  }

  if (!IsGauss) {
    OscFlux = GetHistogram<TH1D>(inpFile, inpHistName);

    if (FitBetween_low == 0xdeadbeef) {
      binLow = 1;
      binHigh = OscFlux->GetXaxis()->GetNbins();
    } else {
      binLow = OscFlux->GetXaxis()->FindFixBin(FitBetween_low);
      binHigh = OscFlux->GetXaxis()->FindFixBin(FitBetween_high);
    }

    std::cout << "[INFO]: Fitting between bins: " << binLow << "--" << binHigh
              << std::endl;

    BuildTargetFlux(OscFlux);
  }

  if (FluxesFile.length()) {
    if (runPlanCfg.length()) {
      detStops = ReadDetectorStopConfig(runPlanCfg);
      referencePOT = detStops[0].POTExposure;

      TFile *ifl = CheckOpenFile(FluxesFile.c_str());

      for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
        detStops[ds_it].Read(ifl);

        std::vector<double> doubleupevrate;
        for (size_t ms_it = 0; ms_it < detStops[ds_it].GetNMeasurementSlices();
             ++ms_it) {
          if (detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it) < 0) {
            doubleupevrate.push_back(
                fabs(detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it)));
            std::cout << "[INFO]: Excluding Stop " << ds_it << ", slice "
                      << ms_it << " ("
                      << detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it)
                      << " m) from fit as it should have a mirrored slice."
                      << std::endl;
            continue;
          }

          if (doubleupevrate.size()) {
            bool found = false;
            for (size_t du_it = 0; du_it < doubleupevrate.size(); ++du_it) {
              if (fabs(doubleupevrate[du_it] -
                       detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it)) <
                  1E-3) {
                std::cout << "[INFO]: Stop " << ds_it << ", slice " << ms_it
                          << ", getting double stats from mirrored measurement "
                             "slice ("
                          << detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it)
                          << " m)." << std::endl;

                doubleupevrate.erase(doubleupevrate.begin() + du_it);
                found = true;
                break;
              }
            }
            if (found) {
              MeasurementFactor.push_back(2.0 * detStops[ds_it].POTExposure /
                                          referencePOT);
              std::cout << "MeasurementFactor " << MeasurementFactor.size()
                        << " = " << MeasurementFactor.back() << std::endl;
            } else {
              MeasurementFactor.push_back(detStops[ds_it].POTExposure /
                                          referencePOT);
              std::cout << "MeasurementFactor " << MeasurementFactor.size()
                        << " = " << MeasurementFactor.back() << std::endl;
            }
          } else {
            MeasurementFactor.push_back(detStops[ds_it].POTExposure /
                                        referencePOT);
            std::cout << "MeasurementFactor " << MeasurementFactor.size()
                      << " = " << MeasurementFactor.back() << std::endl;
          }

          Fluxes.push_back(detStops[ds_it].GetFluxForSpecies(ms_it, 14));
          UsedDetStopSlices.push_back(std::make_pair(ds_it, ms_it));
        }

        if (MeasurementFactor.size() != Fluxes.size()) {
          std::cout << "[ERROR]: Didn't determine coeff POT scale for every "
                       "flux slice."
                    << std::endl;
          throw;
        }

        if (doubleupevrate.size()) {
          std::cout << "[ERROR]: Still had " << doubleupevrate.size()
                    << " ignored negative slices that did not have positive "
                       "offset mirror slices."
                    << std::endl;
          for (size_t du_it = 0; du_it < doubleupevrate.size(); ++du_it) {
            std::cout << "\t\t" << doubleupevrate[du_it];
          }
          throw;
        }
      }
      if (!Fluxes.size()) {
        std::cout << "[ERROR]: Found no input fluxes from input run plan: \""
                  << runPlanCfg << "\"." << std::endl;
        throw;  // exit(1);
      } else {
        std::cout << "[INFO]: Found " << Fluxes.size() << " input fluxes."
                  << std::endl;
      }

      ifl->Close();
      delete ifl;
    } else if (FluxHist2DName.size()) {
      TH2D *Flux2D = GetHistogram<TH2D>(FluxesFile, FluxHist2DName);
      if (!Flux2D) {
        std::cout << "[ERROR]: Found no input flux with name: \""
                  << FluxHist2DName << "\" in file: \"" << FluxesFile << "\"."
                  << std::endl;
        throw;  // exit(1);
      }

      if (MergeOAABins) {
        Flux2D->RebinY(MergeOAABins);
        Flux2D->Scale(1.0 / double(MergeOAABins));
      }
      if (MergeENuBins) {
        Flux2D->RebinX(MergeENuBins);
        Flux2D->Scale(1.0 / double(MergeENuBins));
      }

      std::vector<std::pair<double, TH1D *> > Fluxes_and_OAPs;
      if (InterpolatedOAAValues.size()) {
        Fluxes_and_OAPs =
            InterpolateSplitTH2D(Flux2D, true, InterpolatedOAAValues);

        for (size_t f_it = 0; f_it < Fluxes_and_OAPs.size(); ++f_it) {
          ApplyReg.push_back(true);
          FluxOffaxisPositions_2D_inputs.push_back(Fluxes_and_OAPs[f_it].first);
          Fluxes.push_back(Fluxes_and_OAPs[f_it].second);
        }
      } else {
        Fluxes_and_OAPs = SplitTH2D(Flux2D, true, 0);

        if (!IncludedOffAxisRange_2D_inputs.size()) {
          for (size_t f_it = 0; f_it < Fluxes_and_OAPs.size(); ++f_it) {
            ApplyReg.push_back(true);
            FluxOffaxisPositions_2D_inputs.push_back(
                Fluxes_and_OAPs[f_it].first);
            Fluxes.push_back(Fluxes_and_OAPs[f_it].second);
          }
        } else {
          size_t this_kept_range = -1,
                 last_kept_range = std::numeric_limits<size_t>::max();
          for (size_t f_it = 0; f_it < Fluxes_and_OAPs.size(); ++f_it) {
            bool keep = false;
            for (size_t inc_it = 0;
                 inc_it < IncludedOffAxisRange_2D_inputs.size(); ++inc_it) {
              bool ge_low =
                  ((Fluxes_and_OAPs[f_it].first >
                    IncludedOffAxisRange_2D_inputs[inc_it].first) ||
                   (fabs(Fluxes_and_OAPs[f_it].first -
                         IncludedOffAxisRange_2D_inputs[inc_it].first)) < 1E-5);
              bool le_up =
                  ((Fluxes_and_OAPs[f_it].first <
                    IncludedOffAxisRange_2D_inputs[inc_it].second) ||
                   (fabs(Fluxes_and_OAPs[f_it].first -
                         IncludedOffAxisRange_2D_inputs[inc_it].second)) <
                       1E-5);
              if (ge_low && le_up) {
                keep = true;
                this_kept_range = inc_it;
                break;
              }
            }
            if (!keep) {
              delete Fluxes_and_OAPs[f_it].second;
              continue;
            }

            std::cout << "[INFO]: Keeping flux slice at "
                      << Fluxes_and_OAPs[f_it].first << " m OAP." << std::endl;

            if ((last_kept_range != std::numeric_limits<size_t>::max()) &&
                (this_kept_range != last_kept_range)) {
              ApplyReg.back() = false;
              std::cout << "[INFO]: Ignoring reg factor for flux at "
                        << FluxOffaxisPositions_2D_inputs.back() << " m OAP."
                        << std::endl;
            }

            ApplyReg.push_back(true);
            FluxOffaxisPositions_2D_inputs.push_back(
                Fluxes_and_OAPs[f_it].first);
            Fluxes.push_back(Fluxes_and_OAPs[f_it].second);

            Int_t ybi_it =
                Flux2D->GetYaxis()->FindFixBin(Fluxes_and_OAPs[f_it].first);
            SliceXWidth_m.push_back(Flux2D->GetYaxis()->GetBinWidth(ybi_it));
            last_kept_range = this_kept_range;
          }
        }
      }
      if (!Fluxes_and_OAPs.size()) {
        std::cout << "[ERROR]: Couldn't find any fluxes in split TH2D."
                  << std::endl;
        throw;  // exit(1);
      }
      std::cout << "[INFO]: Found " << Fluxes_and_OAPs.size()
                << " input fluxes." << std::endl;

    } else if (inpFluxHistsPattern.size()) {
      Fluxes = GetHistograms<TH1D>(FluxesFile, inpFluxHistsPattern);
      if (!Fluxes.size()) {
        std::cout << "[ERROR]: Found no input fluxes matching pattern: \""
                  << inpFluxHistsPattern << "\" in file: \"" << FluxesFile
                  << "\"." << std::endl;
        throw;  // exit(1);
      } else {
        std::cout << "[INFO]: Found " << Fluxes.size() << " input fluxes."
                  << std::endl;
      }
    }

  } else {
    std::cout << "[ERROR]: Expected either -f (h) or -r options to be passed."
              << std::endl;
    throw;  // exit(1);
  }

  if (!Fluxes.size()) {
    std::cout << "[ERROR]: Found no input fluxes." << std::endl;
    throw;  // exit(1);
  }

  SummedFlux = static_cast<TH1D *>(Fluxes[0]->Clone());
  SummedFlux->SetDirectory(NULL);
  SummedFlux->Reset();

  if (IsGauss) {
    if (FitBetween_low == 0xdeadbeef) {
      binLow = 1;
      binHigh = SummedFlux->GetXaxis()->GetNbins();
    } else {
      binLow = SummedFlux->GetXaxis()->FindFixBin(FitBetween_low);
      binHigh = SummedFlux->GetXaxis()->FindFixBin(FitBetween_high);
    }
  }

  double *coeffs = new double[Fluxes.size()];
  if (InpCoeffFile.size()) {
    TFile *icf = new TFile(InpCoeffFile.c_str(), "READ");
    if (!icf || !icf->IsOpen()) {
      std::cout << "[ERROR]: Couldn't open coeff input file: " << InpCoeffFile
                << std::endl;
      throw;  // exit(1);
    }

    if (InpCoeffDir.size() && (InpCoeffDir.back() != '/')) {
      InpCoeffDir += "/";
    }

    TGraph *cg =
        dynamic_cast<TGraph *>(icf->Get((InpCoeffDir + "coeffs").c_str()));
    if (!cg) {
      std::cout << "[ERROR]: Couldn't get TGraph \'coeffs\' from output file: "
                << InpCoeffFile << std::endl;
      throw;  // exit(1);
    }

    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      double x, y;
      cg->GetPoint(flux_it, x, y);
      coeffs[flux_it] = y;
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
    TargetFlux->Scale(OnAxisPeak / TargetMax);
    OscFlux->Scale(OnAxisPeak / TargetMax);
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
  }

  SumHistograms(SummedFlux, coeffs, Fluxes);

  TFile *oupF =
      new TFile(oupFile.c_str(), UPDATEOutputFile ? "UPDATE" : "RECREATE");
  if (!oupF || !oupF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open output file: " << oupFile << std::endl;
    throw;  // exit(1);
  }

  TDirectory *oupD = oupF;

  if (oupDir.length()) {
    oupD = oupF->mkdir(oupDir.c_str());
  }

  TDirectory *wD = oupD->mkdir("weighted_fluxes");
  wD->cd();

  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
    cl->Scale(coeffs[flux_it]);
    cl->SetDirectory(wD);
  }

  if (minimizer) {
    TH1D *coeffsH = new TH1D("Coeffs", "Coeffs;Off-axis bin;Weight",
                             Fluxes.size(), 0, Fluxes.size());
    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      coeffsH->SetBinContent(flux_it + 1,
                             coeffs[(Fluxes.size() - 1) - flux_it]);
      coeffsH->SetBinError(
          flux_it + 1, minimizer->GetParError((Fluxes.size() - 1) - flux_it));
    }
  }

  oupD->cd();

  std::stringstream ss("");
  if (detStops.size()) {
    TGraph peak_coeffs(1);
    peak_coeffs.Set(Fluxes.size());
    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      peak_coeffs.SetPoint(
          flux_it,
          detStops[UsedDetStopSlices[flux_it].first].GetAbsoluteOffsetOfSlice(
              UsedDetStopSlices[flux_it].second),
          coeffs[flux_it]);
    }
    peak_coeffs.Write("coeffs");
  } else {
    TGraph peak_coeffs(1);
    peak_coeffs.Set(Fluxes.size());
    for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
      peak_coeffs.SetPoint(flux_it,
                           FluxOffaxisPositions_2D_inputs.size()
                               ? FluxOffaxisPositions_2D_inputs[flux_it]
                               : flux_it,
                           coeffs[flux_it]);
    }
    peak_coeffs.Write("coeffs");
  }

  oupD->cd();

  wD = oupD->mkdir("flux_build_unordered");
  wD->cd();

  TH1D *cl_0 = static_cast<TH1D *>(Fluxes[0]->Clone());
  cl_0->Scale(coeffs[0]);
  cl_0->SetDirectory(wD);
  ss.str("flux_0");
  cl_0->SetName(ss.str().c_str());
  for (size_t flux_it = 1; flux_it < Fluxes.size(); flux_it++) {
    TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
    cl->Scale(coeffs[flux_it]);
    cl->SetDirectory(nullptr);

    cl_0->Add(cl);
    ss.str("");
    ss << "flux_sum_0-" << flux_it;

    cl_0->Write(ss.str().c_str());

    delete cl;
  }

  oupD->cd();

  wD = oupD->mkdir("flux_build_ordered");
  wD->cd();

  std::vector<double> order_coeff;
  std::vector<size_t> order_idxs;

  order_coeff.resize(Fluxes.size());
  order_idxs.resize(Fluxes.size());

  for (size_t fl_it = 0; fl_it < Fluxes.size(); ++fl_it) {
    order_coeff[fl_it] = coeffs[fl_it];
    order_idxs[fl_it] = fl_it;
  }

  // Bubbbbles
  bool IsOrdered = false;
  while (!IsOrdered) {
    IsOrdered = true;
    for (size_t fl_it = 0; fl_it < (Fluxes.size() - 1); ++fl_it) {
      if (fabs(order_coeff[fl_it]) < fabs(order_coeff[fl_it + 1])) {
        std::swap(order_coeff[fl_it], order_coeff[fl_it + 1]);
        std::swap(order_idxs[fl_it], order_idxs[fl_it + 1]);
        IsOrdered = false;
      }
    }
  }

  std::cout << "[INFO]: Ordered by coeff: " << std::endl;
  for (size_t fl_it = 0; fl_it < (Fluxes.size() - 1); ++fl_it) {
    std::cout << "\t" << fl_it << ": idx: " << order_idxs[fl_it]
              << ", coeff: " << order_coeff[fl_it] << std::endl;
  }

  TH1D *clordered_0 = static_cast<TH1D *>(Fluxes[order_idxs[0]]->Clone());
  clordered_0->Scale(order_coeff[0]);
  clordered_0->SetDirectory(wD);
  ss.str("flux_0");
  clordered_0->Write(ss.str().c_str());
  for (size_t flux_it = 1; flux_it < Fluxes.size(); flux_it++) {
    TH1D *clordered = static_cast<TH1D *>(Fluxes[order_idxs[flux_it]]->Clone());
    clordered->Scale(order_coeff[flux_it]);
    clordered->SetDirectory(nullptr);

    clordered_0->Add(clordered);
    ss.str("");
    ss << "flux_order_0-" << flux_it;

    clordered_0->Write(ss.str().c_str());

    delete clordered;
  }

  oupD->cd();

  if (XSecComponentInputs.size() && detStops.size()) {
    for (std::pair<std::string, std::string> hdescript : XSecComponentInputs) {
      TH1D *hist = GetHistogram<TH1D>(hdescript.first, hdescript.second);
      std::cout << "[INFO]: Got XSec component: " << hist->GetName()
                << std::endl;
      XSecComponents[hist->GetName()] = static_cast<TH1D *>(hist->Clone());
      XSecComponents[hist->GetName()]->SetDirectory(nullptr);
    }

    std::vector<TH1D *> EvrRates;

    for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
      detStops[ds_it].PredictEventRates(XSecComponents, 14);
      detStops[ds_it].Write();
    }

    for (std::pair<size_t, size_t> usedSlices : UsedDetStopSlices) {
      EvrRates.push_back(detStops[usedSlices.first].GetTotalPredictedEventRate(
          14, usedSlices.second));
    }

    TH1D *evr = static_cast<TH1D *>(SummedFlux->Clone());

    std::vector<double> evrcoeffs = MeasurementFactor;
    for (size_t co_it = 0; co_it < evrcoeffs.size(); ++co_it) {
      evrcoeffs[co_it] = coeffs[co_it] / evrcoeffs[co_it];
    }

    oupD->cd();

    wD = oupD->mkdir("evr_build_unordered");
    wD->cd();

    TH1D *ev_0 = static_cast<TH1D *>(EvrRates[0]->Clone());
    ev_0->Scale(evrcoeffs[0]);
    ev_0->SetDirectory(wD);
    std::stringstream ss("flux_0");
    ev_0->Write(ss.str().c_str(), TObject::kOverwrite);
    for (size_t evr_it = 1; evr_it < EvrRates.size(); evr_it++) {
      TH1D *ev = static_cast<TH1D *>(EvrRates[evr_it]->Clone());
      ev->Scale(evrcoeffs[evr_it]);
      ev->SetDirectory(nullptr);

      ev_0->Add(ev);
      ss.str("");
      ss << "evr_sum_0-" << evr_it;

      ev_0->Write(ss.str().c_str(), TObject::kOverwrite);

      delete ev;
    }

    oupD->cd();

    SumHistograms(evr, evrcoeffs.data(), EvrRates);
    evr->SetDirectory(nullptr);
    evr->GetYaxis()->SetTitle("Events / GeV");
    evr->Write("PredictedMeasurement");
  }

  if (InpGConfigFile.size() && GXMLFile.size()) {
    if ((detYZ_m == 0xdeadbeef) || (detDensity_kgm3 == 0xdeadbeef) ||
        (POT == 0xdeadbeef)) {
      std::cout << "[WARN]: Cannot produce evrate predictions as detector and "
                   "POT are not set."
                << std::endl;
    } else {
      GENIEXSecReader gxr(InpGConfigFile);

      gxr.Read(GXMLFile);

      GENIEXSecs = gxr.GetTGraphs(Fluxes[0]->GetXaxis()->GetBinCenter(1),
                                  Fluxes[0]->GetXaxis()->GetBinCenter(
                                      Fluxes[0]->GetXaxis()->GetNbins()));

      oupD->cd();

      wD = oupD->mkdir("GENIE_evr");
      wD->cd();

      for (size_t xsec_it = 0; xsec_it < GENIEXSecs.size(); ++xsec_it) {
        GENIEXSecs[xsec_it].second.Write(GENIEXSecs[xsec_it].first.c_str(),
                                         TObject::kOverwrite);
      }

      std::vector<TH1D *> CCIncEvr;

      static double const mass_proton_kg = 1.6727E-27;   // Proton mass in kg
      static double const mass_neutron_kg = 1.6750E-27;  // Neutron mass in kg
      static double const mass_nucleon_kg =
          (mass_proton_kg + mass_neutron_kg) / 2.;

      std::stringstream ss("");
      for (size_t fl_it = 0; fl_it < Fluxes.size(); ++fl_it) {
        CCIncEvr.push_back(static_cast<TH1D *>(Fluxes[fl_it]->Clone()));
        ss.str("");
        ss << Fluxes[fl_it]->GetName() << "_GCCInc";
        CCIncEvr.back()->SetNameTitle(ss.str().c_str(),
                                      "E_{#nu} (GeV);CCInc Event Rate;");
        CCIncEvr.back()->SetDirectory(wD);

        for (Int_t bi_it = 1;
             bi_it < CCIncEvr.back()->GetXaxis()->GetNbins() + 1; ++bi_it) {
          double xsec_cm2_gev = 0;
          double Enu = CCIncEvr.back()->GetXaxis()->GetBinCenter(bi_it);

          for (size_t xsec_it = 0; xsec_it < GENIEXSecs.size(); ++xsec_it) {
            xsec_cm2_gev += GENIEXSecs[xsec_it].second.Eval(Enu);
          }

          double flux_pcm2_pgev_pnuc_pPOT =
              CCIncEvr.back()->GetBinContent(bi_it);

          double nnuc = detYZ_m * SliceXWidth_m[fl_it] * detDensity_kgm3 /
                        mass_nucleon_kg;

          double NEv = flux_pcm2_pgev_pnuc_pPOT * xsec_cm2_gev * nnuc * POT;

          if (fl_it == 0 && bi_it == 80) {
            std::cout << "ENu: " << Enu << ", TXSec: " << xsec_cm2_gev
                      << ", DetMass: "
                      << detYZ_m * SliceXWidth_m[fl_it] * detDensity_kgm3
                      << ", POT: " << POT << ", NEv: " << NEv << std::endl;
          }

          CCIncEvr.back()->SetBinContent(bi_it, NEv);
          double frac_flux_error = Fluxes[fl_it]->GetBinError(bi_it) /
                                   Fluxes[fl_it]->GetBinContent(bi_it);
          double frac_evr_error = 1.0 / sqrt(NEv);

          if ((NEv < 1E-6)) {
            frac_evr_error = 0;
          }

          if (Fluxes[fl_it]->GetBinContent(bi_it) < 1E-15) {
            frac_flux_error = 0;
          }

          CCIncEvr.back()->SetBinError(
              bi_it, NEv * sqrt(frac_evr_error * frac_evr_error +
                                frac_flux_error * frac_flux_error));
        }
      }

      TH1D *evr = static_cast<TH1D *>(SummedFlux->Clone());
      SumHistograms(evr, coeffs, CCIncEvr);
      evr->SetDirectory(wD);
      evr->GetYaxis()->SetTitle("Events / GeV");
      evr->SetName("PredictedMeasurement");

      TH1D *target = static_cast<TH1D *>(SummedFlux->Clone());
      target->SetDirectory(wD);
      target->GetYaxis()->SetTitle("Events / GeV");
      target->SetName("PredictedTargetMeasurement");

      target->Reset();

      for (Int_t bi_it = 1; bi_it < target->GetXaxis()->GetNbins() + 1;
           ++bi_it) {
        if (IsGauss) {
          target->SetBinContent(
              bi_it,
              TargetGauss->Eval(target->GetXaxis()->GetBinCenter(bi_it)));
        } else {
          target->SetBinContent(bi_it, OscFlux->GetBinContent(bi_it));
        }

        double xsec_cm2_gev = 0;
        double Enu = target->GetXaxis()->GetBinCenter(bi_it);

        for (size_t xsec_it = 0; xsec_it < GENIEXSecs.size(); ++xsec_it) {
          xsec_cm2_gev += GENIEXSecs[xsec_it].second.Eval(Enu);
        }

        double flux_pcm2_pgev_pnuc_pPOT = target->GetBinContent(bi_it);

        double nnuc =
            detYZ_m * SliceXWidth_m[0] * detDensity_kgm3 / mass_nucleon_kg;

        double NEv = flux_pcm2_pgev_pnuc_pPOT * xsec_cm2_gev * nnuc * POT;
        target->SetBinContent(bi_it, NEv);
        double frac_flux_error =
            target->GetBinError(bi_it) / target->GetBinContent(bi_it);
        double frac_evr_error = 1.0 / sqrt(NEv);

        if ((NEv < 1E-6)) {
          frac_evr_error = 0;
        }

        if (target->GetBinContent(bi_it) < 1E-15) {
          frac_flux_error = 0;
        }

        target->SetBinError(bi_it,
                            NEv * sqrt(frac_evr_error * frac_evr_error +
                                       frac_flux_error * frac_flux_error));
      }

      oupD->cd();
    }
  }

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

    if (XSecComponents.size() && detStops.size()) {
      TGraph *tGauss_evt = new TGraph(1);
      tGauss_evt->Set(1E4 - 1);

      TH1D *totalXSec =
          static_cast<TH1D *>(XSecComponents.begin()->second->Clone());
      totalXSec->Reset();
      totalXSec->SetDirectory(nullptr);
      for (std::map<std::string, TH1D *>::iterator xs_it =
               XSecComponents.begin();
           xs_it != XSecComponents.end(); ++xs_it) {
        totalXSec->Add(xs_it->second);
      }

      const double mass_proton_kg = 1.6727E-27;   // Proton mass in kg
      const double mass_neutron_kg = 1.6750E-27;  // Neutron mass in kg
      const double mass_nucleon_kg = (mass_proton_kg + mass_neutron_kg) / 2.;
      double Vol = detStops[0].MeasurementRegionWidth *
                   detStops[0].DetectorFiducialDepth *
                   detStops[0].DetectorFiducialHeight;
      double Mass = Vol * detStops[0].FiducialVolumeDensity;
      double NNucleons = Mass / mass_nucleon_kg;

      for (size_t i = 1; i < 1E4; ++i) {
        double enu = min + i * step;
        double val = TargetGauss->Eval(enu);
        double b_evr =
            val * totalXSec->Interpolate(enu) * referencePOT * NNucleons;
        tGauss_evt->SetPoint(i - 1, enu, b_evr);
      }
      tGauss_evt->Write("Target_gauss_EvrRate");
    }

  } else {
    TargetFlux->SetDirectory(oupD);
    TargetFlux->SetName("Target");

    OscFlux->SetDirectory(oupD);
    OscFlux->SetName("InputFlux");
  }

  int dum1, dum4 = 0;
  double dum2, dum3;
  IsGauss ? TargetSumGauss(dum1, &dum2, dum3, coeffs, dum4)
          : TargetSumChi2(dum1, &dum2, dum3, coeffs, dum4);
  std::cout << "Used " << Fluxes.size() << " fluxes to fit "
            << (binHigh - binLow) << " bins." << std::endl;

  if (fitstatus > 0) {
    std::cout << "[WARN]: Failed to find minimum (STATUS: " << fitstatus << ")."
              << std::endl;
  }

  oupF->Write();
  oupF->Close();
}
