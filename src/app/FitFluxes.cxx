#include "DetectorStop.hxx"
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
std::string inpFile, inpHistName;
std::string oupFile;
std::string oupDir;
std::vector<DetectorStop> detStops;
std::vector<std::pair<size_t, size_t> > UsedDetStopSlices;
std::string runPlanCfg;

bool UPDATEOutputFile = false;

std::vector<std::pair<std::string, std::string> > XSecComponentInputs;
std::map<std::string, TH1D *> XSecComponents;
TH1D *TotalXSec;

bool IsGauss = false;
double GaussC, GaussW;
double TargetPeakNorm;

bool UseErrorsInTestStat = true;

double CoeffLimit = 30;
double RegFactor = 0xdeadbeef;

int MaxMINUITCalls = 50000;

double FitBetween_low = 0xdeadbeef, FitBetween_high = 0xdeadbeef;
int binLow, binHigh;

double referencePOT = std::numeric_limits<double>::min();
std::vector<double> MeasurementFactor;
std::vector<TH1D *> Fluxes;

TH1D *SummedFlux;

TH1D *TargetFlux;
TH1D *OscFlux;

TF1 *TargetGauss;

enum OutOfRangeModeEnum { kIgnore = 0, kZero, kExponentialDecay };
/// If using a FitBetween mode:
/// 0: Ignore all bins outside range
/// 1: Try to force bins to 0
/// 2: Exponential decay from target flux at closest kept bin.
int OutOfRangeMode = kZero;
double ExpDecayRate = 3;

void BuildTargetFlux(TH1D *OscFlux) {
  TargetFlux = static_cast<TH1D *>(OscFlux->Clone());
  TargetFlux->SetDirectory(NULL);

  for (Int_t bi_it = 1; bi_it < binLow; ++bi_it) {
    if (OutOfRangeMode) {
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

    TargetFlux->SetBinError(bi_it, TargetFlux->GetBinContent(bi_it) *
                                       (UseErrorsInTestStat ? 0.01 : 0));
  }

  for (Int_t bi_it = binLow; bi_it < (binHigh + 1); ++bi_it) {
    TargetFlux->SetBinContent(bi_it, OscFlux->GetBinContent(bi_it));
    TargetFlux->SetBinError(
        bi_it, UseErrorsInTestStat ? OscFlux->GetBinError(bi_it) : 0);
  }

  for (Int_t bi_it = (binHigh + 1);
       bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
    if (OutOfRangeMode) {
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

    TargetFlux->SetBinError(bi_it, TargetFlux->GetBinContent(bi_it) *
                                       (UseErrorsInTestStat ? 0.01 : 0));
  }
}

void TargetSumChi2(int &nDim, double *gout, double &result, double coeffs[],
                   int flg) {
  SumHistograms(SummedFlux, coeffs, Fluxes);

  double sumdiff = 0;

  static size_t NTargBins = TargetFlux->GetXaxis()->GetNbins();

  for (size_t bi_it = 1; bi_it < NTargBins; ++bi_it) {
    if ((TargetFlux->GetBinContent(bi_it) <
         std::numeric_limits<double>::min()) &&
        (OutOfRangeMode ==
         kZero)) {  // Skip bins if ignoring out of fit range bins.
      continue;
    }
    double err =
        (TargetFlux->GetBinError(bi_it) > std::numeric_limits<double>::min())
            ? (TargetFlux->GetBinError(bi_it) * TargetFlux->GetBinError(bi_it) +
               SummedFlux->GetBinError(bi_it) * SummedFlux->GetBinError(bi_it))
            : 1;

    sumdiff +=
        (TargetFlux->GetBinContent(bi_it) - SummedFlux->GetBinContent(bi_it)) *
        (TargetFlux->GetBinContent(bi_it) - SummedFlux->GetBinContent(bi_it)) /
        err;
  }

  if (RegFactor != 0xdeadbeef) {
    for (size_t i = 0; i < Fluxes.size(); i++) {
      sumdiff += pow((coeffs[i] - coeffs[i + 1]) / RegFactor, 2);
    }
  }

  result = sumdiff;
}

void TargetSumGauss(int &nDim, double *gout, double &result, double coeffs[],
                    int flg) {
  SumHistograms(SummedFlux, coeffs, Fluxes);

  double sumdiff = 0;
  for (Int_t bi_it = binLow; bi_it < binHigh + 1; ++bi_it) {
    double bi_c_E = SummedFlux->GetXaxis()->GetBinCenter(bi_it);

    double GaussEval = TargetGauss->Eval(bi_c_E);
    double SummedBinContent = SummedFlux->GetBinContent(bi_it);
    double Uncert = pow(GaussEval / 20.0, 2) + pow(TargetPeakNorm / 30.0, 2);
    sumdiff += pow(GaussEval - SummedBinContent, 2) / Uncert;
  }

  double reg = 0;
  if (RegFactor != 0xdeadbeef) {
    for (size_t i = 0; i < Fluxes.size(); i++) {
      reg += pow((coeffs[i] - coeffs[i + 1]) / RegFactor, 2);
    }
  }
  sumdiff += reg;

  result = sumdiff;
  std::cout << "[INFO]: Gauss chi2: " << sumdiff << " (reg = " << reg << " )."
            << std::endl;
}

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
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
               "\t-r <RunPlan.XML>                   : An XML file specifying "
               "a run plan  \n"
               "\t                                     to make event rate "
               "predictions and \n"
               "\t                                     measured spectra for. "
               "See          \n"
               "\t                                     documentation for XML "
               "structure.   \n"
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
               "\t-c <CoeffLimit=30>                 : Parameter limits of "
               "flux component \n"
               "\t                                     coefficients.           "
               "           \n"
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
        exit(1);
      }
      GaussC = params[0];
      GaussW = params[1];
      if ((GaussC <= 0) || (GaussW <= 0)) {
        std::cout << "[ERROR]: Recieved " << argv[opt]
                  << " argument for -g, expected \"<GaussMean>,<GaussWidth>\", "
                     "where both values are > 0."
                  << std::endl;
        exit(1);
      }
      IsGauss = true;
      TargetGauss = new TF1("tFunc", "gaus", 0, 20);
    } else if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for l, expected 2." << std::endl;
        exit(1);
      }
      FitBetween_low = params[0];
      FitBetween_high = params[1];
    } else if (std::string(argv[opt]) == "-t") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -i, expected 2." << std::endl;
        exit(1);
      }
      inpFile = params[0];
      inpHistName = params[1];
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
    } else if (std::string(argv[opt]) == "-E") {
      UseErrorsInTestStat = false;
    } else if (std::string(argv[opt]) == "-r") {
      runPlanCfg = argv[++opt];
    } else if (std::string(argv[opt]) == "-x") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() < 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -x, expected 2." << std::endl;
        exit(1);
      }
      for (size_t xs_it = 1; xs_it < params.size(); ++xs_it) {
        XSecComponentInputs.push_back(std::make_pair(params[0], params[xs_it]));
      }
    } else if (std::string(argv[opt]) == "-rg") {
      RegFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
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
    exit(1);
  }

  if (!oupFile.length()) {
    std::cout << "[ERROR]: No output file specified." << std::endl;
    exit(1);
  }

  if (!IsGauss) {
    OscFlux = GetHistogram(inpFile, inpHistName);

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
        exit(1);
      } else {
        std::cout << "[INFO]: Found " << Fluxes.size() << " input fluxes."
                  << std::endl;
      }

      ifl->Close();
      delete ifl;
    } else {
      Fluxes = GetHistograms(FluxesFile, inpFluxHistsPattern);
      if (!Fluxes.size()) {
        std::cout << "[ERROR]: Found no input fluxes matching pattern: \""
                  << inpFluxHistsPattern << "\" in file: \"" << FluxesFile
                  << "\"." << std::endl;
        exit(1);
      } else {
        std::cout << "[INFO]: Found " << Fluxes.size() << " input fluxes."
                  << std::endl;
      }
    }

  } else {
    std::cout << "[ERROR]: Expected either -f or -r options to be passed."
              << std::endl;
    exit(1);
  }

  SummedFlux = static_cast<TH1D *>(Fluxes[0]->Clone());
  SummedFlux->SetDirectory(NULL);
  SummedFlux->Reset();

  if (IsGauss) {
    binLow = 1;
    binHigh = SummedFlux->GetXaxis()->GetNbins();
  }

  TFitter *minimizer = new TFitter(Fluxes.size());

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
  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    std::cout << "[INFO]: fluxpar" << flux_it
              << " uses flux: " << Fluxes[flux_it]->GetName() << std::endl;

    std::stringstream ss("");
    ss << "fluxpar" << flux_it;

    if (abs(flux_it - flux_with_closest_peak_it) ==
        1) {  // Neighbours are negative
      minimizer->SetParameter(flux_it, ss.str().c_str(), -0.35, 0.1,
                              -CoeffLimit, CoeffLimit);
    } else if (abs(flux_it - flux_with_closest_peak_it) == 0) {
      minimizer->SetParameter(flux_it, ss.str().c_str(), 1., 0.1, -CoeffLimit,
                              CoeffLimit);
      TargetPeakNorm = Fluxes[flux_it]->GetMaximum();
    } else {  // Others start free
      minimizer->SetParameter(flux_it, ss.str().c_str(), 0., 0.1, -CoeffLimit,
                              CoeffLimit);
    }
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

  int status = minimizer->ExecuteCommand("MIGRAD", alist, 2);
  if (status) {
    std::cout << "[WARN]: Failed to find minimum (STATUS: " << status << ")."
              << std::endl;
  } else {
    minimizer->PrintResults(1, 0);
  }

  double *coeffs = new double[Fluxes.size()];
  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    coeffs[flux_it] = minimizer->GetParameter(flux_it);
  }

  SumHistograms(SummedFlux, coeffs, Fluxes);

  TFile *oupF =
      new TFile(oupFile.c_str(), UPDATEOutputFile ? "UPDATE" : "RECREATE");
  if (!oupF || !oupF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open output file: " << oupFile << std::endl;
    exit(1);
  }

  TDirectory *oupD = oupF;

  if (oupDir.length()) {
    oupF->mkdir(oupDir.c_str());
    oupF->cd(oupDir.c_str());
    oupD = oupF->GetDirectory(oupDir.c_str());
  }

  oupD->mkdir("weighted_fluxes");
  oupD->cd("weighted_fluxes");
  TDirectory *wD = oupD->GetDirectory("weighted_fluxes");

  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
    cl->Scale(coeffs[flux_it]);
    cl->SetDirectory(wD);
  }

  TH1D *coeffsH = new TH1D("Coeffs", "Coeffs;Off-axis bin;Weight",
                           Fluxes.size(), 0, Fluxes.size());
  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    coeffsH->SetBinContent(flux_it + 1, coeffs[(Fluxes.size() - 1) - flux_it]);
    coeffsH->SetBinError(flux_it + 1,
                         minimizer->GetParError((Fluxes.size() - 1) - flux_it));
  }

  oupD->cd();

  oupD->mkdir("flux_build_unordered");
  oupD->cd("flux_build_unordered");
  wD = oupD->GetDirectory("flux_build_unordered");

  TH1D *cl_0 = static_cast<TH1D *>(Fluxes[0]->Clone());
  cl_0->Scale(coeffs[0]);
  cl_0->SetDirectory(wD);
  std::stringstream ss("flux_0");
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

  oupD->mkdir("flux_build_ordered");
  oupD->cd("flux_build_ordered");
  wD = oupD->GetDirectory("flux_build_ordered");

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
      TH1D *hist = GetHistogram(hdescript.first, hdescript.second);
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

    oupD->mkdir("evr_build_unordered");
    oupD->cd("evr_build_unordered");
    wD = oupD->GetDirectory("evr_build_unordered");

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
  oupF->Write();
  oupF->Close();

  std::cout << "Used " << Fluxes.size() << " fluxes to fit "
            << (binHigh - binLow) << " bins." << std::endl;
}
