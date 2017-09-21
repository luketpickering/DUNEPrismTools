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

std::string FluxesFile;
std::string inpFile, inpHistName;
std::string oupFile;
std::string oupDir;
bool UPDATEOutputFile = false;

bool IsGauss = false;
double GaussC, GaussW;
double TargetPeakNorm;

bool UseErrorsInTestStat = true;

/// If using a FitBetween mode:
/// 0: Ignore all bins outside range
/// 1: Try to force bins to 0
/// 2: Exponential decay from target flux at closest kept bin.
int OutOfRangeMode = 1;

double CoeffLimit = 10;

int MaxMINUITCalls = 5000;

double FitBetween_low = 0xdeadbeef, FitBetween_high = 0xdeadbeef;
int binLow, binHigh;

double MaxGFit = 4;

std::vector<TH1D *> Fluxes;

TH1D *SummedFlux;

TH1D *TargetFlux;
TH1D *OscFlux;

TF1 *TargetGauss;

void SumFluxes(double coeffs[]) {
  size_t nf = Fluxes.size();
  SummedFlux->Reset();
  for (size_t f_it = 0; f_it < nf; ++f_it) {
    SummedFlux->Add(Fluxes[f_it], coeffs[f_it]);
  }
}

void BuildTargetFlux(TH1D *OscFlux) {
  TargetFlux = static_cast<TH1D *>(OscFlux->Clone());
  TargetFlux->SetDirectory(NULL);

  for (size_t bi_it = 1; bi_it < binLow; ++bi_it) {
    if (OutOfRangeMode) {
      double target = 0;
      if (OutOfRangeMode == 2) {
        double enu_first_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(binLow);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_first_counted_bin = TargetFlux->GetBinContent(binLow);
        double enu_bottom_bin = TargetFlux->GetXaxis()->GetBinCenter(1);
        double sigma5_range = enu_first_counted_bin - enu_bottom_bin;
        target = content_first_counted_bin *
                 exp(-3.0 * (enu_first_counted_bin - enu) / sigma5_range);
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, TargetFlux->GetBinContent(bi_it) *
                                       (UseErrorsInTestStat ? 0.01 : 0));
    std::cout << "Setting target bin[" << bi_it << " < " << binLow << "] to "
              << TargetFlux->GetBinContent(bi_it)
              << ", E: " << TargetFlux->GetBinError(bi_it) << std::endl;
  }

  for (size_t bi_it = binLow; bi_it < (binHigh + 1); ++bi_it) {
    TargetFlux->SetBinContent(bi_it, OscFlux->GetBinContent(bi_it));
    TargetFlux->SetBinError(
        bi_it, UseErrorsInTestStat ? OscFlux->GetBinError(bi_it) : 0);

    std::cout << "Setting target bin[" << bi_it
              << "] from osc: " << TargetFlux->GetBinContent(bi_it)
              << ", E: " << TargetFlux->GetBinError(bi_it) << std::endl;
  }

  for (size_t bi_it = (binHigh + 1);
       bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
    if (OutOfRangeMode) {
      double target = 0;
      if (OutOfRangeMode == 2) {
        double enu_last_counted_bin =
            TargetFlux->GetXaxis()->GetBinCenter(binHigh);
        double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
        double content_last_counted_bin = TargetFlux->GetBinContent(binHigh);
        double enu_top_bin = TargetFlux->GetXaxis()->GetBinCenter(
            TargetFlux->GetXaxis()->GetNbins());
        double sigma5_range = enu_top_bin - enu_last_counted_bin;
        target = content_last_counted_bin *
                 exp(-3.0 * (enu - enu_last_counted_bin) / sigma5_range);
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, TargetFlux->GetBinContent(bi_it) *
                                       (UseErrorsInTestStat ? 0.01 : 0));

    std::cout << "Setting target bin[" << bi_it << " > " << binHigh << "] to "
              << TargetFlux->GetBinContent(bi_it)
              << ", E: " << TargetFlux->GetBinError(bi_it) << std::endl;
  }
}

void TargetSumChi2(int &nDim, double *gout, double &result, double coeffs[],
                   int flg) {
  SumFluxes(coeffs);

  double sumdiff = 0;

  static size_t NTargBins = TargetFlux->GetXaxis()->GetNbins();

  for (size_t bi_it = 1; bi_it < NTargBins; ++bi_it) {
    if ((TargetFlux->GetBinContent(bi_it) <
         std::numeric_limits<double>::min()) &&
        (OutOfRangeMode ==
         0)) {  // Skip bins if ignoring out of fit range bins.
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

  result = sumdiff;
}

void TargetSumGauss(int &nDim, double *gout, double &result, double coeffs[],
                    int flg) {
  SumFluxes(coeffs);

  double sumdiff = 0;
  for (size_t bi_it = binLow; bi_it < binHigh + 1; ++bi_it) {
    double bi_c_E = SummedFlux->GetXaxis()->GetBinCenter(bi_it);

    // if (bi_c_E > MaxGFit) {
    //   continue;
    // }

    double GaussEval = TargetGauss->Eval(bi_c_E);
    double SummedBinContent = SummedFlux->GetBinContent(bi_it);
    double Uncert = pow(GaussEval / 20.0, 2) + pow(TargetPeakNorm / 30.0, 2);
    sumdiff += pow(GaussEval - SummedBinContent, 2) / Uncert;
  }
  result = sumdiff;
}

std::vector<int> ParseToInt(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<int> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    std::istringstream stream(temp_string);
    int entry;
    stream >> entry;

    vals.push_back(entry);
  }
  return vals;
}

double str2d(std::string str) {
  std::istringstream stream(str);
  double d;
  stream >> d;
  return d;
}

int str2i(std::string str) {
  std::istringstream stream(str);
  int d;
  stream >> d;
  return d;
}

std::vector<double> ParseToDbl(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<double> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    std::istringstream stream(temp_string);
    double entry;
    stream >> entry;

    vals.push_back(entry);
  }
  return vals;
}

std::vector<std::string> ParseToStr(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<std::string> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    vals.push_back(temp_string);
  }

  return vals;
}

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << " -i <input ROOT file>,<input flux hist name> -o/a <output ROOT "
               "file> -d <output_dir> -g <GaussMean>,<GaussWidth> -f <Fluxes "
               "file ROOT> -n "
               "<MaxMINUITCalls> -l <FitEnuRegion_low>,<FitEnuRegion_low> -c "
               "<CoeffLimit> -E (Use uncertainties in fit target flux fit) -m "
               "0 {ignore outside range, force outside range to 0, force "
               "outside range to exp decay}"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-g") {
      std::vector<double> params = ParseToDbl(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for g, expected 2." << std::endl;
        exit(1);
      }
      GaussC = params[0];
      GaussW = params[1];
      if ((GaussC < 1E-8) || (GaussW < 1E-8)) {
        std::cout << "[ERROR]: Recieved " << argv[opt]
                  << " argument for -g, expected \"<GaussMean>,<GaussWidth>\", "
                     "where both values are > 0."
                  << std::endl;
        exit(1);
      }
      IsGauss = true;
      TargetGauss = new TF1("tFunc", "gaus", 0, 20);
    } else if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToDbl(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for l, expected 2." << std::endl;
        exit(1);
      }
      FitBetween_low = params[0];
      FitBetween_high = params[1];
    } else if (std::string(argv[opt]) == "-i") {
      std::vector<std::string> params = ParseToStr(argv[++opt], ",");
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
      FluxesFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-n") {
      MaxMINUITCalls = str2i(argv[++opt]);
    } else if (std::string(argv[opt]) == "-c") {
      CoeffLimit = str2d(argv[++opt]);
    } else if (std::string(argv[opt]) == "-m") {
      OutOfRangeMode = str2i(argv[++opt]);
    } else if (std::string(argv[opt]) == "-E") {
      UseErrorsInTestStat = false;
    } else if (std::string(argv[opt]) == "-gl") {
      MaxGFit = str2d(argv[++opt]);
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
    TFile *inpF = new TFile(inpFile.c_str());
    if (!inpF || !inpF->IsOpen()) {
      std::cout << "[ERROR]: Couldn't open input file: " << inpFile
                << std::endl;
      exit(1);
    }

    OscFlux = dynamic_cast<TH1D *>(inpF->Get(inpHistName.c_str()));

    if (!OscFlux) {
      std::cout << "[ERROR]: Couldn't get TH1D: " << inpHistName
                << " from input file: " << inpFile << std::endl;
      exit(1);
    }

    OscFlux = dynamic_cast<TH1D *>(OscFlux->Clone());
    OscFlux->SetDirectory(NULL);

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

  TFile *fluxF = new TFile(FluxesFile.c_str());
  if (!fluxF || !fluxF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open fluxes file: " << FluxesFile
              << std::endl;
    exit(1);
  }

  TIter next(fluxF->GetListOfKeys());
  TKey *k;
  while ((k = dynamic_cast<TKey *>(next()))) {
    TObject *obj = k->ReadObj();
    if (std::string(obj->ClassName()) != "TH1D") {
      continue;
    }
    if (std::string(obj->GetName()).find("CalcBinning") != std::string::npos) {
      continue;
    }
    std::cout << "Found input flux histogram: " << obj->GetName() << std::endl;
    Fluxes.push_back(static_cast<TH1D *>(obj->Clone()));
    Fluxes.back()->SetDirectory(NULL);
  }

  fluxF->Close();

  if (!Fluxes.size()) {
    std::cout << "[ERROR]: Found no input fluxes." << std::endl;
    exit(1);
  }

  SummedFlux = static_cast<TH1D *>(Fluxes[0]->Clone());
  SummedFlux->SetDirectory(NULL);
  SummedFlux->Reset();

  if (IsGauss) {
    binLow = 1;
    binHigh = SummedFlux->GetXaxis()->GetNbins();
  }

  std::cout << "Have " << Fluxes.size() << " fluxes to work with." << std::endl;

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

  std::cout << "flux_with_closest_peak_it: " << flux_with_closest_peak_it
            << ", flux_peak_dist = " << flux_peak_dist << std::endl;
  for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
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

  double alist[2] = {MaxMINUITCalls, 1E-5 * TargetPeakNorm};

  int status = minimizer->ExecuteCommand("MIGRAD", alist, 2);
  if (status) {
    std::cout << "[WARN]: Failed to find minimum (STATUS: " << status << ")."
              << std::endl;
  } else {
    minimizer->PrintResults(1, 0);
  }

  double *coeffs = new double[Fluxes.size()];
  for (int flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    coeffs[flux_it] = minimizer->GetParameter(flux_it);
  }

  SumFluxes(coeffs);

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

  for (int flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
    cl->Scale(coeffs[flux_it]);
    cl->SetDirectory(wD);
  }

  TH1D *coeffsH = new TH1D("Coeffs", "Coeffs;Off-axis bin;Weight",
                           Fluxes.size(), 0, Fluxes.size());
  for (int flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
    coeffsH->SetBinContent(flux_it + 1,
                           coeffs[(Fluxes.size() - 1) - flux_it]);
    coeffsH->SetBinError(flux_it + 1,
                         minimizer->GetParError((Fluxes.size() - 1) - flux_it));
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

    OscFlux->SetDirectory(oupD);
    OscFlux->SetName("InputFlux");
  }
  oupF->Write();
  oupF->Close();

  std::cout << "Used " << Fluxes.size() << " fluxes to fit "
            << (binHigh - binLow) << " bins." << std::endl;
}
