#include "GetUsage.hxx"
#include "OscillationHelper.hxx"
#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"

#include "BargerPropagator.h"

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

double DipAngle = 5.8;
double OscParams[6] = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
std::string inpFile, inpHistName;
std::string oupFile, oupHistName;
std::string outputDir;
bool UPDATEOutputFile = false;

int nuPDGFrom, nuPDGTo;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << GetUsageText(argv[0], "flux_tools")
            << std::endl;
}

const static double REarth_cm = 6371.393 * 1.0E5;
const static double ProductionHeight_cm = 0;
static const double deg2rad = asin(1) / 90.0;

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-d") {
      DipAngle = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-p") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 6) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -p, expected 6." << std::endl;
        exit(1);
      }

      OscParams[0] = params[0];
      OscParams[1] = params[1];
      OscParams[2] = params[2];
      OscParams[3] = params[3];
      OscParams[4] = params[4];
      OscParams[5] = params[5];

      std::cout << "Sin^2(Theta_12) = " << OscParams[0] << std::endl;
      std::cout << "Sin^2(Theta_13) = " << OscParams[1] << std::endl;
      std::cout << "Sin^2(Theta_23) = " << OscParams[2] << std::endl;

      std::cout << "Dm^2_21 = " << OscParams[3] << std::endl;
      std::cout << "|Dm^2_Atm| = " << OscParams[4] << std::endl;

      std::cout << "dcp = " << OscParams[5] << std::endl;

    } else if (std::string(argv[opt]) == "-i") {
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
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -o, expected 2." << std::endl;
        exit(1);
      }
      oupFile = params[0];
      oupHistName = params[1];
      UPDATEOutputFile = false;
    } else if (std::string(argv[opt]) == "-a") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -o, expected 2." << std::endl;
        exit(1);
      }
      oupFile = params[0];
      oupHistName = params[1];
      UPDATEOutputFile = true;
    } else if (std::string(argv[opt]) == "-D") {
      outputDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-n") {
      std::vector<int> params = ParseToVect<int>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -n, expected 2." << std::endl;
        exit(1);
      }
      nuPDGFrom = params[0];
      nuPDGTo = params[1];
    } else if (std::string(argv[opt]) == "-L") {
      double baseline_cm = str2T<double>(argv[++opt]) * 1E5;

      DipAngle = asin(baseline_cm / (2.0 * REarth_cm)) / deg2rad;

    } else if (std::string(argv[opt]) == "-?" ||
               std::string(argv[opt]) == "--help") {
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

enum nuTypes {
  kNuebarType = -1,
  kNumubarType = -2,
  kNutaubarType = -3,
  kNueType = 1,
  kNumuType = 2,
  kNutauType = 3,
};

nuTypes GetNuType(int pdg) {
  switch (pdg) {
  case 16:
    return kNutauType;
  case 14:
    return kNumuType;
  case 12:
    return kNueType;
  case -16:
    return kNutaubarType;
  case -14:
    return kNumubarType;
  case -12:
    return kNuebarType;
  default: {
    std::cout << "[ERROR]: Attempting to convert \"neutrino pdg\": " << pdg
              << std::endl;
    exit(1);
  }
  }
}

double OscWeight(double enu) {
  static bool first = true;
  BargerPropagator bp;
  int NuType = GetNuType(nuPDGFrom);
  bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
            OscParams[4], OscParams[5], enu, true, NuType);

  double lengthParam = cos((90.0 + DipAngle) * deg2rad);
  if (first) {
    double PathLength =
        sqrt((REarth_cm + ProductionHeight_cm) *
                 (REarth_cm + ProductionHeight_cm) -
             (REarth_cm * REarth_cm) * (1 - lengthParam * lengthParam)) -
        REarth_cm * lengthParam;

    std::cout << "Calculated path length: " << (PathLength / 1.0E5) << " km."
              << std::endl;
    first = false;
  }

  bp.DefinePath(lengthParam, 0);
  bp.propagate(NuType);
  return bp.GetProb(NuType, GetNuType(nuPDGTo));
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!inpFile.length() || !inpHistName.length()) {
    std::cout << "[ERROR]: No input file or input histogram name specified."
              << std::endl;
    exit(1);
  }

  if (!oupFile.length() || !oupHistName.length()) {
    std::cout << "[ERROR]: No output file or input histogram name specified."
              << std::endl;
    exit(1);
  }

  if (nuPDGTo == 0 || nuPDGFrom == 0) {
    std::cout << "[ERROR]: Must specify-n <nuPDGFrom>,<nuPDGTo> option"
              << std::endl;
    exit(1);
  }

  TFile *inpF = CheckOpenFile(inpFile, "READ");

  TH1 *inpH = GetHistogram<TH1D>(inpF, inpHistName);

  inpH = static_cast<TH1 *>(inpH->Clone());

  TFile *oupF =
      CheckOpenFile(oupFile.c_str(), UPDATEOutputFile ? "UPDATE" : "RECREATE");
  TDirectory *oupD = oupF;

  if (outputDir.length()) {
    oupD = oupF->mkdir(outputDir.c_str());
  }
  oupD->cd();

  inpH->Write();
  inpH->SetDirectory(NULL);

  OscillationHelper oh;
  oh.Setup(OscParams, DipAngle);
  oh.SetOscillationChannel(nuPDGFrom, nuPDGTo);

  for (Int_t bi_it = 1; bi_it < inpH->GetXaxis()->GetNbins() + 1; ++bi_it) {
    double ow = oh.GetWeight(inpH->GetXaxis()->GetBinCenter(bi_it));
    inpH->SetBinContent(bi_it, inpH->GetBinContent(bi_it) * ow);
  }
  inpH->SetName(oupHistName.c_str());
  inpH->SetDirectory(oupD);

  TGraph *POsc = new TGraph();

  POsc->Set(1E4 - 1);

  double min = inpH->GetXaxis()->GetBinLowEdge(1);
  double step = (inpH->GetXaxis()->GetBinUpEdge(inpH->GetXaxis()->GetNbins()) -
                 inpH->GetXaxis()->GetBinLowEdge(1)) /
                double(1E4);
  for (size_t i = 1; i < 1E4; ++i) {
    double enu = min + i * step;
    double ow = oh.GetWeight(enu);
    if (ow != ow) {
      std::cout << "Bad osc weight for ENu: " << enu << std::endl;
    }
    POsc->SetPoint(i - 1, enu, ow);
  }

  if (!UPDATEOutputFile || outputDir.size()) {
    POsc->Write("POsc");
    oh.WriteConfigTree(oupD);
  }

  oupF->Write();
  oupF->Close();
}
