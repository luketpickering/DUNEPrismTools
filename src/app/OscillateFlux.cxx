#include "Utils.hxx"

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

int nuPDGFrom, nuPDGTo;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << " -p "
               "<sin2(theta12)=0.825>,<sin2(theta13)=0.10>,<sin2(theta23)=1.0>,"
               "<dm12=7.9e-5>,<"
               "dm23=2.5e-3>,<dcp=0.0> -d <dipangle=5.8> -i <input ROOT "
               "file>,<input flux hist name> -o <output ROOT file>,<output "
               "flux hist name> -n <nuPDGFrom>,<nuPDGTo>"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-d") {
      DipAngle = str2T<double>(argv[++opt]);

      const static double REarth_cm = 6371.393 * 1.0E5;
      const static double ProductionHeight_cm = 0;
      static const double deg2rad = asin(1) / 90.0;
      double cz = cos((90.0 + DipAngle) * deg2rad);
      double PathLength = sqrt((REarth_cm + ProductionHeight_cm) *
                                   (REarth_cm + ProductionHeight_cm) -
                               (REarth_cm * REarth_cm) * (1 - cz * cz)) -
                          REarth_cm * cz;

      std::cout << "Calculated path length: " << (PathLength / 1.0E5) << " km."
                << std::endl;

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
    } else if (std::string(argv[opt]) == "-n") {
      std::vector<int> params = ParseToVect<int>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -n, expected 2." << std::endl;
        exit(1);
      }
      nuPDGFrom = params[0];
      nuPDGTo = params[1];
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
  BargerPropagator bp;
  int NuType = GetNuType(nuPDGFrom);
  bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
            OscParams[4], OscParams[5], enu, true, NuType);

  static const double deg2rad = asin(1) / 90.0;
  double lengthParam = cos((90.0 + DipAngle) * deg2rad);
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

  TFile *inpF = new TFile(inpFile.c_str());
  if (!inpF || !inpF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open input file: " << inpFile << std::endl;
    exit(1);
  }

  TH1D *inpH = dynamic_cast<TH1D *>(inpF->Get(inpHistName.c_str()));

  if (!inpH) {
    std::cout << "[ERROR]: Couldn't get TH1D: " << inpHistName
              << " from input file: " << inpFile << std::endl;
    exit(1);
  }

  inpH = static_cast<TH1D *>(inpH->Clone());

  TFile *oupF = new TFile(oupFile.c_str(), "RECREATE");
  if (!oupF || !oupF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open output file: " << oupFile << std::endl;
    exit(1);
  }

  inpH->Write();
  inpH->SetDirectory(NULL);

  for (Int_t bi_it = 1; bi_it < inpH->GetXaxis()->GetNbins() + 1; ++bi_it) {
    double ow = OscWeight(inpH->GetXaxis()->GetBinCenter(bi_it));
    inpH->SetBinContent(bi_it, inpH->GetBinContent(bi_it) * ow);
  }
  inpH->SetName(oupHistName.c_str());
  inpH->SetDirectory(oupF);

  TGraph *POsc = new TGraph();

  POsc->Set(1E4 - 1);

  double min = inpH->GetXaxis()->GetBinLowEdge(1);
  double step = (inpH->GetXaxis()->GetBinUpEdge(inpH->GetXaxis()->GetNbins()) -
                 inpH->GetXaxis()->GetBinLowEdge(1)) /
                double(1E4);
  for (size_t i = 1; i < 1E4; ++i) {
    double enu = min + i * step;
    double ow = OscWeight(enu);
    if (ow != ow) {
      std::cout << "Bad osc weight for ENu: " << enu << std::endl;
    }
    POsc->SetPoint(i - 1, enu, ow);
  }

  POsc->Write("POsc");
  oupF->Write();
  oupF->Close();
}
