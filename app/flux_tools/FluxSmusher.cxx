#include "GetUsage.hxx"
#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <string>
#include <vector>

struct HistDescriptor {
  std::string file;
  std::string histname;
};

std::vector<HistDescriptor> Input1DFluxes;
HistDescriptor Input2DFlux;
HistDescriptor OutputFlux;

bool UPDATEOutputFile = false;
std::vector<std::pair<double, double>> XRanges;
std::vector<TH1D *> Fluxes;
TH1D *SummedFlux = nullptr;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << GetUsageText(argv[0], "flux_tools")
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-i") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -i, expected 2." << std::endl;
        exit(1);
      }
      Input1DFluxes.push_back({params[0], params[1]});
    } else if (std::string(argv[opt]) == "-f") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -f, expected 2." << std::endl;
        throw;
      }
      Input2DFlux.file = params[0];
      Input2DFlux.histname = params[1];
    } else if (std::string(argv[opt]) == "-o") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -f, expected 2." << std::endl;
        throw;
      }
      OutputFlux.file = params[0];
      OutputFlux.histname = params[1];
      UPDATEOutputFile = false;
    } else if (std::string(argv[opt]) == "-a") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -f, expected 2." << std::endl;
        throw;
      }
      OutputFlux.file = params[0];
      OutputFlux.histname = params[1];
      UPDATEOutputFile = true;
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

  if ((!OutputFlux.file.length()) || (!OutputFlux.histname.length())) {
    std::cout << "[ERROR]: No output flux file or histogram name specified."
              << std::endl;
    throw;
  }

  if (!Input1DFluxes.size() &&
      ((!Input2DFlux.file.length()) || (!Input2DFlux.histname.length()))) {
    std::cout << "[ERROR]: No input flux file or histogram name specified."
              << std::endl;
    throw;
  }

  if (Input1DFluxes.size()) {
    for (HistDescriptor const &hd : Input1DFluxes) {
      Fluxes.push_back(GetHistogram<TH1D>(hd.file, hd.histname));
    }
  } else if (XRanges.size()) {
    Fluxes = MergeSplitTH2D(
        GetHistogram<TH2D>(Input2DFlux.file, Input2DFlux.histname), true,
        XRanges);
  } else {
    std::cout << "[ERROR]: Invalid configuration." << std::endl;
    exit(1);
  }

  for (TH1D *flux : Fluxes) {
    if (!SummedFlux) {
      SummedFlux =
          dynamic_cast<TH1D *>(flux->Clone(OutputFlux.histname.c_str()));
    } else {
      SummedFlux->Add(flux);
    }
  }

  TFile *oupF =
      CheckOpenFile(OutputFlux.file, UPDATEOutputFile ? "UPDATE" : "RECREATE");
  SummedFlux->SetDirectory(oupF);
  oupF->Write();
  oupF->Close();
}
