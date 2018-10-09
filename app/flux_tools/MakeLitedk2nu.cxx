#include <iostream>
#include <string>

#include "GetUsage.hxx"
#include "StringParserUtility.hxx"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include "dk2nu_TreeReader.hxx"

std::string inpDescriptor = "";
std::string ppfx_inpDescriptor = "";
UInt_t ppfx_NUniverses = 100;
std::string outputFile;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << GetUsageText(argv[0], "flux_tools")
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpDescriptor = argv[++opt];
    } else if (std::string(argv[opt]) == "-p") {
      ppfx_inpDescriptor = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      ppfx_NUniverses = str2T<UInt_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
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
  if (!inpDescriptor.size()) {
    std::cout << "[ERROR]: Expected -i option to be passed." << std::endl;
    SayUsage(argv);
    return 1;
  }

  DK2NuReader *dk2nuRdr = new DK2NuReader("dk2nuTree", inpDescriptor, false);
  if (!dk2nuRdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }

  if (ppfx_inpDescriptor.size()) {
    dk2nuRdr->AddPPFXFriend("PPFX_weights", ppfx_inpDescriptor,
                            ppfx_NUniverses);
  }

  DKMetaReader *dkmRdr = new DKMetaReader("dkmetaTree", inpDescriptor, false);

  TFile *of = new TFile(outputFile.c_str(), "RECREATE");

  TTree *dk2numetalite = new TTree("dkmetaTree_lite", "dk2nu lite meta tree");

  std::cout << "[INFO]: Cloning " << dkmRdr->GetEntries()
            << " meta tree entries." << std::endl;

  dkmRdr->WriteOutLiteTree(dk2numetalite);

  std::cout << "[INFO]: Done." << std::endl;
  dk2numetalite->Write();

  TTree *dk2nulite = new TTree("dk2nuTree_lite", "dk2nu lite tree");

  std::cout << "[INFO]: Cloning " << dk2nuRdr->GetEntries()
            << " dk2nu tree entries." << std::endl;
  dk2nuRdr->WriteOutLiteTree(dk2nulite);

  std::cout << "[INFO]: Done." << std::endl;

  dk2nulite->Write();

  of->Write();
  of->Close();
}
