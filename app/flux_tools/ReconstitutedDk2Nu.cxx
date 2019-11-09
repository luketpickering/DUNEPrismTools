#include <iostream>
#include <string>

#include "GetUsage.hxx"
#include "StringParserUtility.hxx"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include "dk2nu_TreeReader.hxx"

#include "dk2nu/tree/dkmeta.h"

std::string inpDescriptor = "";
UInt_t ppfx_NUniverses = 0;
std::string outputFile;
int job = 0;

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
    } else if (std::string(argv[opt]) == "-N") {
      ppfx_NUniverses = str2T<UInt_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-J") {
      job = str2T<int>(argv[++opt]);
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

  handleOpts(argc, argv);
  if (!inpDescriptor.size()) {
    std::cout << "[ERROR]: Expected -i option to be passed." << std::endl;
    SayUsage(argv);
    return 1;
  }

  DK2NuReader *dk2nuRdr =
      new DK2NuReader("dk2nuTree_lite", inpDescriptor, true);
  if (!dk2nuRdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }
  dk2nuRdr->SetPPFXBranchAddresses(ppfx_NUniverses);

  DKMetaReader *dkmetaRdr =
      new DKMetaReader("dkmetaTree_lite", inpDescriptor, true);
  if (!dkmetaRdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }

  TFile *of = new TFile(outputFile.c_str(), "RECREATE");

  TTree *dk2nuTree = new TTree("dk2nuTree", "dk2nu tree");
  bsim::Dk2Nu *dk2nuobj = new bsim::Dk2Nu();
  dk2nuTree->Branch<bsim::Dk2Nu>("dk2nu", dk2nuobj);

  TTree *dkMetaTree = new TTree("dkmetaTree", "dk2nu meta tree");
  bsim::DkMeta *DkMetaobj = new bsim::DkMeta();
  dkMetaTree->Branch<bsim::DkMeta>("dkmeta", DkMetaobj);

  DkMetaobj->job = job;
  DkMetaobj->pots = 0;

  for (size_t it = 0; it < dkmetaRdr->GetEntries(); ++it) {
    dkmetaRdr->GetEntry(it);
    DkMetaobj->pots += dkmetaRdr->pots;
  }

  dkMetaTree->Fill();

  TTree *ppfxWeightFriendTree = nullptr;
  int nppfxu = ppfx_NUniverses;
  Int_t potnum;
  Int_t jobindx;

  if (ppfx_NUniverses) {
    ppfxWeightFriendTree = new TTree("ppfxweights", "ppfxweights");
    ppfxWeightFriendTree->Branch("nppfxu", &nppfxu, "nppfxu/I");
    ppfxWeightFriendTree->Branch("ppfx_cvwgt", &dk2nuRdr->ppfx_cvwgt,
                                 "ppfx_cvwgt/D");
    ppfxWeightFriendTree->Branch("ppfx_vwgt_tot", dk2nuRdr->ppfx_vwgt_tot,
                                 "ppfx_vwgt_tot[nppfxu]/D");
    ppfxWeightFriendTree->Branch("job", &job, "job/I");
    ppfxWeightFriendTree->Branch("potnum", &potnum, "potnum/I");
    ppfxWeightFriendTree->Branch("jobindx", &jobindx, "jobindx/I");
  }

  size_t NEnts = dk2nuRdr->GetEntries();
  size_t NShout = (NEnts / 100);

  for (size_t it = 0; it < dk2nuRdr->GetEntries(); ++it) {
    dk2nuRdr->GetEntry(it);

    if (NShout && !(it % NShout)) {
      std::cout << "\rReconstituted " << it << "/" << NEnts;
    }

    potnum = it;
    jobindx = it;
    *dk2nuobj = dk2nuRdr->GetDk2Nu();
    dk2nuobj->job = job;
    dk2nuobj->potnum = potnum;
    dk2nuobj->jobindx = jobindx;

    dk2nuTree->Fill();

    if (ppfx_NUniverses) {
      ppfxWeightFriendTree->Fill();
    }
  }
  std::cout << "\rReconstituted " << NEnts << "/" << NEnts << std::endl;

  of->Write();
  of->Close();
}
