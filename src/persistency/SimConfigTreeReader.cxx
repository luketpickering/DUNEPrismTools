#include "SimConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

SimConfig::SimConfig() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
SimConfig::SimConfig(std::string const &treeName, std::string const &inputFile) {

    tree = OpenTChainWithFileList(treeName, inputFile);

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[SimConfig]: Loaded TChain: " << NFiles
              << " files and " << NEntries << " entries." << std::endl;
    GetEntry(0);
  }

  void SimConfig::Reset() {
    NXSteps = 0;
    std::fill_n(DetMin,3,0);
    std::fill_n(DetMax,3,0);
    std::fill_n(FVGap,3,0);
    NMaxTrackSteps = 0;
    POTPerFile = 0xdeadbeef;
    timesep_us = 0xdeadbeef;
  }

  void SimConfig::Copy(SimConfig const &other) {
    NXSteps = other.NXSteps;
    std::copy_n(other.DetMin,3,DetMin);
    std::copy_n(other.DetMax,3,DetMax);
    std::copy_n(other.FVGap,3,FVGap);
    NMaxTrackSteps = other.NMaxTrackSteps;
    POTPerFile = other.POTPerFile;
    timesep_us = other.timesep_us;
  }

  void SimConfig::SetBranchAddresses() {
    tree->SetBranchAddress("NXSteps", &NXSteps);
    tree->SetBranchAddress("DetMin", &DetMin);
    tree->SetBranchAddress("DetMax", &DetMax);
    tree->SetBranchAddress("FVGap", &FVGap);
    tree->SetBranchAddress("NMaxTrackSteps", &NMaxTrackSteps);
    tree->SetBranchAddress("POTPerFile", &POTPerFile);
    tree->SetBranchAddress("timesep_us", &timesep_us);
  }

  void SimConfig::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t SimConfig::GetEntry() { return CEnt; }
  UInt_t SimConfig::GetEntries() { return NEntries; }

  SimConfig *SimConfig::MakeTreeWriter(TTree *tree) {
    SimConfig *fdr = new SimConfig();
    tree->Branch("NXSteps", &fdr->NXSteps, "NXSteps/I");
    tree->Branch("DetMin", &fdr->DetMin, "DetMin[3]/D");
    tree->Branch("DetMax", &fdr->DetMax, "DetMax[3]/D");
    tree->Branch("FVGap", &fdr->FVGap, "FVGap[3]/D");
    tree->Branch("NMaxTrackSteps", &fdr->NMaxTrackSteps, "NMaxTrackSteps/I");
    tree->Branch("POTPerFile", &fdr->POTPerFile, "POTPerFile/D");
    tree->Branch("timesep_us", &fdr->timesep_us, "timesep_us/D");
    return fdr;
  }

  void SimConfig::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  SimConfig::~SimConfig() {
    delete tree;
  }
