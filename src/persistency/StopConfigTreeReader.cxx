#include "StopConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

StopConfig::StopConfig() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
StopConfig::StopConfig(std::string const &treeName, std::string const &inputFile) {

    tree = OpenTChainWithFileList(treeName, inputFile);

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[StopConfig]: Loaded TChain: " << NFiles
              << " files and " << NEntries << " entries." << std::endl;
    GetEntry(0);
  }

  void StopConfig::Reset() {
    std::fill_n(ActiveMin,3,0);
  std::fill_n(ActiveMax,3,0);
  std::fill_n(VetoGap,3,0);
  std::fill_n(CenterPosition,3,0);
  POTExposure = 0xdeadbeef;
  }

  void StopConfig::Copy(StopConfig const &other) {
    std::copy_n(other.ActiveMin,3,ActiveMin);
    std::copy_n(other.ActiveMax,3,ActiveMax);
    std::copy_n(other.VetoGap,3,VetoGap);
    std::copy_n(other.CenterPosition,3,CenterPosition);
    POTExposure = other.POTExposure;
  }

  void StopConfig::SetBranchAddresses() {
    tree->SetBranchAddress("ActiveMin", &ActiveMin);
    tree->SetBranchAddress("ActiveMax", &ActiveMax);
    tree->SetBranchAddress("VetoGap", &VetoGap);
    tree->SetBranchAddress("CenterPosition", &CenterPosition);
    tree->SetBranchAddress("POTExposure", &POTExposure);
  }

  void StopConfig::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t StopConfig::GetEntry() { return CEnt; }
  UInt_t StopConfig::GetEntries() { return NEntries; }

  std::vector<BoundingBox> GetStopBoundingBoxes(bool RemoveVeto,
    std::array<double,3> FVReduction){
    if(!tree){
      std::cout << "[ERROR]: Attempted to get stop bounding boxes from "
        "StopConfig, but tree not initialized." << std::endl;
      throw;
    }
    std::vector<BoundingBox> BBs;
    std::cout << "[INFO]: Using reduced fiducial volumes:"  << std::endl;
    for(Long64_t cs_it = 0; cs_it < GetEntries(); ++cs_it){
      GetEntry(cs_it);
      std::array<Double_t,3> FVMin, FVMax;
      for(size_t dim_it = 0; dim_it < 3; ++dim_it){
        FVMin[dim_it] = ActiveMin[dim_it] + (RemoveVeto?VetoGap[dim_it]:0) +
          FVReduction[dim_it];
        FVMax[dim_it] = ActiveMax[dim_it] - (RemoveVeto?VetoGap[dim_it]:0) -
          FVReduction[dim_it];
      }
      BBs.emplace_back(FVMax, FVMin);
    }
    return BBs;
  }


  StopConfig *StopConfig::MakeTreeWriter(TTree *tree) {
    StopConfig *fdr = new StopConfig();
    tree->Branch("ActiveMin", &fdr->ActiveMin,"ActiveMin[3]/D");
    tree->Branch("ActiveMax", &fdr->ActiveMax,"ActiveMax[3]/D");
    tree->Branch("VetoGap", &fdr->VetoGap,"VetoGap[3]/D");
    tree->Branch("CenterPosition", &fdr->CenterPosition,"CenterPosition[3]/D");
    tree->Branch("POTExposure", &fdr->POTExposure,"POTExposure/D");
    fdr->Reset();
    return fdr;
  }

  void StopConfig::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  StopConfig::~StopConfig() {
    delete tree;
  }
