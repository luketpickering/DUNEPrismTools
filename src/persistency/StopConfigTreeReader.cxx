#include "StopConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

  StopConfig::StopConfig(std::string const &inputFile) {

    LoadTree(inputFile);
  }

  std::string StopConfig::TreeName(){
    return "StopConfigTree";
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

  void StopConfig::DetermineNStops(){
    if(!TreeOwned){
      return;
    }
    TChain * chain = dynamic_cast<TChain*>(tree);
    if(!chain){
      std::cout << "[ERROR]: Failed to cast input TTree to TChain in StopConfig"
        " when reading inputs." << std::endl;
      throw;
    }
    UInt_t FFUUID = std::numeric_limits<UInt_t>::max();
    NStops = 0;
    for(Long64_t cs_it = 0; cs_it < GetEntries(); ++cs_it){
      GetEntry(cs_it);
      if(FFUUID == std::numeric_limits<UInt_t>::max()){
        FFUUID = chain->GetFile()->GetUniqueID();
        NStops++;
      } else {
        // Only read the first file's worth.
        if(FFUUID != chain->GetFile()->GetUniqueID()){
          break;
        } else {
          NStops++;
        }
      }
    }
  }

  void StopConfig::SetBranchAddresses() {
    tree->SetBranchAddress("ActiveMin", &ActiveMin);
    tree->SetBranchAddress("ActiveMax", &ActiveMax);
    tree->SetBranchAddress("VetoGap", &VetoGap);
    tree->SetBranchAddress("CenterPosition", &CenterPosition);
    tree->SetBranchAddress("POTExposure", &POTExposure);
  }

  std::vector<BoundingBox> StopConfig::GetStopBoundingBoxes(bool RemoveVeto,
    std::array<double,3> FVReduction){
    if(!tree){
      std::cout << "[ERROR]: Attempted to get stop bounding boxes from "
        "StopConfig, but tree not initialized." << std::endl;
      throw;
    }
    std::vector<BoundingBox> BBs;
    std::cout << "[INFO]: Using stop volumes:"  << std::endl;
    DetermineNStops();
    for(Long64_t cs_it = 0; cs_it < NStops; ++cs_it){
      GetEntry(cs_it);
      std::array<Double_t,3> FVMin, FVMax;
      for(size_t dim_it = 0; dim_it < 3; ++dim_it){
        FVMin[dim_it] = ActiveMin[dim_it] + (RemoveVeto?VetoGap[dim_it]:0) +
          FVReduction[dim_it];
        FVMax[dim_it] = ActiveMax[dim_it] - (RemoveVeto?VetoGap[dim_it]:0) -
          FVReduction[dim_it];
      }
      BBs.emplace_back(FVMax, FVMin);
      std::cout << "\t[" << cs_it << "]" << BBs.back().Print() << std::endl;
    }
    return BBs;
  }


  StopConfig *StopConfig::MakeTreeWriter() {
    StopConfig *fdr = new StopConfig();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->tree->Branch("ActiveMin", &fdr->ActiveMin,"ActiveMin[3]/D");
    fdr->tree->Branch("ActiveMax", &fdr->ActiveMax,"ActiveMax[3]/D");
    fdr->tree->Branch("VetoGap", &fdr->VetoGap,"VetoGap[3]/D");
    fdr->tree->Branch("CenterPosition", &fdr->CenterPosition,"CenterPosition[3]/D");
    fdr->tree->Branch("POTExposure", &fdr->POTExposure,"POTExposure/D");
    fdr->Reset();
    return fdr;
  }

  StopConfig::~StopConfig() {}
