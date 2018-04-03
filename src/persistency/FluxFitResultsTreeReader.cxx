#include "FluxFitResultsTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

FluxFitResultsTreeReader::FluxFitResultsTreeReader() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
FluxFitResultsTreeReader::FluxFitResultsTreeReader(std::string const &treeName, std::string const &inputFile) : FluxFitResultsTreeReader() {

    tree = OpenTChainWithFileList(treeName, inputFile);

    if(!tree){
      std::cout << "[FluxFitResultsTreeReader]: Failed to read input tree from "
        "file." << std::endl;
        throw;
    }

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[FluxFitResultsTreeReader]: Loaded TChain with " << NEntries
      << " entries." << std::endl;
    GetEntry(0);
  }

  void FluxFitResultsTreeReader::Reset() {
    NFluxes = 0;
    NIterations = 0;
    Chi2 = 0;
    RegularisationPenalty = 0;
    OutOfRangePenalty = 0xdeadbeef;
    std::fill_n(FitRange,2,0);
    NDOverFDFitScaleFactor = 0xdeadbeef;
    IsGaussFit = 0;
    GaussCenter_GeV = 0xdeadbeef;
    GaussWidth_GeV = 0xdeadbeef;
  }

  void FluxFitResultsTreeReader::Copy(FluxFitResultsTreeReader const &other) {
    NFluxes = other.NFluxes;
    NIterations = other.NIterations;
    Chi2 = other.Chi2;
    RegularisationPenalty = other.RegularisationPenalty;
    OutOfRangePenalty = other.OutOfRangePenalty;
    std::copy_n(other.FitRange,2,FitRange);
    NDOverFDFitScaleFactor = other.NDOverFDFitScaleFactor;
    IsGaussFit = other.IsGaussFit;
    GaussCenter_GeV = other.GaussCenter_GeV;
    GaussWidth_GeV = other.GaussWidth_GeV;
  }

  void FluxFitResultsTreeReader::SetBranchAddresses() {
    tree->SetBranchAddress("NFluxes", &NFluxes);
    tree->SetBranchAddress("NIterations", &NIterations);
    tree->SetBranchAddress("Chi2", &Chi2);
    tree->SetBranchAddress("RegularisationPenalty", &RegularisationPenalty);
    tree->SetBranchAddress("FitRange", &FitRange);

    IsGaussFit = CheckTTreeHasBranch(tree, "GaussCenter_GeV");
    if(IsGaussFit){
      tree->SetBranchAddress("GaussCenter_GeV", &GaussCenter_GeV);
      tree->SetBranchAddress("GaussWidth_GeV", &GaussWidth_GeV);
    } else {
      tree->SetBranchAddress("OutOfRangePenalty", &OutOfRangePenalty);
      tree->SetBranchAddress("NDOverFDFitScaleFactor", &NDOverFDFitScaleFactor);
    }

  }

  void FluxFitResultsTreeReader::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t FluxFitResultsTreeReader::GetEntry() { return CEnt; }
  UInt_t FluxFitResultsTreeReader::GetEntries() { return NEntries; }

  FluxFitResultsTreeReader *FluxFitResultsTreeReader::MakeTreeWriter(TTree *tree, bool IsGaussFit) {
    FluxFitResultsTreeReader *fdr = new FluxFitResultsTreeReader();
    tree->Branch("NFluxes", &fdr->NFluxes,"NFluxes/I");
    tree->Branch("NIterations", &fdr->NIterations,"NIterations/I");
    tree->Branch("Chi2", &fdr->Chi2,"Chi2/D");
    tree->Branch("RegularisationPenalty", &fdr->RegularisationPenalty,
      "RegularisationPenalty/D");
    tree->Branch("FitRange", &fdr->FitRange,"FitRange[2]/D");

    if(IsGaussFit){
      tree->Branch("GaussCenter_GeV", &fdr->GaussCenter_GeV,
        "GaussCenter_GeV/D");
      tree->Branch("GaussWidth_GeV", &fdr->GaussWidth_GeV,"GaussWidth_GeV/D");
    } else {
      tree->Branch("OutOfRangePenalty", &fdr->OutOfRangePenalty,
        "OutOfRangePenalty/D");
      tree->Branch("NDOverFDFitScaleFactor", &fdr->NDOverFDFitScaleFactor,
        "NDOverFDFitScaleFactor/D");
    }

    fdr->Reset();
    return fdr;
  }

  void FluxFitResultsTreeReader::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  FluxFitResultsTreeReader::~FluxFitResultsTreeReader() {
    delete tree;
  }
