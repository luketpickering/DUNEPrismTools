#include "FluxFitResultsTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

  FluxFitResultsTreeReader::FluxFitResultsTreeReader(
    std::string const &inputFile) {

    LoadTree(inputFile);
  }

  std::string FluxFitResultsTreeReader::TreeName(){
    return "FluxFitResultsTree";
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

  FluxFitResultsTreeReader *FluxFitResultsTreeReader::MakeTreeWriter(
    bool IsGaussFit) {

    FluxFitResultsTreeReader *fdr = new FluxFitResultsTreeReader();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->tree->Branch("NFluxes", &fdr->NFluxes,"NFluxes/I");
    fdr->tree->Branch("NIterations", &fdr->NIterations,"NIterations/I");
    fdr->tree->Branch("Chi2", &fdr->Chi2,"Chi2/D");
    fdr->tree->Branch("RegularisationPenalty", &fdr->RegularisationPenalty,
      "RegularisationPenalty/D");
    fdr->tree->Branch("FitRange", &fdr->FitRange,"FitRange[2]/D");

    if(IsGaussFit){
      fdr->tree->Branch("GaussCenter_GeV", &fdr->GaussCenter_GeV,
        "GaussCenter_GeV/D");
      fdr->tree->Branch("GaussWidth_GeV", &fdr->GaussWidth_GeV,"GaussWidth_GeV/D");
    } else {
      fdr->tree->Branch("OutOfRangePenalty", &fdr->OutOfRangePenalty,
        "OutOfRangePenalty/D");
      fdr->tree->Branch("NDOverFDFitScaleFactor", &fdr->NDOverFDFitScaleFactor,
        "NDOverFDFitScaleFactor/D");
    }

    fdr->Reset();
    return fdr;
  }

  FluxFitResultsTreeReader::~FluxFitResultsTreeReader() {}
