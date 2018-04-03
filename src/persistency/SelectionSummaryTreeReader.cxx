#include "SelectionSummaryTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

SelectionSummary::SelectionSummary() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
SelectionSummary::SelectionSummary(std::string const &treeName, std::string const &inputFile) : SelectionSummary() {

    tree = OpenTChainWithFileList(treeName, inputFile);
    if(!tree){
      std::cout << "[SelectionSummary]: Failed to read input tree from file."
        << std::endl;
        throw;
    }

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[SelectionSummary]: Loaded TChain with " << NEntries
      << " entries." << std::endl;
    GetEntry(0);
  }

  void SelectionSummary::Reset() {
    NTotal = 0;
    NInStops = 0;
    NNumuCC = 0;
    NNumuNC = 0;
    NNumubCC = 0;
    NNumubNC = 0;
    NNueCC = 0;
    NNueNC = 0;
    NNuebCC = 0;
    NNuebNC = 0;
    NOOFV = 0;
    NSel = 0;
    SelectOnMuonExit = 0;
    MuonExitKECut_MeV = 0;
    HadronicShowerVetoCut_MeV = 0;
    std::fill_n(VertexSelectionFV,3,0);
  }

  void SelectionSummary::Copy(SelectionSummary const &other) {
    NTotal  = other.NTotal;
    NInStops  = other.NInStops;
    NNumuCC  = other.NNumuCC;
    NNumuNC  = other.NNumuNC;
    NNumubCC  = other.NNumubCC;
    NNumubNC  = other.NNumubNC;
    NNueCC  = other.NNueCC;
    NNueNC  = other.NNueNC;
    NNuebCC  = other.NNuebCC;
    NNuebNC  = other.NNuebNC;
    NOOFV = other.NOOFV;
    NSel  = other.NSel;
    SelectOnMuonExit  = other.SelectOnMuonExit;
    MuonExitKECut_MeV  = other.MuonExitKECut_MeV;
    HadronicShowerVetoCut_MeV  = other.HadronicShowerVetoCut_MeV;
    std::copy_n(other.VertexSelectionFV,3,VertexSelectionFV);
  }

  void SelectionSummary::SetBranchAddresses() {
    tree->Branch("NTotal", &NTotal);
    tree->Branch("NInStops", &NInStops);
    tree->Branch("NNumuCC", &NNumuCC);
    tree->Branch("NNumuNC", &NNumuNC);
    tree->Branch("NNumubCC", &NNumubCC);
    tree->Branch("NNumubNC", &NNumubNC);
    tree->Branch("NNueCC", &NNueCC);
    tree->Branch("NNueNC", &NNueNC);
    tree->Branch("NNuebCC", &NNuebCC);
    tree->Branch("NNuebNC", &NNuebNC);
    tree->Branch("NOOFV", &NOOFV);
    tree->Branch("NSel", &NSel);
    tree->Branch("SelectOnMuonExit", &SelectOnMuonExit);
    tree->Branch("MuonExitKECut_MeV", &MuonExitKECut_MeV);
    tree->Branch("HadronicShowerVetoCut_MeV",
    &HadronicShowerVetoCut_MeV);
    tree->Branch("VertexSelectionFV", &VertexSelectionFV);

  }

  void SelectionSummary::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t SelectionSummary::GetEntry() { return CEnt; }
  UInt_t SelectionSummary::GetEntries() { return NEntries; }

  SelectionSummary *SelectionSummary::MakeTreeWriter(TTree *tree) {
    SelectionSummary *fdr = new SelectionSummary();
    tree->Branch("NTotal", &fdr->NTotal, "NTotal/I");
    tree->Branch("NInStops", &fdr->NInStops, "NInStops/I");
    tree->Branch("NNumuCC", &fdr->NNumuCC, "NNumuCC/I");
    tree->Branch("NNumuNC", &fdr->NNumuNC, "NNumuNC/I");
    tree->Branch("NNumubCC", &fdr->NNumubCC, "NNumubCC/I");
    tree->Branch("NNumubNC", &fdr->NNumubNC, "NNumubNC/I");
    tree->Branch("NNueCC", &fdr->NNueCC, "NNueCC/I");
    tree->Branch("NNueNC", &fdr->NNueNC, "NNueNC/I");
    tree->Branch("NNuebCC", &fdr->NNuebCC, "NNuebCC/I");
    tree->Branch("NNuebNC", &fdr->NNuebNC, "NNuebNC/I");
    tree->Branch("NOOFV", &fdr->NOOFV, "NOOFV/I");
    tree->Branch("NSel", &fdr->NSel, "NSel/I");
    tree->Branch("SelectOnMuonExit", &fdr->SelectOnMuonExit,
    "SelectOnMuonExit/O");
    tree->Branch("MuonExitKECut_MeV", &fdr->MuonExitKECut_MeV,
    "MuonExitKECut_MeV/D");
    tree->Branch("HadronicShowerVetoCut_MeV",
    &fdr->HadronicShowerVetoCut_MeV, "HadronicShowerVetoCut_MeV/D");
    tree->Branch("VertexSelectionFV", &fdr->VertexSelectionFV,
    "VertexSelectionFV[3]/D");
    fdr->Reset();
    return fdr;
  }

  void SelectionSummary::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  SelectionSummary::~SelectionSummary() {
    delete tree;
  }
