#include "SelectionSummaryTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

SelectionSummary::SelectionSummary(std::string const &inputFile) {

    LoadTree(inputFile);
  }

  std::string SelectionSummary::TreeName(){
    return "SelectionSummaryTree";
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
    NOOAcceptance = 0;
    NOOFV = 0;
    NSel = 0;
    SelectOnMuonExit = 0;
    MuonExitKECut_MeV = 0;
    HadronicShowerVetoCut_MeV = 0;
    std::fill_n(VertexSelectionFV,3,0);
    TotalPOT = 0;
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
    NOOAcceptance = other.NOOAcceptance;
    NOOFV = other.NOOFV;
    NSel  = other.NSel;
    SelectOnMuonExit  = other.SelectOnMuonExit;
    MuonExitKECut_MeV  = other.MuonExitKECut_MeV;
    HadronicShowerVetoCut_MeV  = other.HadronicShowerVetoCut_MeV;
    std::copy_n(other.VertexSelectionFV,3,VertexSelectionFV);
    TotalPOT = other.TotalPOT;
  }

  void SelectionSummary::SetBranchAddresses() {
    tree->SetBranchAddress("NTotal", &NTotal);
    tree->SetBranchAddress("NInStops", &NInStops);
    tree->SetBranchAddress("NNumuCC", &NNumuCC);
    tree->SetBranchAddress("NNumuNC", &NNumuNC);
    tree->SetBranchAddress("NNumubCC", &NNumubCC);
    tree->SetBranchAddress("NNumubNC", &NNumubNC);
    tree->SetBranchAddress("NNueCC", &NNueCC);
    tree->SetBranchAddress("NNueNC", &NNueNC);
    tree->SetBranchAddress("NNuebCC", &NNuebCC);
    tree->SetBranchAddress("NNuebNC", &NNuebNC);
    tree->SetBranchAddress("NOOAcceptance", &NOOAcceptance);
    tree->SetBranchAddress("NOOFV", &NOOFV);
    tree->SetBranchAddress("NSel", &NSel);
    tree->SetBranchAddress("SelectOnMuonExit", &SelectOnMuonExit);
    tree->SetBranchAddress("MuonExitKECut_MeV", &MuonExitKECut_MeV);
    tree->SetBranchAddress("HadronicShowerVetoCut_MeV",
    &HadronicShowerVetoCut_MeV);
    tree->SetBranchAddress("VertexSelectionFV", &VertexSelectionFV);
    tree->SetBranchAddress("TotalPOT", &TotalPOT);

  }

  SelectionSummary *SelectionSummary::MakeTreeWriter() {
    SelectionSummary *fdr = new SelectionSummary();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->tree->Branch("NTotal", &fdr->NTotal, "NTotal/I");
    fdr->tree->Branch("NInStops", &fdr->NInStops, "NInStops/I");
    fdr->tree->Branch("NNumuCC", &fdr->NNumuCC, "NNumuCC/I");
    fdr->tree->Branch("NNumuNC", &fdr->NNumuNC, "NNumuNC/I");
    fdr->tree->Branch("NNumubCC", &fdr->NNumubCC, "NNumubCC/I");
    fdr->tree->Branch("NNumubNC", &fdr->NNumubNC, "NNumubNC/I");
    fdr->tree->Branch("NNueCC", &fdr->NNueCC, "NNueCC/I");
    fdr->tree->Branch("NNueNC", &fdr->NNueNC, "NNueNC/I");
    fdr->tree->Branch("NNuebCC", &fdr->NNuebCC, "NNuebCC/I");
    fdr->tree->Branch("NNuebNC", &fdr->NNuebNC, "NNuebNC/I");
    fdr->tree->Branch("NOOAcceptance", &fdr->NOOAcceptance, "NOOAcceptance/I");
    fdr->tree->Branch("NOOFV", &fdr->NOOFV, "NOOFV/I");
    fdr->tree->Branch("NSel", &fdr->NSel, "NSel/I");
    fdr->tree->Branch("SelectOnMuonExit", &fdr->SelectOnMuonExit,
    "SelectOnMuonExit/O");
    fdr->tree->Branch("MuonExitKECut_MeV", &fdr->MuonExitKECut_MeV,
    "MuonExitKECut_MeV/D");
    fdr->tree->Branch("HadronicShowerVetoCut_MeV",
    &fdr->HadronicShowerVetoCut_MeV, "HadronicShowerVetoCut_MeV/D");
    fdr->tree->Branch("VertexSelectionFV", &fdr->VertexSelectionFV,
    "VertexSelectionFV[3]/D");
    fdr->tree->Branch("TotalPOT", &fdr->TotalPOT, "TotalPOT/D");
    fdr->Reset();
    return fdr;
  }

  SelectionSummary::~SelectionSummary() {}
