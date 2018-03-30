#include "OscillationParametersTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

OscillationParameters::OscillationParameters() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
OscillationParameters::OscillationParameters(std::string const &treeName, std::string const &inputFile) {

    tree = OpenTChainWithFileList(treeName, inputFile);

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[OscillationParameters]: Loaded TChain: " << NFiles
              << " files and " << NEntries << " entries." << std::endl;
    GetEntry(0);
  }

  void OscillationParameters::Reset() {
    DipAngle_degrees = 0xdeadbeef;
    std::fill_n(OscParams,6,0xdeadbeef);
    FromNuPDG = 0;
    ToNuPDG = 0;
  }

  void OscillationParameters::Copy(OscillationParameters const &other) {
    DipAngle_degrees = other.DipAngle_degrees;
    std::copy_n(other.OscParams,6,OscParams);
    FromNuPDG = other.FromNuPDG;
    ToNuPDG = other.ToNuPDG;
  }

  void OscillationParameters::SetBranchAddresses() {
    tree->SetBranchAddress("DipAngle_degrees", &DipAngle_degrees);
    tree->SetBranchAddress("OscParams", &OscParams);
    tree->SetBranchAddress("FromNuPDG", &FromNuPDG);
    tree->SetBranchAddress("ToNuPDG", &ToNuPDG);
  }

  void OscillationParameters::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t OscillationParameters::GetEntry() { return CEnt; }
  UInt_t OscillationParameters::GetEntries() { return NEntries; }

  OscillationParameters *OscillationParameters::MakeTreeWriter(TTree *tree) {
    OscillationParameters *fdr = new OscillationParameters();
    tree->Branch("DipAngle_degrees", &fdr->DipAngle_degrees, "DipAngle_degrees/D");
    tree->Branch("OscParams", &fdr->OscParams, "OscParams[6]/D");
    tree->Branch("FromNuPDG", &fdr->FromNuPDG, "FromNuPDG/I");
    tree->Branch("ToNuPDG", &fdr->ToNuPDG, "ToNuPDG/I");
    fdr->Reset();
    return fdr;
  }

  void OscillationParameters::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  OscillationParameters::~OscillationParameters() {
    delete tree;
  }
