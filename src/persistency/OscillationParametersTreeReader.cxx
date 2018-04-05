#include "OscillationParametersTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

OscillationParameters::OscillationParameters(std::string const &inputFile) {
    LoadTree(inputFile);
  }

  std::string OscillationParameters::TreeName(){
    return "OscConfigTree";
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

  OscillationParameters *OscillationParameters::MakeTreeWriter() {
    OscillationParameters *fdr = new OscillationParameters();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->tree->Branch("DipAngle_degrees", &fdr->DipAngle_degrees, "DipAngle_degrees/D");
    fdr->tree->Branch("OscParams", &fdr->OscParams, "OscParams[6]/D");
    fdr->tree->Branch("FromNuPDG", &fdr->FromNuPDG, "FromNuPDG/I");
    fdr->tree->Branch("ToNuPDG", &fdr->ToNuPDG, "ToNuPDG/I");
    fdr->Reset();
    return fdr;
  }

  OscillationParameters::~OscillationParameters() {}
