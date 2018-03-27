#include "SliceConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

SliceConfig::SliceConfig() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
SliceConfig::SliceConfig(std::string const &treeName, std::string const &inputFile) {

    tree = OpenTChainWithFileList(treeName, inputFile);

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[SliceConfig]: Loaded TChain: " << NFiles
              << " files and " << NEntries << " entries." << std::endl;
    GetEntry(0);
  }

  void SliceConfig::Reset() {
    std::fill_n(XRange,2,0);
    Coeff = 0;
  }

  void SliceConfig::SetBranchAddresses() {
    tree->SetBranchAddress("XRange", &XRange);
    tree->SetBranchAddress("Coeff", &Coeff);
  }

  void SliceConfig::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t SliceConfig::GetEntry() { return CEnt; }
  UInt_t SliceConfig::GetEntries() { return NEntries; }

  void SliceConfig::ReadTree(){
    if(!tree){
      std::cout << "[ERROR]: SliceConfig attempted to read input tree "
      "but it was not initialized. " << std::endl;
      throw;
    }

    XRangeBins.clear();
    CoeffsVector.clear();

    std::cout << "[INFO]: XRange bins: " << std::flush;
    tree->GetEntry(0);
    XRangeBins.push_back(XRange[0]);
    CoeffsVector.push_back(Coeff);
    XRangeBins.push_back(XRange[1]);
    std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
    for (Long64_t i = 1; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      // If non-contiguous, must push an empty bit between.
      if (fabs(XRangeBins.back() - XRange[0]) > 1E-5) {
        CoeffsVector.push_back(0);
        XRangeBins.push_back(XRange[0]);
        std::cout << ", " << XRangeBins.back() << std::flush;
      }

      CoeffsVector.push_back(Coeff);
      XRangeBins.push_back(XRange[1]);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }
    std::cout << std::endl;

  }

  std::vector<double> SliceConfig::GetXRangeBins() {
    if(!XRangeBins.size()){
      ReadTree();
    }
    return XRangeBins;
  }

  std::vector<double> SliceConfig::GetCoeffs() {
    if(!CoeffsVector.size()){
      ReadTree();
    }
    return CoeffsVector;
  }

  TH1D *SliceConfig::BuildSliceBinningHelper(std::string const &histName){

    TH1D * BinningHelper = new TH1D(histName.c_str(), "",
                                    (GetXRangeBins().size() - 1), GetXRangeBins().data());
    BinningHelper->SetDirectory(nullptr);

    for (size_t bin_it = 1; bin_it < GetXRangeBins().size(); ++bin_it) {
      BinningHelper->SetBinContent(bin_it, GetCoeffs()[bin_it - 1]);
    }
    return BinningHelper;
  }


  SliceConfig *SliceConfig::MakeTreeWriter(TTree *tree) {
    SliceConfig *fdr = new SliceConfig();
    tree->Branch("XRange", &fdr->XRange, "XRange[2]/D");
    tree->Branch("Coeff", &fdr->Coeff, "Coeff/D");
    return fdr;
  }

  void SliceConfig::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  SliceConfig::~SliceConfig() {
    delete tree;
  }
