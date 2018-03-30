#include "SliceConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

SliceConfig::SliceConfig() : tree(nullptr), NFiles(0), NEntries(0), CEnt(0) {}
SliceConfig::SliceConfig(std::string const &treeName, std::string const &inputFile) {

    tree = OpenTChainWithFileList(treeName, inputFile);

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    std::cout << "[SliceConfig]: Loaded TChain with " << NEntries 
      << " entries." << std::endl;
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

  std::pair< std::vector<double>, std::vector<double> >
    SliceConfig::BuildXRangeBinsCoeffs(
      std::vector< std::pair<double, double> > const &XRanges,
      double const *Coeffs){

    std::vector<double> XRangeBins;
    std::vector<double> CoeffVector;

    std::cout << "[INFO]: XRange bins: " << std::flush;
    XRangeBins.push_back(XRanges[0].first);
    CoeffVector.push_back(Coeffs[0]);
    XRangeBins.push_back(XRanges[0].second);
    std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
    for (size_t i = 1; i < XRanges.size(); ++i) {
      // If non-contiguous, must push an empty bit between.
      if (fabs(XRangeBins.back() - XRanges[i].first) > 1E-5) {
        CoeffVector.push_back(0);
        XRangeBins.push_back(XRanges[i].first);
        std::cout << ", " << XRangeBins.back() << std::flush;
      }
      CoeffVector.push_back(Coeffs[i]);
      XRangeBins.push_back(XRanges[i].second);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }
    std::cout << std::endl;
    return std::make_pair(XRangeBins, CoeffVector);
  }

  void SliceConfig::ReadTree(){
    if(!tree){
      std::cout << "[ERROR]: SliceConfig attempted to read input tree "
      "but it was not initialized. " << std::endl;
      throw;
    }

    XRanges.clear();
    Coeffs.clear();

    XRangeBins.push_back(XRange[1]);
    std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      XRanges.push_back(std::make_pair(XRange[0],XRange[1]));
      Coeffs.push_back(Coeff);
    }

    std::pair< std::vector<double>, std::vector<double> > XRBC =
      BuildXRangeBinsCoeffs(XRanges, Coeffs.data());
    XRangeBins = XRBC.first;
    BinCoeffsVector = XRBC.second;
  }

  std::vector< std::pair<double, double> > SliceConfig::GetXRanges() {
    if(!XRanges.size()){
      ReadTree();
    }
    return XRanges;
  }

  std::vector<double> SliceConfig::GetCoeffs() {
    if(!Coeffs.size()){
      ReadTree();
    }
    return Coeffs;
  }

  std::vector<double> SliceConfig::GetXRangeBins() {
    if(!XRangeBins.size()){
      ReadTree();
    }
    return XRangeBins;
  }

  std::vector<double> SliceConfig::GetBinCoeffs() {
    if(!BinCoeffsVector.size()){
      ReadTree();
    }
    return BinCoeffsVector;
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
    fdr->Reset();
    return fdr;
  }

  void SliceConfig::ReleaseInputFile(){
    delete tree;
    tree = nullptr;
  }

  SliceConfig::~SliceConfig() {
    delete tree;
  }
