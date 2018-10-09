#include "SliceConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

SliceConfig::SliceConfig(std::string const &inputFile,
                         std::string const &inputDir) {
  dir = inputDir;
  if (dir.size() && (dir.back() != '/')) {
    dir += "/";
  }

  LoadTree(inputFile);
}

std::string SliceConfig::TreeName() { return dir + "SliceConfigTree"; }

void SliceConfig::Reset() {
  std::fill_n(XRange, 2, 0);
  Coeff = 0;
}

void SliceConfig::SetBranchAddresses() {
  tree->SetBranchAddress("XRange", &XRange);
  tree->SetBranchAddress("Coeff", &Coeff);
}

std::pair<std::vector<double>, std::vector<double>>
SliceConfig::BuildXRangeBinsCoeffs(
    std::vector<std::pair<double, double>> const &XRanges, double const *Coeffs,
    bool SignFlipX) {

  std::vector<double> XRangeBins;
  std::vector<double> CoeffVector;

  std::cout << "[INFO]: XRange bins: " << std::flush;
  XRangeBins.push_back((SignFlipX ? -1 : 1) * XRanges[0].first);
  CoeffVector.push_back(Coeffs[0]);
  XRangeBins.push_back((SignFlipX ? -1 : 1) * XRanges[0].second);
  std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
  for (size_t i = 1; i < XRanges.size(); ++i) {
    // If non-contiguous, must push an empty bit between.
    if (fabs(XRangeBins.back() - ((SignFlipX ? -1 : 1) * XRanges[i].first)) >
        1E-5) {
      CoeffVector.push_back(0);
      XRangeBins.push_back((SignFlipX ? -1 : 1) * XRanges[i].first);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }
    CoeffVector.push_back(Coeffs[i]);
    XRangeBins.push_back((SignFlipX ? -1 : 1) * XRanges[i].second);
    std::cout << ", " << XRangeBins.back() << std::flush;
  }
  std::cout << std::endl;

  if (SignFlipX) {
    std::vector<double> XRangeBins_rev;
    XRangeBins_rev.resize(XRangeBins.size());
    std::vector<double> CoeffVector_rev;
    CoeffVector_rev.resize(CoeffVector.size());

    size_t ctr = 0;
    for (std::vector<double>::reverse_iterator xrb_it = XRangeBins.rbegin();
         xrb_it < XRangeBins.rend(); ++xrb_it) {
      XRangeBins_rev[ctr++] = (*xrb_it);
    }

    ctr = 0;
    for (std::vector<double>::reverse_iterator cvr_it = CoeffVector.rbegin();
         cvr_it < CoeffVector.rend(); ++cvr_it) {
      CoeffVector_rev[ctr++] = (*cvr_it);
    }

    XRangeBins.swap(XRangeBins_rev);
    CoeffVector.swap(CoeffVector_rev);
  }

  return std::make_pair(XRangeBins, CoeffVector);
}

void SliceConfig::ReadTree() {
  if (!tree) {
    std::cout << "[ERROR]: SliceConfig attempted to read input tree "
                 "but it was not initialized. "
              << std::endl;
    throw;
  }

  XRanges.clear();
  Coeffs.clear();

  XRangeBins.push_back(XRange[1]);
  std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    XRanges.push_back(std::make_pair(XRange[0], XRange[1]));
    Coeffs.push_back(Coeff);
  }

  std::pair<std::vector<double>, std::vector<double>> XRBC =
      BuildXRangeBinsCoeffs(XRanges, Coeffs.data());
  XRangeBins = XRBC.first;
  BinCoeffsVector = XRBC.second;
}

std::vector<std::pair<double, double>> SliceConfig::GetXRanges() {
  if (!XRanges.size()) {
    ReadTree();
  }
  return XRanges;
}

std::vector<double> SliceConfig::GetCoeffs() {
  if (!Coeffs.size()) {
    ReadTree();
  }
  return Coeffs;
}

std::vector<double> SliceConfig::GetXRangeBins() {
  if (!XRangeBins.size()) {
    ReadTree();
  }
  return XRangeBins;
}

std::vector<double> SliceConfig::GetBinCoeffs() {
  if (!BinCoeffsVector.size()) {
    ReadTree();
  }
  return BinCoeffsVector;
}

TH1D *SliceConfig::BuildSliceBinningHelper(std::string const &histName) {

  TH1D *BinningHelper =
      new TH1D(histName.c_str(), "", (GetXRangeBins().size() - 1),
               GetXRangeBins().data());
  BinningHelper->SetDirectory(nullptr);

  for (size_t bin_it = 1; bin_it < GetXRangeBins().size(); ++bin_it) {
    BinningHelper->SetBinContent(bin_it, GetCoeffs()[bin_it - 1]);
  }
  return BinningHelper;
}

SliceConfig *SliceConfig::MakeTreeWriter() {
  SliceConfig *fdr = new SliceConfig();
  fdr->tree = new TTree(fdr->TreeName().c_str(), "");
  fdr->TreeOwned = false;

  fdr->tree->Branch("XRange", &fdr->XRange, "XRange[2]/D");
  fdr->tree->Branch("Coeff", &fdr->Coeff, "Coeff/D");
  fdr->Reset();
  return fdr;
}

SliceConfig::~SliceConfig() {}
