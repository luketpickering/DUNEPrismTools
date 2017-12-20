#include "Utils.hxx"

#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"

#include <iostream>
#include <map>
#include <utility>

/// -n Nominal.root -v Name,Up.root,Down.root -E 3 -X 35 -M 1 -o Output.root

struct FlatHistTMatrixD {
  TMatrixD* BinContent;
  TMatrixD* BinErrors;
  bool IsTranspose;

  FlatHistTMatrixD() {
    BinContent = nullptr;
    BinErrors = nullptr;
    IsTranspose = false;
  }
  FlatHistTMatrixD(size_t NRows, bool Transpose = false) {
    BinContent =
        Transpose ? new TMatrixD(0, NRows, 0, 1) : new TMatrixD(0, 1, 0, NRows);
    BinErrors =
        Transpose ? new TMatrixD(0, NRows, 0, 1) : new TMatrixD(0, 1, 0, NRows);
    IsTranspose = Transpose;
  }

  FlatHistTMatrixD& operator=(FlatHistTMatrixD const & other) {
    IsTranspose = other.IsTranspose;

    BinContent = IsTranspose ? new TMatrixD(0, other.Size(), 0, 1)
                             : new TMatrixD(0, 1, 0, other.Size());
    BinErrors = IsTranspose ? new TMatrixD(0, other.Size(), 0, 1)
                            : new TMatrixD(0, 1, 0, other.Size());

    for (size_t i = 0; i < Size(); ++i) {
      Set(i, 0);
      ESet(i, 0);
    }

    return *this;
  }

  FlatHistTMatrixD(FlatHistTMatrixD const & other) {
    IsTranspose = other.IsTranspose;

    BinContent = IsTranspose ? new TMatrixD(0, other.Size(), 0, 1)
                             : new TMatrixD(0, 1, 0, other.Size());
    BinErrors = IsTranspose ? new TMatrixD(0, other.Size(), 0, 1)
                            : new TMatrixD(0, 1, 0, other.Size());

    for (size_t i = 0; i < Size(); ++i) {
      Set(i, 0);
      ESet(i, 0);
    }

  }

  ~FlatHistTMatrixD() {
    delete BinContent;
    delete BinErrors;
  }

  size_t Size() const {
    return IsTranspose ? BinContent->GetNcols() : BinContent->GetNrows();
  }

  double At(size_t i) const {
    return IsTranspose ? (*BinContent)[0][i] : (*BinContent)[i][0];
  }
  double Set(size_t i, double val) {
    return (IsTranspose ? (*BinContent)[0][i] : (*BinContent)[i][0]) = val;
  }

  double EAt(size_t i) const {
    return IsTranspose ? (*BinErrors)[0][i] : (*BinErrors)[i][0];
  }
  double ESet(size_t i, double val) {
    return (IsTranspose ? (*BinErrors)[0][i] : (*BinErrors)[i][0]) = val;
  }

  FlatHistTMatrixD EmptyClone() {
    FlatHistTMatrixD cl(Size(), IsTranspose);

    for (size_t i = 0; i < Size(); ++i) {
      cl.Set(i, 0);
      cl.ESet(i, 0);
    }
    return cl;
  }

  FlatHistTMatrixD GetTranspose() {
    FlatHistTMatrixD cl(Size(), !IsTranspose);

    for (size_t i = 0; i < Size(); ++i) {
      cl.Set(i, At(i));
      cl.ESet(i, EAt(i));
    }
    return cl;
  }
};

FlatHistTMatrixD Nominal;
std::map<std::string, std::pair<FlatHistTMatrixD, FlatHistTMatrixD> >
    VariedFluxes;
std::map<std::string, FlatHistTMatrixD> FractionalVariations;
std::map<std::string, TMatrixD*> IndividualMatrices;
TMatrixD* TotalMatrix;

std::string NominalInput;
std::map<std::string, std::pair<std::string, std::string> > VariationalInputs;

std::string OutputFile;

double MaxEnergy = std::numeric_limits<double>::max();
double MaxAbsoluteOffset = std::numeric_limits<double>::max();
int BinMergeScheme = 0;

size_t NNominalBins;

void SayUsage(char const* argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
}

void handleOpts(int argc, char const* argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-n") {
      NominalInput = argv[++opt];
      std::cout << "--Reading nominal flux from " << NominalInput << std::endl;
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
      std::cout << "--Reading nominal flux from " << NominalInput << std::endl;
    } else if (std::string(argv[opt]) == "-v") {
      std::vector<std::string> VariationDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");

      if (VariationDescriptors.size() != 3) {
        std::cout << "[ERROR]: Expected to find an option like: -v "
                     "<VariationName>,<OneSigUp.root>,<OneSigDown.root>, but "
                     "found \""
                  << argv[opt] << "\"." << std::endl;
        exit(1);
      }
      VariationalInputs[VariationDescriptors[0]] =
          std::make_pair(VariationDescriptors[1], VariationDescriptors[2]);

      std::cout << "--Reading " << VariationDescriptors[0]
                << " varied fluxes from " << VariationDescriptors[1] << " and "
                << VariationDescriptors[2] << std::endl;

    } else if (std::string(argv[opt]) == "-E") {
      MaxEnergy = str2T<double>(argv[++opt]);
      std::cout << "--Maximum energy included in histograms: " << MaxEnergy
                << " GeV." << std::endl;
    } else if (std::string(argv[opt]) == "-X") {
      MaxAbsoluteOffset = str2T<double>(argv[++opt]);
      std::cout << "--Maximum offset included in histograms: "
                << MaxAbsoluteOffset << " m." << std::endl;
    } else if (std::string(argv[opt]) == "-M") {
      BinMergeScheme = str2T<int>(argv[++opt]);
      std::cout << "--Merging " << (BinMergeScheme + 1) << " bins."
                << std::endl;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

FlatHistTMatrixD BuildFluxVector(std::string const& InputFileName) {
  TH2D* fluxhist = GetHistogram<TH2D>(InputFileName, "numu_flux_2D");

  if (BinMergeScheme) {
    fluxhist->RebinX(BinMergeScheme + 1);
  }

  std::vector<double> content;
  std::vector<double> errors;

  Int_t dummy = 0;

  for (Int_t x_it = 1; x_it < fluxhist->GetXaxis()->GetNbins() + 1; ++x_it) {
    double ble = fluxhist->GetXaxis()->GetBinLowEdge(x_it);
    if ((fabs(ble - MaxEnergy) < 1E-5) || (ble > MaxEnergy)) {
      break;
    }

    for (Int_t y_it = 1; y_it < fluxhist->GetYaxis()->GetNbins() + 1; ++y_it) {
      double ble = fluxhist->GetYaxis()->GetBinLowEdge(y_it);
      if ((fabs(ble - MaxAbsoluteOffset) < 1E-5) || (ble > MaxAbsoluteOffset)) {
        break;
      }

      Int_t GBin = fluxhist->GetBin(x_it, y_it, dummy);
      content.push_back(fluxhist->GetBinContent(GBin));
      errors.push_back(fluxhist->GetBinError(GBin));
    }
  }

  FlatHistTMatrixD fh(content.size());

  for (size_t i = 0; i < content.size(); ++i) {
    fh.Set(i, content[i]);
    fh.ESet(i, content[i]);
  }

  return fh;
}

int main(int argc, char const* argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  // Build Nominal
  Nominal = BuildFluxVector(NominalInput);
  NNominalBins = Nominal.Size();

  // Build Variational
  for (std::map<std::string, std::pair<std::string, std::string> >::iterator
           var_it = VariationalInputs.begin();
       var_it != VariationalInputs.end(); ++var_it) {
    VariedFluxes[var_it->first] =
        std::make_pair(BuildFluxVector(var_it->second.first),
                       BuildFluxVector(var_it->second.second));

    if (VariedFluxes[var_it->first].first.Size() != NNominalBins) {
      std::cout << "[ERROR]: Flux histogram found in " << var_it->second.first
                << " had " << VariedFluxes[var_it->first].first.Size()
                << " bins after processing, but the nominal flux histogram has "
                << NNominalBins << std::endl;
      exit(2);
    }
    if (VariedFluxes[var_it->first].second.Size() != NNominalBins) {
      std::cout << "[ERROR]: Flux histogram found in " << var_it->second.second
                << " had " << VariedFluxes[var_it->first].second.Size()
                << " bins after processing, but the nominal flux histogram has "
                << NNominalBins << std::endl;
      exit(3);
    }
  }

  // Build Averaged Fractions
  for (std::map<std::string,
                std::pair<FlatHistTMatrixD, FlatHistTMatrixD> >::iterator
           var_it = VariedFluxes.begin();
       var_it != VariedFluxes.end(); ++var_it) {
    FractionalVariations[var_it->first] = Nominal.EmptyClone();

    for (size_t bi_it = 0; bi_it < NNominalBins; ++bi_it) {
      double UpFracVar =
          fabs(var_it->second.first.At(bi_it) - Nominal.At(bi_it)) /
          Nominal.At(bi_it);
      double DownFracVar =
          fabs(var_it->second.second.At(bi_it) - Nominal.At(bi_it)) /
          Nominal.At(bi_it);

      double UpUncert = var_it->second.first.EAt(bi_it);
      double DownUncert = var_it->second.second.EAt(bi_it);

      FractionalVariations[var_it->first].Set(
          bi_it, (UpFracVar * (1.0 / (UpUncert * UpUncert)) +
                  DownFracVar * (1.0 / (DownUncert * DownUncert))) /
                     ((1.0 / (UpUncert * UpUncert)) +
                      (1.0 / (DownUncert * DownUncert))));
    }
  }

  // Build Indivudial Matrices
  for (std::map<std::string, FlatHistTMatrixD>::iterator var_it =
           FractionalVariations.begin();
       var_it != FractionalVariations.end(); ++var_it) {
    FlatHistTMatrixD transp = var_it->second.GetTranspose();
    (*transp.BinContent) * (*var_it->second.BinContent);

    if (transp.BinContent->GetNrows() != int(NNominalBins)) {
      std::cout << "[ERROR]: Expected post multiplcation matrix to be square. "
                   "NRows = "
                << transp.BinContent->GetNrows() << "!= " << NNominalBins
                << std::endl;
      exit(4);
    }
    if (transp.BinContent->GetNcols() != int(NNominalBins)) {
      std::cout << "[ERROR]: Expected post multiplcation matrix to be square. "
                   "NCols = "
                << transp.BinContent->GetNrows() << "!= " << NNominalBins
                << std::endl;
      exit(5);
    }

    IndividualMatrices[var_it->first] = transp.BinContent;
    transp.BinContent = nullptr;
  }

  // Build total matrix
  TotalMatrix = new TMatrixD(0, NNominalBins, 0, NNominalBins);
  for (std::map<std::string, TMatrixD*>::iterator var_it =
           IndividualMatrices.begin();
       var_it != IndividualMatrices.end(); ++var_it) {
    for (size_t x_it = 0; x_it < NNominalBins; ++x_it) {
      for (size_t y_it = 0; y_it < NNominalBins; ++y_it) {
        (*TotalMatrix)[y_it][x_it] += (*var_it->second)[y_it][x_it];
      }
    }
  }
  // Write out

  TFile* of = new TFile(OutputFile.c_str(), "RECREATE");

  for (std::map<std::string, TMatrixD*>::iterator var_it =
           IndividualMatrices.begin();
       var_it != IndividualMatrices.end(); ++var_it) {
    var_it->second->Write((var_it->first + "_Uncertainty").c_str(),
                          TObject::kOverwrite);
  }
  TotalMatrix->Write("Total_Uncertainty", TObject::kOverwrite);

  of->Write();
  of->Close();
  delete of;

  for (std::map<std::string, TMatrixD*>::iterator var_it =
           IndividualMatrices.begin();
       var_it != IndividualMatrices.end(); ++var_it) {
    delete var_it->second;
  }

  delete TotalMatrix;
}
