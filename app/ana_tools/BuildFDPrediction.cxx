#include "CovarianceHelper.h"
#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "BargerPropagator.h"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::string inpfile;
std::string fluxthrowfile;
std::string xsecthrowfile;
std::string efffile;
std::string outputfile = "FarDetectorPrediction.root";

double DipAngle = 5.8;
double OscParams[6] = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
static const double deg2rad = asin(1) / 90.0;

std::vector<double> ERecBinning;

double HadrVeto = 10E-3;  // GeV

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>             : TChain descriptor for"
         " input tree. \n"
         "\t-uF <fluxthrowfriend.root>           : Friend tree for -i "
         "option containing flux throws.\n"
         "\t-uX <xsecthrowfriend.root>           : Friend tree for -i "
         "option containing xsec throws.\n"
         "\t-E <efficiencyfriend.root>           : Friend tree for -i "
         "option containing efficiency correction weights.\n"
         "\t-p "
         "<sin2(theta12)=0.825>,<sin2(theta13)=0.10>,<sin2(theta23)=1.0>,"
         "<dm12=7.9e-5>,<dm23=2.5e-3>,<dcp=0.0>\n"
         "\t-d <dipangle=5.8>\n"
         "\t-o <outputfile.root>                 : Output file name.\n"
         "\t-b <low_up:width[,low_up:width...]>  : ENuBinning descriptor.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
      << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-uF") {
      fluxthrowfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-uX") {
      xsecthrowfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-E") {
      efffile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-d") {
      DipAngle = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-p") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 6) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -p, expected 6." << std::endl;
        exit(1);
      }

      OscParams[0] = params[0];
      OscParams[1] = params[1];
      OscParams[2] = params[2];
      OscParams[3] = params[3];
      OscParams[4] = params[4];
      OscParams[5] = params[5];

      std::cout << "Sin^2(Theta_12) = " << OscParams[0] << std::endl;
      std::cout << "Sin^2(Theta_13) = " << OscParams[1] << std::endl;
      std::cout << "Sin^2(Theta_23) = " << OscParams[2] << std::endl;

      std::cout << "Dm^2_21 = " << OscParams[3] << std::endl;
      std::cout << "|Dm^2_Atm| = " << OscParams[4] << std::endl;

      std::cout << "dcp = " << OscParams[5] << std::endl;

    } else if (std::string(argv[opt]) == "-b") {
      std::vector<std::string> binDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");
      ERecBinning.clear();
      for (size_t vbd_it = 0; vbd_it < binDescriptors.size(); ++vbd_it) {
        AppendVect(ERecBinning, BuildDoubleList(binDescriptors[vbd_it]));
      }

      for (size_t bin_it = 1; bin_it < ERecBinning.size(); ++bin_it) {
        if (ERecBinning[bin_it] == ERecBinning[bin_it - 1]) {
          std::cout << "[INFO]: Removing duplciate bin low edge " << bin_it
                    << " low edge: " << ERecBinning[bin_it] << std::endl;
          ERecBinning.erase(ERecBinning.begin() + bin_it);
        }
      }

      for (size_t bin_it = 1; bin_it < ERecBinning.size(); ++bin_it) {
        if (ERecBinning[bin_it] < ERecBinning[bin_it - 1]) {
          std::cout << "[ERROR]: Bin " << bin_it
                    << " low edge: " << ERecBinning[bin_it]
                    << " is smaller than bin " << (bin_it - 1)
                    << " low edge: " << ERecBinning[bin_it - 1] << std::endl;
          exit(1);
        }
      }
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

enum nuTypes {
  kNuebarType = -1,
  kNumubarType = -2,
  kNutaubarType = -3,
  kNueType = 1,
  kNumuType = 2,
  kNutauType = 3,
};

nuTypes GetNuType(int pdg) {
  switch (pdg) {
    case 16:
      return kNutauType;
    case 14:
      return kNumuType;
    case 12:
      return kNueType;
    case -16:
      return kNutaubarType;
    case -14:
      return kNumubarType;
    case -12:
      return kNuebarType;
    default: {
      std::cout << "[ERROR]: Attempting to convert \"neutrino pdg\": " << pdg
                << std::endl;
      exit(1);
    }
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!inpfile.size()) {
    std::cout << "[ERROR]: Was not passed both an input analysis file and a "
                 "flux fit file."
              << std::endl;
  }

  EDep edr("EDeps", inpfile);
  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  BargerPropagator bp;
  double lengthParam = cos((90.0 + DipAngle) * deg2rad);

  if (!ERecBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.25"};
    handleOpts(argc_dum, argv_dum);
  }

  TFile *outfile = CheckOpenFile(outputfile, "RECREATE");

  TH1D *FDPrediction = new TH1D("FDPrediction_ERec", ";E_{Dep} (GeV);Count",
                                (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDPrediction_True =
      new TH1D("FDPrediction_True", ";E_{#nu} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());

  TH1D *FDPredictionOsc =
      new TH1D("FDPredictionOsc_ERec", ";E_{Dep} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDPredictionOsc_True =
      new TH1D("FDPredictionOsc_True", ";E_{#nu} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());

  TH2D *ERecFD = new TH2D("ERecFD", ";ETrue;ERec;Count",
                          (ERecBinning.size() - 1), ERecBinning.data(),
                          (ERecBinning.size() - 1), ERecBinning.data());
  TH2D *ERecFDosc = new TH2D("ERecFDOsc", ";ETrue;ERec;Count",
                             (ERecBinning.size() - 1), ERecBinning.data(),
                             (ERecBinning.size() - 1), ERecBinning.data());
  TH2D *EProxyFD = new TH2D("EProxyFD", ";ETrue;EProxy;Count",
                            (ERecBinning.size() - 1), ERecBinning.data(),
                            (ERecBinning.size() - 1), ERecBinning.data());
  TH2D *EProxyFDosc = new TH2D("EProxyFDOsc", ";ETrue;EProxy;Count",
                               (ERecBinning.size() - 1), ERecBinning.data(),
                               (ERecBinning.size() - 1), ERecBinning.data());

  nuTypes numu = GetNuType(14);

  // Set up EffTree
  TChain *EffFriendTree = nullptr;
  TH1D *FDPrediction_True_EffWeighted = nullptr;
  TH1D *FDPredictionOsc_True_EffWeighted = nullptr;
  double EffWeight = 1;

  if (efffile.size()) {
    EffFriendTree = new TChain("EffWeights");

    EffFriendTree->Add(efffile.c_str());

    EffFriendTree->SetBranchAddress("EffWeight", &EffWeight);

    FDPrediction_True_EffWeighted =
        new TH1D("FDPrediction_True_EffWeighted", "", (ERecBinning.size() - 1),
                 ERecBinning.data());
    FDPredictionOsc_True_EffWeighted =
        new TH1D("FDPredictionOsc_True_EffWeighted", "",
                 (ERecBinning.size() - 1), ERecBinning.data());
  }

  TChain *XSecThrowFriendTree = nullptr;
  TChain *XSecThrowFriendTree_config = nullptr;
  Int_t NThrows = 0;
  TH1D *MeanThrow = nullptr;
  std::vector<TH1D *> ThrowDistributions;
  double *xsecweights = nullptr;
  if (xsecthrowfile.size()) {
    XSecThrowFriendTree_config = new TChain("ConfigTree");
    XSecThrowFriendTree_config->Add(xsecthrowfile.c_str());
    XSecThrowFriendTree_config->SetBranchAddress("NThrows", &NThrows);
    XSecThrowFriendTree_config->GetEntry(0);

    std::cout << "[INFO]: Found XSec throw tree with " << NThrows
              << " parameter throws." << std::endl;

    XSecThrowFriendTree = new TChain("XSecWeights");
    XSecThrowFriendTree->Add(xsecthrowfile.c_str());

    xsecweights = new double[NThrows];
    XSecThrowFriendTree->SetBranchAddress("XSecWeights", xsecweights);

    MeanThrow = new TH1D("CentralERecDist", ";ERec;Count",
                         (ERecBinning.size() - 1), ERecBinning.data());
    for (int t_it = 0; t_it < NThrows; ++t_it) {
      ThrowDistributions.push_back(new TH1D(
          "throw_dist", "", (ERecBinning.size() - 1), ERecBinning.data()));
      ThrowDistributions.back()->SetDirectory(nullptr);
    }
  }

  Long64_t NSel = 0;
  // Iterate through entries
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    if ((edr.stop < 0) || (edr.TotalNonlep_Dep_veto > HadrVeto) ||
        (edr.PrimaryLepPDG != 13)) {
      continue;
    }
    NSel++;

    bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
              OscParams[4], OscParams[5], edr.nu_4mom[3], true, numu);
    bp.DefinePath(lengthParam, 0);
    bp.propagate(numu);

    double oweight = bp.GetProb(numu, numu);

    FDPrediction->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                       edr.TotalNonlep_Dep_veto);
    FDPrediction_True->Fill(edr.nu_4mom[3]);
    FDPredictionOsc->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                              edr.TotalNonlep_Dep_veto,
                          oweight);
    FDPredictionOsc_True->Fill(edr.nu_4mom[3], oweight);

    ERecFD->Fill(edr.nu_4mom[3], edr.PrimaryLep_4mom[3] +
                                     edr.TotalNonlep_Dep_FV +
                                     edr.TotalNonlep_Dep_veto);
    ERecFDosc->Fill(edr.nu_4mom[3],
                    edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                        edr.TotalNonlep_Dep_veto,
                    oweight);

    EProxyFD->Fill(edr.nu_4mom[3], edr.ERecProxy_True);
    EProxyFDosc->Fill(edr.nu_4mom[3], edr.ERecProxy_True, oweight);

    if (EffFriendTree) {
      EffFriendTree->GetEntry(e_it);
      FDPrediction_True_EffWeighted->Fill(edr.nu_4mom[3], EffWeight);
      FDPredictionOsc_True_EffWeighted->Fill(edr.nu_4mom[3],
                                             oweight * EffWeight);
    }

    if (NThrows) {
      XSecThrowFriendTree->GetEntry(e_it);
      for (Int_t t_it = 0; t_it < NThrows; ++t_it) {
        ThrowDistributions[t_it]->Fill(edr.PrimaryLep_4mom[3] +
                                           edr.TotalNonlep_Dep_FV +
                                           edr.TotalNonlep_Dep_veto,
                                       oweight * EffWeight);
        MeanThrow->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                            edr.TotalNonlep_Dep_veto,
                        oweight * EffWeight);
      }
    }
  }

  if (NThrows) {
    CovarianceBuilder cb(MeanThrow->GetXaxis()->GetNbins());

    MeanThrow->Scale(1.0 / double(NThrows));

    cb.SetMean(MeanThrow);

    for (int t_it = 0; t_it < NThrows; ++t_it) {
      cb.AddThrow_CovMatCalc(ThrowDistributions[t_it]);
    }
    cb.FinalizeCovMatCalc();

    cb.GetCovMatrix()->Clone()->Write("ERec_CovMat");
    cb.GetCorrMatrix()->Clone()->Write("ERec_CorrMat");

    for (Int_t i = 0; i < MeanThrow->GetXaxis()->GetNbins(); ++i) {
      MeanThrow->SetBinError(i + 1, sqrt((*cb.GetCovMatrix())[i][i]));
    }

    double MaxThrow = -std::numeric_limits<double>::max();
    int MaxIndex = 0;
    double MinThrow = -std::numeric_limits<double>::max();
    int MinIndex = 0;

    for (int t_it = 0; t_it < NThrows; ++t_it) {
      if (ThrowDistributions[t_it]->GetMaximum() > MaxThrow) {
        MaxIndex = t_it;
        MaxThrow = ThrowDistributions[t_it]->GetMaximum();
      }
      if (ThrowDistributions[t_it]->GetMaximum() < MinThrow) {
        MinIndex = t_it;
        MinThrow = ThrowDistributions[t_it]->GetMaximum();
      }
    }

    ThrowDistributions[MaxIndex]->SetName(
        (std::string("MaxERecDist_throw") + to_str(MaxIndex)).c_str());
    ThrowDistributions[MinIndex]->SetName(
        (std::string("MinERecDist_throw") + to_str(MinIndex)).c_str());
    ThrowDistributions[MaxIndex]->SetDirectory(outfile);
    ThrowDistributions[MinIndex]->SetDirectory(outfile);
  }

  std::cout << "[INFO]: Used " << NSel << "/" << NEntries
            << " entries to build far detector prediction." << std::endl;

  outfile->Write();
  outfile->Close();
}
