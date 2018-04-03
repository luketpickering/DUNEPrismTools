#include "CovarianceHelper.h"
#include "EDepTreeReader.h"
#include "VALORModelClassifier.h"

#include "Utils.hxx"

#include "omp.h"

#include <string>
#include <vector>

std::string InputDepostisSummaryFile;
std::string OutputFile;

std::string CovmatFile, CovmatHist;
Long64_t NThrows = 1;
UInt_t Seed = 1;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
            "\t-i <stopprocessor.root>     : TChain descriptor for"
            " input tree. \n"
            "\t                              Should be the output of RunSelection."
            "\t-o <outputfile.root>        : Output file to write "
            "efficiency friend tree to.\n"
            "\t-C <covmatfile.root,hname>  : Input covariance matrix to "
            "throw xsec weights from.\n"
            "\t-N <NThrows>                : Number of throws to make. Special "
            "cases are "
            "\t                              1: will produce +1 sigma variation"
            "\t                              2: will produce +1 and -1 sigma "
            "variation.\n"
            "\t-S <Seed>                   : Used to make reproducible "
            "throws. {Default = 1}"
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
      InputDepostisSummaryFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-C") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -C, expected 2." << std::endl;
        throw;
      }
      CovmatFile = params[0];
      CovmatHist = params[1];
    } else if (std::string(argv[opt]) == "-N") {
      NThrows = str2T<Long64_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-S") {
      Seed = str2T<UInt_t>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  TMatrixDSym *UncertMatrix;
  if (CovmatFile.size() && CovmatHist.size()) {
  } else {
    UncertMatrix =
        new TMatrixDSym(static_cast<int>(VALORModel::TrueClass::kNVALORDials));
    for (int d_it = 0;
         d_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++d_it) {
      (*UncertMatrix)[d_it][d_it] =
          pow(GetDefaultDialTweak(static_cast<VALORModel::TrueClass>(d_it)), 2);
    }
    for (int d_it = 0;
         d_it < (static_cast<int>(VALORModel::TrueClass::kNVALORDials) - 1);
         ++d_it) {
      (*UncertMatrix)[d_it][d_it + 1] = (*UncertMatrix)[d_it][d_it] * 0.1;
      (*UncertMatrix)[d_it + 1][d_it] = (*UncertMatrix)[d_it][d_it] * 0.1;
    }
  }

  EDep edr("EDeps", InputDepostisSummaryFile);

  TFile *of = new TFile(OutputFile.c_str(), "RECREATE");
  TTree *friendtree = new TTree("XSecWeights", "");
  TTree *configtree = new TTree("ConfigTree", "");
  TTree *throwConfigtree = new TTree("ThrowConfigTree", "");

  configtree->Branch("NThrows", &NThrows, "NThrows/I");
  Int_t NParams = static_cast<int>(VALORModel::TrueClass::kNVALORDials);
  configtree->Branch("NParams", &NParams, "NParams/I");
  configtree->Fill();

  double *ParamValues = new double[NParams];
  throwConfigtree->Branch(
      "ParamValues", ParamValues,
      (std::string("ParamValues[") + to_str(NParams) + "]/D").c_str());

  std::vector<std::vector<double> > ParameterThrows;
  ParameterThrows.resize(NThrows);
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    ParameterThrows[t_it].resize(
        static_cast<int>(VALORModel::TrueClass::kNVALORDials));
  }
  CovarianceThrower ct(*UncertMatrix, Seed);

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    TMatrixD const *CovVector = ct.Throw();
    for (int p_it = 0;
         p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {
      ParameterThrows[t_it][p_it] = (*CovVector)[p_it][0];
      ParamValues[p_it] = ParameterThrows[t_it][p_it];
    }
    throwConfigtree->Fill();
  }

  for (int p_it = 0;
       p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {
    std::cout << "[INFO]: Param " << static_cast<VALORModel::TrueClass>(p_it)
              << ", Uncert "
              << GetDefaultDialTweak(static_cast<VALORModel::TrueClass>(p_it))
              << ", Throws: " << std::flush;
    for (int t_it = 0; t_it < std::min(NThrows, Long64_t(5)); ++t_it) {
      std::cout << ParameterThrows[t_it][p_it] << " " << std::flush;
    }
    std::cout << std::endl;
  }

  double *xsecweights = new double[NThrows];
  friendtree->Branch(
      "XSecWeights", xsecweights,
      (std::string("XSecWeights[") + to_str(NThrows) + "]/D").c_str());

  size_t loud_every = edr.GetEntries() / 10;
  size_t NFills = 0;
  Long64_t NEntries = edr.GetEntries();

  TH1D *MeanThrow = new TH1D("CentralERecDist", ";ERec;Count", 100, 0, 10);
  std::vector<TH1D *> ThrowDistributions;
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    ThrowDistributions.push_back(new TH1D("throw_dist", "", 100, 0, 10));
    ThrowDistributions.back()->SetDirectory(nullptr);
  }

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    if (edr.stop == -1) {
      xsecweights[0] = 1;
      NFills++;
      friendtree->Fill();
      continue;
    }

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    std::vector<VALORModel::TrueClass> appdial = GetApplicableDials(edr);

    if (NThrows < 3) {
#ifdef DEBUG
      std::string evc = edr.EventCode->GetString().Data();

      GENIECodeStringParser gcp(evc);
      std::cout << "[INFO]: event " << e_it << ", code: " << evc << std::endl;
      size_t d_it = 0;
#endif
      xsecweights[0] = 1;
      for (auto &d : appdial) {
#ifdef DEBUG
        std::cout << "\t [" << d_it << "]: " << d << ", +1 sigma weight: "
                  << GetVALORWeight(d, GetDefaultDialTweak(d), edr)
                  << std::endl;
        d_it++;
#endif
        xsecweights[0] *= GetVALORWeight(d, 1, edr);
      }

#ifdef DEBUG
      d_it = 0;
#endif
      if (NThrows == 2) {
        xsecweights[1] = 1;
        for (auto &d : appdial) {
#ifdef DEBUG
          std::cout << "\t [" << d_it << "]: " << d << ", +1 sigma weight: "
                    << GetVALORWeight(d, -GetDefaultDialTweak(d), edr)
                    << std::endl;
          d_it++;
#endif
          xsecweights[1] *= GetVALORWeight(d, -1, edr);
        }
      }
    } else {  // Do throws;

#pragma omp parallel for
      for (int t_it = 0; t_it < NThrows; ++t_it) {
        xsecweights[t_it] = 1;
        for (VALORModel::TrueClass d : appdial) {
          xsecweights[t_it] *= GetVALORWeight(
              d, ParameterThrows[t_it][static_cast<int>(d)], edr);
        }
        ThrowDistributions[t_it]->Fill(edr.PrimaryLep_4mom[3] +
                                           edr.TotalNonlep_Dep_FV +
                                           edr.TotalNonlep_Dep_veto,
                                       xsecweights[t_it]);
        MeanThrow->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                            edr.TotalNonlep_Dep_veto,
                        xsecweights[t_it]);
      }
    }

    NFills++;
    friendtree->Fill();
  }

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
  ThrowDistributions[MaxIndex]->SetDirectory(of);
  ThrowDistributions[MinIndex]->SetDirectory(of);

  of->Write();
  of->Close();
}
