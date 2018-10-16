#include "DepositsSummaryTreeReader.hxx"

#include "CovarianceHelper.hxx"
#include "VALORModelClassifier.hxx"

#include "GetUsage.hxx"
#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include <string>
#include <vector>
#include <memory>

std::string InputDepostisSummaryFile;
std::string OutputFile;

std::string CovmatFile, CovmatHist;
UInt_t NThrows = 1;
UInt_t NEvents = std::numeric_limits<UInt_t>::max();
UInt_t Seed = 1;
double DiagonalErrorInflataion = 0.0001;
double InversionTolerance = 1E-31;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "ana_tools") << std::endl;
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
      NThrows = str2T<UInt_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-D") {
      DiagonalErrorInflataion = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-T") {
      InversionTolerance = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-n") {
      NEvents = str2T<UInt_t>(argv[++opt]);
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
    std::unique_ptr<TH2> covhist(GetHistogram<TH2>(CovmatFile, CovmatHist));
    UncertMatrix =
        new TMatrixDSym(static_cast<int>(VALORModel::TrueClass::kNVALORDials));
    for (int x_it = 0;
         x_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++x_it) {
      for (int y_it = x_it;
           y_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials);
           ++y_it) {

        (*UncertMatrix)[x_it][y_it] =
            covhist->GetBinContent(x_it + 1, y_it + 1);

        if (y_it == x_it) {

          std::cout << "[INFO]: Reading covariance " << x_it << ", " << y_it
                    << " = " << covhist->GetBinContent(x_it + 1, y_it + 1)
                    << std::endl;
          (*UncertMatrix)[x_it][y_it] += DiagonalErrorInflataion;
        }
      }
    }
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

  DepositsSummary edr(InputDepostisSummaryFile);

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

  std::vector<std::vector<double>> ParameterThrows;
  ParameterThrows.resize(NThrows);
  for (UInt_t t_it = 0; t_it < NThrows; ++t_it) {
    ParameterThrows[t_it].resize(
        static_cast<int>(VALORModel::TrueClass::kNVALORDials));
  }
  CovarianceDecompThrower ct(*UncertMatrix, Seed, InversionTolerance);
  CovarianceBuilder paramThrowsCB(
      static_cast<int>(VALORModel::TrueClass::kNVALORDials));
  paramThrowsCB.SetZeroMean();

  TRandom3 *RNJesus = new TRandom3(Seed);
  for (Int_t t_it = 0; t_it < Int_t(NThrows); ++t_it) {
    TMatrixD const *CovVector = ct.Throw();
    for (int p_it = 0;
         p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {

      // input matrix has no MEC uncertainty, throw flat weight 0--2
      if ((static_cast<VALORModel::TrueClass>(p_it) ==
           VALORModel::TrueClass::kNu_MEC) ||
          (static_cast<VALORModel::TrueClass>(p_it) ==
           VALORModel::TrueClass::kNuBar_MEC)) {
        ParameterThrows[t_it][p_it] = RNJesus->Uniform(-1, 1);
        ParamValues[p_it] = ParameterThrows[t_it][p_it];
      } else {
        ParameterThrows[t_it][p_it] = (*CovVector)[p_it][0];
        ParamValues[p_it] = ParameterThrows[t_it][p_it];
      }
      //Lets not allow any negative values
      if (ParameterThrows[t_it][p_it] < -1) {
        t_it = t_it - 1;
        continue;
      }
    }
    paramThrowsCB.AddThrow_CovMatCalc(CovVector);
    throwConfigtree->Fill();
  }

  for (int p_it = 0;
       p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {
    std::cout << "[INFO]: Param " << static_cast<VALORModel::TrueClass>(p_it)
              << ", Uncert "
              << GetDefaultDialTweak(static_cast<VALORModel::TrueClass>(p_it))
              << ", Throws: " << std::flush;
    for (UInt_t t_it = 0; t_it < std::min(NThrows, UInt_t(5)); ++t_it) {
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
  UInt_t NEntries = std::min(edr.GetEntries(), NEvents);

  TH1D *MeanThrow = new TH1D("NominalERecDist", ";ERec;Count", 100, 0, 10);
  std::vector<TH1D *> ThrowDistributions;
  for (UInt_t t_it = 0; t_it < NThrows; ++t_it) {
    ThrowDistributions.push_back(new TH1D("throw_dist", "", 100, 0, 10));
    ThrowDistributions.back()->SetDirectory(nullptr);
  }

  for (UInt_t e_it = 0; e_it < NEntries; ++e_it) {
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
    } else { // Do throws;

      for (UInt_t t_it = 0; t_it < NThrows; ++t_it) {
        xsecweights[t_it] = 1;
        for (VALORModel::TrueClass d : appdial) {
          double w = GetVALORWeight(
              d, ParameterThrows[t_it][static_cast<int>(d)], edr);
          if (w < 0) {
            std::cout << "[INFO]: Saw negative weight for dial: " << d
                      << ", set at "
                      << ParameterThrows[t_it][static_cast<int>(d)]
                      << " gave a weight = " << w << std::endl;
          }
          xsecweights[t_it] *= w;
          if (w != w || (w && !std::isnormal(w))) {
            std::cout << "[ERROR]: Found bad weight " << w << ", throw " << t_it
                      << ", dial " << d << ", value "
                      << ParameterThrows[t_it][static_cast<int>(d)]
                      << std::endl;
          }
        }
        ThrowDistributions[t_it]->Fill(
            edr.GetProjection(DepositsSummary::kERec), xsecweights[t_it]);
        MeanThrow->Fill(edr.GetProjection(DepositsSummary::kERec),
                        xsecweights[t_it]);
      }
    }

    NFills++;
    friendtree->Fill();
  }

  CovarianceBuilder cb(MeanThrow->GetXaxis()->GetNbins());

  MeanThrow->Scale(1.0 / double(NThrows));
  cb.SetMean(MeanThrow);

  for (UInt_t t_it = 0; t_it < NThrows; ++t_it) {
    cb.AddThrow_CovMatCalc(ThrowDistributions[t_it]);
  }
  cb.FinalizeCovMatCalc();

  cb.GetCovMatrix()->Clone()->Write("ERec_CovMat");
  cb.GetCorrMatrix()->Clone()->Write("ERec_CorrMat");

  for (Int_t i = 0; i < MeanThrow->GetXaxis()->GetNbins(); ++i) {
    double binval = MeanThrow->GetBinContent(i + 1);
    if (binval != binval || (binval && !std::isnormal(binval))) {
      std::cout << "[ERROR]: Attempted to add mean vector with bad value: "
                << binval << " at index " << i << std::endl;
    }
    double cval = (*cb.GetCovMatrix())[i][i];
    if (cval < 0) {
      std::cout
          << "[ERROR]: Found negative diagonal value in covariance matrix."
          << std::endl;
      throw;
    }
    if (cval && std::isnormal(cval)) {
      MeanThrow->SetBinError(i + 1, sqrt(cval));
    } else {
      if (cval) {
        std::cout << "[WARN]: Found odd covmat value for bin " << i << " = "
                  << cval << std::endl;
      }
      MeanThrow->SetBinError(i + 1, 0);
    }
  }

  double MaxThrow = -std::numeric_limits<double>::max();
  int MaxIndex = 0;
  double MinThrow = -std::numeric_limits<double>::max();
  int MinIndex = 0;

  for (UInt_t t_it = 0; t_it < NThrows; ++t_it) {
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
      (std::string("MaxERecDist_throw_") + to_str(MaxIndex)).c_str());
  ThrowDistributions[MinIndex]->SetName(
      (std::string("MinERecDist_throw_") + to_str(MinIndex)).c_str());
  ThrowDistributions[MaxIndex]->SetDirectory(of);
  ThrowDistributions[MinIndex]->SetDirectory(of);

  UncertMatrix->Write("PriorUncerts");

  paramThrowsCB.FinalizeCovMatCalc();
  paramThrowsCB.GetCovMatrix()->Clone()->Write("ParamThrows_CovMat");
  paramThrowsCB.GetCorrMatrix()->Clone()->Write("ParamThrows_CorrMat");

  of->Write();
  of->Close();
}
