#include "EDepTreeReader.h"
#include "VALORModelClassifier.h"
#include "CovarianceHelper.h"

#include "Utils.hxx"

#include "omp.h"

#include <string>
#include <vector>

std::string inpfile;
std::string oupfile;

std::string CovmatFile, CovmatHist;
Long64_t NThrows = 1;


void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>   : TChain descriptor for"
               " input tree. \n"
               "\t-o <friendoutput.root>     : TFile for output friend tree. \n"
               "\t-C <covmatfile.root,hname> : Input covariance matrix to "
               "throw xsec weights from.\n"
               "\t-N <NThrows>     : Number of throws to make. Special cases "
               "are \n\t\t\t  1: will produce +1 sigma variation\n\t\t\t  2: "
               "will produce +1 and -1 sigma variation.\n"
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
    } else if (std::string(argv[opt]) == "-o") {
      oupfile = argv[++opt];
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

  std::vector<std::vector<double> > ParameterThrows;
  ParameterThrows.resize(NThrows);
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    ParameterThrows[t_it].resize(
        static_cast<int>(VALORModel::TrueClass::kNVALORDials));
  }

  CovarianceThrower ct(*UncertMatrix);

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    TMatrixD const *CovVector = ct.Throw();
    for (int p_it = 0;
         p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {
      ParameterThrows[t_it][p_it] = (*CovVector)[p_it][0];
    }
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

  EDep edr("EDeps", inpfile);

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");
  TTree *friendtree = new TTree("XSecWeights", "");
  TTree *configtree = new TTree("ConfigTree", "");

  configtree->Branch("NThrows", &NThrows, "NThrows/I");
  configtree->Fill();

  double *xsecweights = new double[NThrows];
  friendtree->Branch(
      "XSecWeights", xsecweights,
      (std::string("XSecWeights[") + to_str(NThrows) + "]/D").c_str());

  size_t loud_every = edr.GetEntries() / 10;
  size_t NFills = 0;
  Long64_t NEntries = edr.GetEntries();

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
      }
    }

    NFills++;
    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
