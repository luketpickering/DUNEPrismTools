#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "BargerPropagator.h"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::string inpfile;
std::string outputfile = "FarDetectorPrediction.root";
bool UseTrueENu = true;
double DipAngle = 5.8;
double OscParams[6] = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};
static const double deg2rad = asin(1) / 90.0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>             : TChain descriptor for"
               " input tree. \n"
               "\t-p "
               "<sin2(theta12)=0.825>,<sin2(theta13)=0.10>,<sin2(theta23)=1.0>,"
               "<dm12=7.9e-5>,<dm23=2.5e-3>,<dcp=0.0>\n"
               "\t-d <dipangle=5.8>\n"
               "\t-o <outputfile.root>                 : Output file name.\n"
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

  TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");

  TH1D *FDPrediction = new TH1D("FDPrediction", "", 150, 0, 15);
  TH1D *FDPrediction_True = new TH1D("FDPrediction_True", "", 150, 0, 15);

  TH1D *FDPredictionOsc = new TH1D("FDPredictionOsc", "", 150, 0, 15);
  TH1D *FDPredictionOsc_True = new TH1D("FDPredictionOsc_True", "", 150, 0, 15);

  TH2D *ERecFD =
      new TH2D("ERecFD", ";ETrue;ERec;Count", 150, 0, 15, 150, 0, 15);
  TH2D *ERecFDosc =
      new TH2D("ERecFDOsc", ";ETrue;ERec;Count", 150, 0, 15, 150, 0, 15);

  nuTypes numu = GetNuType(14);

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    if ((edr.stop < 0) || (!edr.HadrShowerContainedInFV) ||
        (edr.PrimaryLepPDG != 13)) {
      continue;
    }

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

    ERecFD->Fill(edr.nu_4mom[3],
                 edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                     edr.TotalNonlep_Dep_veto);
    ERecFDosc->Fill(edr.nu_4mom[3],
                    edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                        edr.TotalNonlep_Dep_veto,
                    oweight);
  }

  outfile->Write();
  outfile->Close();
}
