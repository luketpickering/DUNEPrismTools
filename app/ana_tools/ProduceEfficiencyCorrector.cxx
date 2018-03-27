#include "BoundingBox.h"
#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::vector<DetectorStop> DetectorStops;
std::vector<BoundingBox> BBs;

std::string inpfile;
std::string oupfile;

double HadrVeto = 10E-3;     // GeV
double SelMuExitKE = 50E-3;  // GeV

std::vector<double> ERecBinning;
std::vector<double> XRangeBins;
bool HaveXRangeBins = false;
std::vector<double> Coeffs;
bool BuildMissingProtonEFakeData = false;

std::string FluxFitFile = "";
void GetXRangeBins() {
  TChain *CoeffTree = OpenTChainWithFileList("CoeffTree", FluxFitFile);

  if (!CoeffTree) {
    exit(1);
  }

  double XRange[2];
  double Coeff;

  CoeffTree->SetBranchAddress("XRange", &XRange);
  CoeffTree->SetBranchAddress("Coeff", &Coeff);

  std::cout << "[INFO]: XRange bins: " << std::flush;
  CoeffTree->GetEntry(0);
  XRangeBins.push_back(XRange[0]);
  Coeffs.push_back(Coeff);
  XRangeBins.push_back(XRange[1]);
  std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
  for (Long64_t i = 1; i < CoeffTree->GetEntries(); ++i) {
    CoeffTree->GetEntry(i);

    // If non-contiguous, must push an empty bit between.
    if (fabs(XRangeBins.back() - XRange[0]) > 1E-5) {
      Coeffs.push_back(0);
      XRangeBins.push_back(XRange[0]);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }

    Coeffs.push_back(Coeff);
    XRangeBins.push_back(XRange[1]);
    std::cout << ", " << XRangeBins.back() << std::flush;
  }
  std::cout << std::endl;

  delete CoeffTree;
  HaveXRangeBins = true;
}

bool OneBin = false;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>    : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV {default: 10}.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV "
         "{default: 0}.\n"
         "\t-F <FluxFitFile.root>       : Results of a flux fit that can be "
         "used to produce a perfect efficiency correction in relevant off-axis "
         "position bins.\n"
         "\t-OB                         : Efficiency correct in a single "
         "absolute position bin. Used to produce perfect efficiency correction "
         "at the far detector.\n"
         "\t-b <binning descriptor>     : Energy binning descriptor that can "
         "be used to produce a perfect efficiency correction in the relevant "
         "off-axis bins.\n"
         "\t-FDproton                   : Build missing proton energ fake "
         "data. This is very hard coded, if you don't know what it is, don't "
         "use it.\n"
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
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-F") {
      FluxFitFile = argv[++opt];
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
    } else if (std::string(argv[opt]) == "-FDproton") {
      BuildMissingProtonEFakeData = true;
    } else if (std::string(argv[opt]) == "-OB") {
      OneBin = true;
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

  if (FluxFitFile.size()) {
    GetXRangeBins();
  }

  if (!ERecBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.1"};
    handleOpts(argc_dum, argv_dum);
  }

  TChain *config_in = OpenTChainWithFileList("configTree", inpfile);
  Int_t NDets;
  Double_t FVGap[3];
  config_in->SetBranchAddress("NStops", &NDets);
  config_in->SetBranchAddress("FVGap", &FVGap);
  config_in->GetEntry(0);

  TChain *stopConfig_in = OpenTChainWithFileList("stopConfigTree", inpfile);
  Double_t StopMin[3], StopMax[3];
  Double_t Offset;
  stopConfig_in->SetBranchAddress("Min", &StopMin);
  stopConfig_in->SetBranchAddress("Max", &StopMax);
  stopConfig_in->SetBranchAddress("Offset", &Offset);

  Double_t MaxStopWidth = -std::numeric_limits<double>::max();

  Double_t MaxToWall = -std::numeric_limits<double>::max();

  for (Int_t d_it = 0; d_it < NDets; ++d_it) {
    stopConfig_in->GetEntry(d_it);
    TVector3 Min, Max;

    Min[0] = StopMin[0];
    Min[1] = StopMin[1];
    Min[2] = StopMin[2];

    Max[0] = StopMax[0];
    MaxStopWidth = std::max(MaxStopWidth, (Max[0] - Min[0]) - 2 * FVGap[0]);
    Max[1] = StopMax[1];
    Max[2] = StopMax[2];

    BBs.emplace_back(Max, Min);

    MaxToWall =
        std::max((Max - Min - 2 * TVector3(FVGap[0], FVGap[1], FVGap[2])).Mag(),
                 MaxToWall);

    std::cout << "[INFO]: Read det stop " << d_it << " {" << Min[0] << ", "
              << Min[1] << ", " << Min[2] << "} -- {" << Max[0] << ", "
              << Max[1] << ", " << Max[2] << "}" << std::endl;
  }

  if (OneBin) {
    XRangeBins.push_back(-MaxStopWidth / 2.0);
    XRangeBins.push_back(MaxStopWidth / 2.0);
    HaveXRangeBins = true;
  }

  std::cout << "[INFO]: Max DetXRange = " << MaxStopWidth
            << ", Max ToWall = " << (MaxToWall * 1E-2) << std::endl;

  EDep edr("EDeps", inpfile);

  Double_t ProtonFakeDataWeight = 1;
  if (BuildMissingProtonEFakeData) {
    // Should probably check that this exists, but leaving hardcoded and
    // presumptive for now as CheckTTreeHasBranch doesn't seem to work on
    // pre-friended trees.
    edr.tree->SetBranchAddress("EnuTp", &ProtonFakeDataWeight);
  }

  TFile *of = CheckOpenFile(oupfile, "RECREATE");

  TH2D *MuonKinematics_all =
      new TH2D("MuonKinematics_all", ";#it{E}_{#mu};ToWall_{#mu} (m);Count",
               100, 0, 20, 60, 0, MaxToWall * 1E-2);
  TH2D *MuonKinematics_musel =
      new TH2D("MuonKinematics_musel", ";#it{E}_{#mu};ToWall_{#mu} (m);Count",
               100, 0, 20, 60, 0, MaxToWall * 1E-2);
  TH2D *MuonKinematics_seleff =
      new TH2D("MuonKinematics_eff", ";#it{E}_{#mu};ToWall_{#mu} (m);#epsilon",
               100, 0, 20, 60, 0, MaxToWall * 1E-2);

  TH2D *Ehadr_FV_detpos_all = new TH2D(
      "Ehadr_FV_detpos_all", ";#it{E}_{Hadr};Detector X position (cm);Count",
      100, 0, 20, 100, -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));
  TH2D *Ehadr_FV_detpos_hadrsel =
      new TH2D("Ehadr_FV_detpos_hadrsel",
               ";#it{E}_{Hadr};Detector X position (cm);Count", 100, 0, 20, 100,
               -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));
  TH2D *Ehadr_FV_detpos_seleff =
      new TH2D("Ehadr_FV_detpos_seleff",
               ";#it{E}_{Hadr};Detector X position (cm);#epsilon", 100, 0, 20,
               100, -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));

  TH2D *Eproxy_FV_abspos_all = nullptr;
  TH2D *Eproxy_FV_abspos_hadrsel = nullptr;
  TH2D *Eproxy_FV_abspos_seleff = nullptr;

  if (HaveXRangeBins) {
    Eproxy_FV_abspos_all =
        new TH2D("Eproxy_FV_abspos_all",
                 ";#it{E}_{#nu,Proxy};Off-axis position (cm);Count",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    Eproxy_FV_abspos_hadrsel =
        new TH2D("Eproxy_FV_abspos_hadrsel",
                 ";#it{E}_{#nu,Proxy};Off-axis position (cm);Count",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    Eproxy_FV_abspos_seleff =
        new TH2D("Eproxy_FV_abspos_seleff",
                 ";#it{E}_{#nu,Proxy};Off-axis position (cm);#epsilon",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
  }

  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);
    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;

      if (BuildMissingProtonEFakeData) {
        std::cout << "\tProton fake data weight = " << ProtonFakeDataWeight
                  << std::endl;
      }
    }

    if ((edr.stop == -1) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    if (BuildMissingProtonEFakeData) {
      edr.TotalNonlep_Dep_veto -= edr.ProtonDep_veto * 0.2;
      edr.ERecProxy_True -= edr.EKinProton_True * 0.2;
    }

    TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
    TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                   edr.PrimaryLep_4mom[2]);
    TrDir = TrDir.Unit();

    double ToWall = CalculateToWall(BBs[edr.stop], TrStr, TrDir) * 1E-2;

    MuonKinematics_all->Fill(edr.PrimaryLep_4mom[3], ToWall,
                             ProtonFakeDataWeight);
    MuonKinematics_musel->Fill(
        edr.PrimaryLep_4mom[3], ToWall,
        ProtonFakeDataWeight * double(edr.LepExitKE > SelMuExitKE));

    Ehadr_FV_detpos_all->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                              edr.vtxInDetX, ProtonFakeDataWeight);
    Ehadr_FV_detpos_hadrsel->Fill(
        edr.ERecProxy_True - edr.PrimaryLep_4mom[3], edr.vtxInDetX,
        ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));

    if (HaveXRangeBins) {
      Eproxy_FV_abspos_all->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                                 ProtonFakeDataWeight);
      Eproxy_FV_abspos_hadrsel->Fill(
          edr.ERecProxy_True, -1 * edr.vtx[0],
          ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));
    }
  }

  MuonKinematics_seleff->Divide(MuonKinematics_musel, MuonKinematics_all);
  Ehadr_FV_detpos_seleff->Divide(Ehadr_FV_detpos_hadrsel, Ehadr_FV_detpos_all);
  if (HaveXRangeBins) {
    Eproxy_FV_abspos_seleff->Divide(Eproxy_FV_abspos_hadrsel,
                                    Eproxy_FV_abspos_all);
  }

  of->Write();
  of->Close();
}
