#include "FullDetTreeReader.h"
#include "G4ArReader.h"
#include "Utils.hxx"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

std::vector<double> detmin;
std::vector<double> detmax;
std::vector<double> fvgap;
int nxsteps = 400;
int ntrackingsteps = 1000;
std::string inputG4ArFileName;
std::string inputRooTrackerFileName;
std::string outputFileName;
Long64_t nmaxevents = std::numeric_limits<int>::max();

double IntegratedPOT = 0;

// #define DEBUG

void SayUsage(char* argv[]) {
  std::cout << "[INFO]: Use like: " << argv[0]
            << " -i <inputg4arbofile> -ir <inputGENIERooTrackerfile> -dmn "
               "<detxmin,ymin,zmin> -dmx <detxmax,ymax,zmax> -fv <fidgapx,y,z> "
               "-o <outputfile> [-P <integratedPOT> -nx <nxsteps> -nt "
               "<ntrackingsteps> -n <nmaxevents>]"
            << std::endl;
}

void handleOpts(int argc, char* argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-nx") {
      nxsteps = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-nt") {
      ntrackingsteps = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-n") {
      nmaxevents = str2T<Long64_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-dmn") {
      detmin = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-dmx") {
      detmax = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-fv") {
      fvgap = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-i") {
      inputG4ArFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-ir") {
      inputRooTrackerFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-P") {
      IntegratedPOT = str2T<double>(argv[++opt]);
    } else if ((std::string(argv[opt]) == "-?") ||
               std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char* argv[]) {
  handleOpts(argc, argv);

  DetectorAndFVDimensions detdims;
  detdims.NXSteps = nxsteps;
  for (size_t dim_it = 0; dim_it < 3; ++dim_it) {
    detdims.DetMin[dim_it] = detmin[dim_it];
    detdims.DetMax[dim_it] = detmax[dim_it];
    detdims.FVGap[dim_it] = fvgap[dim_it];
  }

  G4ArReader g4ar(inputG4ArFileName, detdims, inputRooTrackerFileName);

  // WriteConfigTree
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");

  Int_t NMaxTrackSteps = ntrackingsteps;
  TTree* config = new TTree("configTree", "Run configuration tree");
  config->Branch("NXSteps", &detdims.NXSteps, "NXSteps/I");
  config->Branch("DetMin", &detdims.DetMin, "DetMin[3]/D");
  config->Branch("DetMax", &detdims.DetMax, "DetMax[3]/D");
  config->Branch("FVGap", &detdims.FVGap, "FVGap[3]/D");
  config->Branch("NMaxTrackSteps", &NMaxTrackSteps, "NMaxTrackSteps/I");
  config->Branch("IntegratedPOT", &IntegratedPOT, "IntegratedPOT/D");
  config->Fill();

  g4ar.SetNMaxTrackSteps(NMaxTrackSteps);

  // Build coalescer deposits and initialize detector histograms
  DepoParticle LepDep;
  LepDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle HadDep;
  HadDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle ProtonDep;
  ProtonDep.LendDepositMaps(detdims.BuildDetectorMap(),
                            detdims.BuildDetectorMap());
  DepoParticle NeutronDep;
  NeutronDep.SetTrackTime();
  NeutronDep.LendDepositMaps(detdims.BuildDetectorMap(),
                             detdims.BuildDetectorMap());
  DepoParticle PiCDep;
  PiCDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle Pi0Dep;
  Pi0Dep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle OtherDep;
  OtherDep.LendDepositMaps(detdims.BuildDetectorMap(),
                           detdims.BuildDetectorMap());

  TTree* fdTree =
      new TTree("fulldetTree", "G4 and GENIE passthrough information");

  FullDetTreeReader* fdw =
      FullDetTreeReader::MakeTreeWriter(fdTree, nxsteps, NMaxTrackSteps);

  g4ar.ResetCurrentEntry();
  g4ar.TrackTimeForPDG(2112);
  int evnum = 0;
  int loudevery = std::min(nmaxevents, g4ar.NInputEntries) / 10;
  int nfills = 0;
  do {
    Event ev = g4ar.BuildEvent();
    fdw->Reset();

    LepDep.Reset();
    HadDep.Reset();
    ProtonDep.Reset();
    NeutronDep.Reset();
    PiCDep.Reset();
    Pi0Dep.Reset();
    OtherDep.Reset();

    // ====================== START GENIE Pass through ======================

    // Fill primary interaction info
    fdw->EventNum = ev.ev_id;

    (*fdw->EventCode) = (*ev.RooTrackerInteractionCode);

    fdw->VertexPosition[0] = ev.VertexPosition[0];
    fdw->VertexPosition[1] = ev.VertexPosition[1];
    fdw->VertexPosition[2] = ev.VertexPosition[2];

    int nu_pdgs[] = {12, -12, 14, -14};
    PrimaryParticle nu = ev.GetFirstPrimaryWithPDG(nu_pdgs, false);
    fdw->nu_4mom[0] = nu.ThreeMom[0];
    fdw->nu_4mom[1] = nu.ThreeMom[1];
    fdw->nu_4mom[2] = nu.ThreeMom[2];
    fdw->nu_4mom[3] = nu.EKin;

    fdw->nu_PDG = nu.PDG;

    int fslep_pdgs[] = {12, -12, 14, -14, 11, -11, 13, -13};
    PrimaryParticle fslep = ev.GetFirstPrimaryWithPDG(fslep_pdgs, true);

    fdw->PrimaryLepPDG = fslep.PDG;

    fdw->PrimaryLep_4mom[0] = fslep.ThreeMom[0];
    fdw->PrimaryLep_4mom[1] = fslep.ThreeMom[1];
    fdw->PrimaryLep_4mom[2] = fslep.ThreeMom[2];
    fdw->PrimaryLep_4mom[3] = (fslep.EKin + fslep.EMass);

    fdw->FourMomTransfer_True[0] = (nu.ThreeMom[0] - fslep.ThreeMom[0]);
    fdw->FourMomTransfer_True[1] = (nu.ThreeMom[1] - fslep.ThreeMom[1]);
    fdw->FourMomTransfer_True[2] = (nu.ThreeMom[2] - fslep.ThreeMom[2]);
    fdw->FourMomTransfer_True[3] = (nu.EKin - (fslep.EKin + fslep.EMass));

    TLorentzVector FourMomTransf(
        fdw->FourMomTransfer_True[0], fdw->FourMomTransfer_True[1],
        fdw->FourMomTransfer_True[2], fdw->FourMomTransfer_True[3]);

    fdw->Q2_True = -FourMomTransf.Mag2();

    fdw->y_True = 1 - ((fslep.EKin + fslep.EMass) / nu.EKin);
    double nucleon_mass_GeV = .93827208;
    fdw->W_Rest =
        sqrt(fdw->Q2_True +
             2.0 * nucleon_mass_GeV * (nu.EKin - (fslep.EKin + fslep.EMass)) +
             nucleon_mass_GeV * nucleon_mass_GeV);

    for (PrimaryParticle& p : ev.PrimaryParticles) {
      if (!p.IsFinalState) {
        continue;
      }
      switch (abs(p.PDG)) {
        case 11:
        case 12:
        case 13:
        case 14: {
          fdw->NLep++;
          break;
        }
        case 111: {
          fdw->NPi0++;
          fdw->EKinPi0_True += p.EKin;
          fdw->EMassPi0_True += p.EMass;
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
          break;
        }
        case 211: {
          fdw->NPiC++;
          fdw->EKinPiC_True += p.EKin;
          fdw->EMassPiC_True += p.EMass;
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
          break;
        }
        case 2212: {
          fdw->NProton++;
          fdw->EKinProton_True += p.EKin;
          fdw->EMassProton_True += p.EMass;
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
          break;
        }
        case 2112: {
          fdw->NNeutron++;
          fdw->EKinNeutron_True += p.EKin;
          fdw->EMassNeutron_True += p.EMass;
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
          break;
        }
        case 22: {
          fdw->NGamma++;
          fdw->EGamma_True += p.EKin;
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
          break;
        }
        default: {
          fdw->NOther++;
          fdw->EOther_True += (p.EKin + p.EMass);
          fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
#ifdef DEBUG
          std::cout << "[INFO]: NOther PDG = " << p.PDG << std::endl;
#endif
        }
      }
    }

    // ====================== END GENIE Pass through ========================

    // ===================== START GEANT4 Pass through ======================

    for (DepoTracked& td : ev.TrackedDeposits) {
      if (abs(td.PDG) == 13) {
        LepDep.AddDeposit(td);

        // Fill tracking info
        fdw->NMuonTrackSteps = td.NSteps;
        std::copy_n(td._Position, td.NSteps * 3, fdw->MuonTrackPos_1D);
        std::copy_n(td._Momentum, td.NSteps * 3, fdw->MuonTrackMom_1D);

#ifdef DEBUG
        std::cout << "[INFO]: Event had " << fdw->NMuonTrackSteps
                  << " tracked muon steps." << std::endl;
        for (Int_t i = 0; i < fdw->NMuonTrackSteps; ++i) {
          std::cout << "[INFO]: Step [" << i << "] at {"
                    << fdw->MuonTrackPos[i][0] << ", "
                    << fdw->MuonTrackPos[i][1] << ", "
                    << fdw->MuonTrackPos[i][2] << " } with 3-mom {"
                    << fdw->MuonTrackMom[i][0] << ", "
                    << fdw->MuonTrackMom[i][1] << ", "
                    << fdw->MuonTrackMom[i][2] << " }." << std::endl;
        }
#endif
      }
    }

    for (DepoParticle& td : ev.TotalDeposits) {
      switch (abs(td.PDG)) {
        case 13:
        case 11: {
          LepDep.AddDeposit(td);
          break;
        }
        case 2212: {
          ProtonDep.AddDeposit(td);
          HadDep.AddDeposit(td);
          break;
        }
        case 2112: {
          NeutronDep.AddDeposit(td);
          HadDep.AddDeposit(td);
          break;
        }
        case 211: {
          PiCDep.AddDeposit(td);
          HadDep.AddDeposit(td);
          break;
        }
        case 111: {
          Pi0Dep.AddDeposit(td);
          HadDep.AddDeposit(td);
          break;
        }
        default: {
          OtherDep.AddDeposit(td);
          HadDep.AddDeposit(td);
          break;
        }
      }
    }

    // ====================== END GEANT4 Pass through =======================

    // Fill output branches
    for (Int_t x_it = 0; x_it < detdims.NXSteps; ++x_it) {
      for (size_t y_it = 0; y_it < 3; ++y_it) {
        for (size_t z_it = 0; z_it < 3; ++z_it) {
          Int_t gbin = LepDep.Deposits->GetBin(x_it + 1, y_it + 1, z_it + 1);

#ifdef DEBUG
          if (LepDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", LepDep content = "
                      << LepDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (HadDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", HadDep content = "
                      << HadDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (ProtonDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", ProtonDep content = "
                      << ProtonDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (NeutronDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", NeutronDep content = "
                      << NeutronDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (PiCDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", PiCDep content = "
                      << PiCDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (Pi0Dep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", Pi0Dep content = "
                      << Pi0Dep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (OtherDep.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", OtherDep content = "
                      << OtherDep.Deposits->GetBinContent(gbin) << std::endl;
          }
          if (LepDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", LepDaughtDep content = "
                      << LepDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (HadDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", HadDaughtDep content = "
                      << HadDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (ProtonDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", ProtonDaughtDep content = "
                      << ProtonDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (NeutronDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", NeutronDaughtDep content = "
                      << NeutronDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (PiCDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", PiCDaughtDep content = "
                      << PiCDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (Pi0Dep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", Pi0DaughtDep content = "
                      << Pi0Dep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }
          if (OtherDep.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", OtherDaughtDep content = "
                      << OtherDep.DaughterDeposits->GetBinContent(gbin)
                      << std::endl;
          }

#endif

          fdw->LepDep[x_it][y_it][z_it] = LepDep.Deposits->GetBinContent(gbin);
          fdw->HadDep[x_it][y_it][z_it] = HadDep.Deposits->GetBinContent(gbin);
          fdw->ProtonDep[x_it][y_it][z_it] =
              ProtonDep.Deposits->GetBinContent(gbin);
          fdw->NeutronDep[x_it][y_it][z_it] =
              NeutronDep.Deposits->GetBinContent(gbin);
          fdw->NeutronDep_ChrgWSumTime[x_it][y_it][z_it] =
              NeutronDep.Deposits_ChrgWSumTime->GetBinContent(gbin);
          fdw->PiCDep[x_it][y_it][z_it] = PiCDep.Deposits->GetBinContent(gbin);
          fdw->Pi0Dep[x_it][y_it][z_it] = Pi0Dep.Deposits->GetBinContent(gbin);
          fdw->OtherDep[x_it][y_it][z_it] =
              OtherDep.Deposits->GetBinContent(gbin);

          fdw->LepDaughterDep[x_it][y_it][z_it] =
              LepDep.DaughterDeposits->GetBinContent(gbin);
          fdw->HadDaughterDep[x_it][y_it][z_it] =
              HadDep.DaughterDeposits->GetBinContent(gbin);
          fdw->ProtonDaughterDep[x_it][y_it][z_it] =
              ProtonDep.DaughterDeposits->GetBinContent(gbin);
          fdw->NeutronDaughterDep[x_it][y_it][z_it] =
              NeutronDep.DaughterDeposits->GetBinContent(gbin);
          fdw->NeutronDaughterDep_ChrgWSumTime[x_it][y_it][z_it] =
              NeutronDep.DaughterDeposits_ChrgWSumTime->GetBinContent(gbin);
          fdw->PiCDaughterDep[x_it][y_it][z_it] =
              PiCDep.DaughterDeposits->GetBinContent(gbin);
          fdw->Pi0DaughterDep[x_it][y_it][z_it] =
              Pi0Dep.DaughterDeposits->GetBinContent(gbin);
          fdw->OtherDaughterDep[x_it][y_it][z_it] =
              OtherDep.DaughterDeposits->GetBinContent(gbin);
        }
      }
    }

    fdTree->Fill();
    nfills++;

    if (loudevery && !(evnum % loudevery)) {
      std::cout << "[INFO]: Processed " << evnum << " entries." << std::endl;
    }

    evnum++;
  } while (g4ar.GetNextEvent() && (evnum < nmaxevents));

  std::cout << "[INFO]: Filled the output tree " << nfills << " times."
            << std::endl;

  LepDep.DeleteDepositMaps();
  HadDep.DeleteDepositMaps();
  ProtonDep.DeleteDepositMaps();
  NeutronDep.DeleteDepositMaps();
  PiCDep.DeleteDepositMaps();
  Pi0Dep.DeleteDepositMaps();
  OtherDep.DeleteDepositMaps();

  outfile->Write();
  outfile->Close();
}
