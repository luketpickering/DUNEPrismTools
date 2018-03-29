#include "CondensedDepositsTreeReader.hxx"
#include "SimConfigTreeReader.hxx"

#include "G4ArReader.hxx"
#include "StopDimensions.hxx"

#include "ROOTUtility.hxx"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

std::vector<double> detmin;
std::vector<double> detmax;
std::vector<double> stopVetoGap;
int NXSteps = 0;
int NMaxTrackSteps = 1000;
std::string InputG4ArFileName;
std::string InputRooTrackerFileName;
std::string OutputFileName;
Long64_t NMaxEvents = std::numeric_limits<int>::max();
bool KeepOOFVYZEvents = false;
double timesep_us = 0xdeadbeef;
double POTPerFile = 0xdeadbeef;

// #define DEBUG

void SayUsage(char* argv[]) {
  std::cout << "[INFO]: Use like: " << argv[0]
            << "\n\t-i <inputg4arbofile>     : Input G4Py root file."
               "\n\t-ir <GENIERTFile>        : Input GENIE rootracker file "
               "(index sync'd with -i argument)"
               "\n\t-o <outputfile>          : File to write output to"
               "\n\t-dmn <detxmin,ymin,zmin> : Active detector minimum (cm)"
               "\n\t-dmx <detxmax,ymax,zmax> : Active detector maximum (cm)"
               "\n\t-V <vetogap x,y,z>       : Active veto region to pad each "
               "corresponding face of the "
               "\n\t                           active volume with."
               "\n\t                          -- N.B. If -dmn -1,-1,-1 "
               "-dmx 1,1,1 -V 0.5,0.5,0.5 "
               "\n\t                           then the non-veto volume will be"
               " 1x1x1 cm^{3} centered on "
               "\n\t                           the origin."
               "\n\t-P <POTPerFile>          : Adds POTPerFile information to "
               "the metadata, used for "
               "\n\t                           POT-normalising predicted event "
               "rates and Near/Far "
               "\n\t                           comparisons downstream."
               "\n\t-n <NMaxEvents>          : Run no more than -n events."
               "\n\t-A                       : Output all events even if they "
               "occured outside of the "
               "\n\t                           non-veto active volume."
               "\n\t-T                       : Will add timesep branches to the"
               " output that "
               "\n\t                           contain all deposits ocurring "
               "more than -T <timesep>"
               "\n\t                           microseconds after the neutrino "
               "interaction."
               "\n\t-nx <NXSteps>            : Number of x slices to break up "
               "total non-veto active "
               "\n\t                           region into. "
               "\n\t                          -- N.B. two steps will be added "
               "for the X veto gap "
               "\n\t                             passed to -V. "
               "\n\t                          -- N.B. The non-veto active "
               "region X dimension, "
               "\n\t                             i.e. -dmx less -dmn less times"
               " the -V should be "
               "\n\t                             easily divisible by this "
               "number: e.g. 10 cm steps."
               "\n\t-nt <NMaxTrackSteps>     : Track final state charged lepton"
               " through up to -nt GEANT4"
               "\n\t                           steps."
            << std::endl;
}

void handleOpts(int argc, char* argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-nx") {
      NXSteps = str2T<int>(argv[++opt]) + 2;
    } else if (std::string(argv[opt]) == "-nt") {
      NMaxTrackSteps = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-n") {
      NMaxEvents = str2T<Long64_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-dmn") {
      detmin = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-dmx") {
      detmax = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-V") {
      stopVetoGap = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-i") {
      InputG4ArFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-ir") {
      InputRooTrackerFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-P") {
      POTPerFile = str2T<double>(argv[++opt]);
      std::cout << "[INFO]: Using " << POTPerFile << " POT as the POTPerFile"
                << std::endl;
    } else if (std::string(argv[opt]) == "-T") {
      timesep_us = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-A") {
      KeepOOFVYZEvents = true;
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

  if(NXSteps == 0){
    std::cout << "[ERROR]: Must specify -nx option. For the standard "
    "configuration of: -dmn -3800,-150,-250 -dmx 200,150,250 -V 50,50,50, using"
    " -nx 390 is recommended, which corresponds to 10cm deposit slices."
    << std::endl;
    throw;
  }

  StopDimensions detdims;
  detdims.NXSteps = NXSteps;
  std::copy_n(detmin.data(),3,detdims.DetMin);
  std::copy_n(detmax.data(),3,detdims.DetMax);
  std::copy_n(stopVetoGap.data(),3,detdims.VetoGap);

  G4ArReader g4ar(InputG4ArFileName, detdims, InputRooTrackerFileName,
                  timesep_us);

  // WriteConfigTree
  TFile* outfile = CheckOpenFile(OutputFileName.c_str(), "RECREATE");

  TTree* SimConfigTree = new TTree("SimConfigTree", "Run configuration tree");

  SimConfig *sc = SimConfig::MakeTreeWriter(SimConfigTree);
  std::copy_n(detdims.DetMin,3,sc->DetMin);
  std::copy_n(detdims.DetMax,3,sc->DetMax);
  std::copy_n(detdims.VetoGap,3,sc->VetoGap);
  sc->NXSteps = detdims.NXSteps;
  sc->NMaxTrackSteps = NMaxTrackSteps;
  sc->POTPerFile = POTPerFile;
  sc->timesep_us = timesep_us;
  SimConfigTree->Fill();

  g4ar.SetNMaxTrackSteps(NMaxTrackSteps);

  // Build coalescer deposits and initialize detector histograms
  DepoParticle LepDep;
  LepDep.timesep_us = timesep_us;
  LepDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle HadDep;
  HadDep.timesep_us = timesep_us;
  HadDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle ProtonDep;
  ProtonDep.timesep_us = timesep_us;
  ProtonDep.LendDepositMaps(detdims.BuildDetectorMap(),
                            detdims.BuildDetectorMap());
  DepoParticle NeutronDep;
  NeutronDep.timesep_us = timesep_us;
  NeutronDep.SetTrackTime();
  NeutronDep.LendDepositMaps(detdims.BuildDetectorMap(),
                             detdims.BuildDetectorMap());
  DepoParticle PiCDep;
  PiCDep.timesep_us = timesep_us;
  PiCDep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle Pi0Dep;
  Pi0Dep.timesep_us = timesep_us;
  Pi0Dep.LendDepositMaps(detdims.BuildDetectorMap(),
                         detdims.BuildDetectorMap());
  DepoParticle OtherDep;
  OtherDep.timesep_us = timesep_us;
  OtherDep.LendDepositMaps(detdims.BuildDetectorMap(),
                           detdims.BuildDetectorMap());
  DepoParticle NuclRemDep;
  NuclRemDep.timesep_us = timesep_us;
  NuclRemDep.LendDepositMaps(detdims.BuildDetectorMap(),
                             detdims.BuildDetectorMap());

  TTree* fdTree =
      new TTree("CondensedDepositsTree",
      "G4 and GENIE passthrough information");

  CondensedDeposits* fdw = CondensedDeposits::MakeTreeWriter(
      fdTree, NXSteps, NMaxTrackSteps, timesep_us);

  g4ar.ResetCurrentEntry();
  g4ar.TrackTimeForPDG(2112);
  int evnum = 0;
  int loudevery = std::min(NMaxEvents, g4ar.NInputEntries) / 10;
  int nfills = 0;
  do {
    DepoEvent ev = g4ar.BuildEvent();

    if (!KeepOOFVYZEvents) {
      int ybin = LepDep.Deposits->GetYaxis()->FindFixBin(ev.VertexPosition[1]);
      int zbin = LepDep.Deposits->GetZaxis()->FindFixBin(ev.VertexPosition[2]);
      // Skips event if is not in YZ non-veto volume.
      if ((ybin != 2) || (zbin != 2)) {
        evnum++;
        continue;
      }
    }

    fdw->Reset();

    LepDep.Reset();
    HadDep.Reset();
    ProtonDep.Reset();
    NeutronDep.Reset();
    PiCDep.Reset();
    Pi0Dep.Reset();
    OtherDep.Reset();
    NuclRemDep.Reset();

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
        sqrt(2.0 * nucleon_mass_GeV * (nu.EKin - (fslep.EKin + fslep.EMass)) +
             nucleon_mass_GeV * nucleon_mass_GeV - fdw->Q2_True);

    double p4[4];
#ifdef DEBUG
    double TEnergy = 0;
    std::vector<std::string> ss;
#endif
    for (PrimaryParticle& p : ev.PrimaryParticles) {
      if (!p.IsFinalState) {
        continue;
      }
      p4[0] = p.ThreeMom[0];
      p4[1] = p.ThreeMom[1];
      p4[2] = p.ThreeMom[2];
      p4[3] = (p.EKin + p.EMass);
      fdw->AddPassthroughPart(p.PDG, p4);

#ifdef DEBUG

      TEnergy += ((abs(p.PDG) == 2112) || (abs(p.PDG) == 2212))
                     ? p.EKin
                     : (p.EKin + p.EMass);
      ss.emplace_back("");
      ss.back() += "\t" + to_str(p.PDG) + ", E addr = " +
                   to_str(((abs(p.PDG) == 2112) || (abs(p.PDG) == 2212))
                              ? p.EKin
                              : (p.EKin + p.EMass));
#endif

      fdw->TotalFS_3mom[0] += p.ThreeMom[0];
      fdw->TotalFS_3mom[1] += p.ThreeMom[1];
      fdw->TotalFS_3mom[2] += p.ThreeMom[2];

      if (G4ArReader::IsNuclearPDG(abs(p.PDG))) {
        fdw->KENuclearRemnant_True += p.EKin;
      } else {
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
            if (p.PDG < 0) {
              fdw->NAntiNucleons++;
            }
            fdw->NProton++;
            fdw->EKinProton_True += p.EKin;
            fdw->EMassProton_True += p.EMass;
            fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
            break;
          }
          case 2112: {
            if (p.PDG < 0) {
              fdw->NAntiNucleons++;
            }
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
            if ((abs(p.PDG) > 1000) && (abs(p.PDG) < 9999)) {
              fdw->NBaryonicRes++;
            } else {
              fdw->NOther++;
            }
            fdw->EOther_True += (p.EKin + p.EMass);
            fdw->ENonPrimaryLep_True += (p.EKin + p.EMass);
#ifdef DEBUG
            std::cout << "[INFO]: NOther PDG = " << p.PDG << std::endl;
#endif
          }
        }
      }
    }

#ifdef DEBUG
    if ((fdw->nu_4mom[3] - 1) > TEnergy) {
      g4ar.ShoutRooTracker();
      std::cout << fdw->EventCode->GetString() << std::endl;
      std::cout << "[INFO]: Neutrino E = " << fdw->nu_4mom[3]
                << ", total FS = " << TEnergy << std::endl;
      for (auto& s : ss) {
        std::cout << s << std::endl;
      }
    }

    ss.clear();
#endif

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
      if (G4ArReader::IsNuclearPDG(abs(td.PDG))) {
        NuclRemDep.AddDeposit(td);
      } else {
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
          case 22: {
            OtherDep.AddDeposit(td);
            break;
          }
          default: {
            OtherDep.AddDeposit(td);
            HadDep.AddDeposit(td);
            break;
          }
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

          if (timesep_us != 0xdeadbeef) {
            if (LepDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", LepDep_timesep content = "
                        << LepDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (HadDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", HadDep_timesep content = "
                        << HadDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (ProtonDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", ProtonDep_timesep content = "
                        << ProtonDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (NeutronDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", NeutronDep_timesep content = "
                        << NeutronDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (PiCDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", PiCDep_timesep content = "
                        << PiCDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (Pi0Dep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", Pi0Dep_timesep content = "
                        << Pi0Dep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (OtherDep.Deposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", OtherDep_timesep content = "
                        << OtherDep.Deposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (LepDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", LepDaughtDep_timesep content = "
                        << LepDep.DaughterDeposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (HadDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", HadDaughtDep_timesep content = "
                        << HadDep.DaughterDeposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (ProtonDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", ProtonDaughtDep_timesep content = "
                        << ProtonDep.DaughterDeposits_timesep->GetBinContent(
                               gbin)
                        << std::endl;
            }
            if (NeutronDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", NeutronDaughtDep_timesep content = "
                        << NeutronDep.DaughterDeposits_timesep->GetBinContent(
                               gbin)
                        << std::endl;
            }
            if (PiCDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", PiCDaughtDep_timesep content = "
                        << PiCDep.DaughterDeposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (Pi0Dep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", Pi0DaughtDep_timesep content = "
                        << Pi0Dep.DaughterDeposits_timesep->GetBinContent(gbin)
                        << std::endl;
            }
            if (OtherDep.DaughterDeposits_timesep->GetBinContent(gbin)) {
              std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", "
                        << z_it << ", OtherDaughtDep_timesep content = "
                        << OtherDep.DaughterDeposits_timesep->GetBinContent(
                               gbin)
                        << std::endl;
            }
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
          fdw->NuclRemDep[x_it][y_it][z_it] =
              NuclRemDep.Deposits->GetBinContent(gbin) +
              NuclRemDep.DaughterDeposits->GetBinContent(gbin);

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

          if (timesep_us != 0xdeadbeef) {
            fdw->LepDep_timesep[x_it][y_it][z_it] =
                LepDep.Deposits_timesep->GetBinContent(gbin);
            fdw->HadDep_timesep[x_it][y_it][z_it] =
                HadDep.Deposits_timesep->GetBinContent(gbin);
            fdw->ProtonDep_timesep[x_it][y_it][z_it] =
                ProtonDep.Deposits_timesep->GetBinContent(gbin);
            fdw->NeutronDep_timesep[x_it][y_it][z_it] =
                NeutronDep.Deposits_timesep->GetBinContent(gbin);
            fdw->PiCDep_timesep[x_it][y_it][z_it] =
                PiCDep.Deposits_timesep->GetBinContent(gbin);
            fdw->Pi0Dep_timesep[x_it][y_it][z_it] =
                Pi0Dep.Deposits_timesep->GetBinContent(gbin);
            fdw->OtherDep_timesep[x_it][y_it][z_it] =
                OtherDep.Deposits_timesep->GetBinContent(gbin);
            fdw->NuclRemDep[x_it][y_it][z_it] +=
                NuclRemDep.Deposits_timesep->GetBinContent(gbin) +
                NuclRemDep.DaughterDeposits_timesep->GetBinContent(gbin);

            fdw->LepDaughterDep_timesep[x_it][y_it][z_it] =
                LepDep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->HadDaughterDep_timesep[x_it][y_it][z_it] =
                HadDep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->ProtonDaughterDep_timesep[x_it][y_it][z_it] =
                ProtonDep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->NeutronDaughterDep_timesep[x_it][y_it][z_it] =
                NeutronDep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->PiCDaughterDep_timesep[x_it][y_it][z_it] =
                PiCDep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->Pi0DaughterDep_timesep[x_it][y_it][z_it] =
                Pi0Dep.DaughterDeposits_timesep->GetBinContent(gbin);
            fdw->OtherDaughterDep_timesep[x_it][y_it][z_it] =
                OtherDep.DaughterDeposits_timesep->GetBinContent(gbin);

#ifdef DEBUG
            if (fdw->ProtonDep[x_it][y_it][z_it] ||
                fdw->ProtonDep_timesep[x_it][y_it][z_it]) {
              std::cout << fdw->ProtonDep[x_it][y_it][z_it] << ", "
                        << fdw->ProtonDep_timesep[x_it][y_it][z_it]
                        << std::endl;
              if (fdw->ProtonDep[x_it][y_it][z_it] ==
                  fdw->ProtonDep_timesep[x_it][y_it][z_it]) {
                std::cout << "[INFO]: Found identical timesep deposits: "
                          << fdw->LepDep_timesep[x_it][y_it][z_it] << std::endl;
              }
            }
#endif
          }
        }
      }
    }

    fdTree->Fill();
    nfills++;

    if (loudevery && !(evnum % loudevery)) {
      std::cout << "[INFO]: Processed " << evnum << "/" << NMaxEvents
        << " entries." << std::endl;
    }

    evnum++;
  } while (g4ar.GetNextEvent() && (evnum < NMaxEvents));

  std::cout << "[INFO]: Filled the output tree " << nfills << " times."
            << std::endl;

  LepDep.DeleteDepositMaps();
  HadDep.DeleteDepositMaps();
  ProtonDep.DeleteDepositMaps();
  NeutronDep.DeleteDepositMaps();
  PiCDep.DeleteDepositMaps();
  Pi0Dep.DeleteDepositMaps();
  OtherDep.DeleteDepositMaps();
  NuclRemDep.DeleteDepositMaps();

  outfile->Write();
  outfile->Close();
}
