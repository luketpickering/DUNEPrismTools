#include "G4ArReader.h"

#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"

#include <csignal>
#include <iostream>
#include <string>

// #define DEBUG

TH3D* DetectorAndFVDimensions::BuildDetectorMap() {
  std::vector<double> YBins = {DetMin[1], DetMin[1] + FVGap[1],
                               DetMax[1] - FVGap[1], DetMax[1]};
  std::vector<double> ZBins = {DetMin[2], DetMin[2] + FVGap[2],
                               DetMax[2] - FVGap[2], DetMax[2]};

  std::vector<double> XBins;
  double step =
      ((DetMax[0] - FVGap[0]) - (DetMin[0] + FVGap[0])) / double(NXSteps - 2);
  XBins.push_back(DetMin[0]);
  XBins.push_back(DetMin[0] + FVGap[0]);

#ifdef DEBUG
  std::cout << "Building detector map: X { " << XBins.back() << ", "
            << std::flush;
#endif

  for (int i = 0; i < (NXSteps-2); ++i) {
    XBins.push_back(XBins.back() + step);
#ifdef DEBUG
    std::cout << XBins.back() << ((i != (NXSteps - 1)) ? ", " : "")
              << std::flush;
#endif
  }
  XBins.push_back(XBins.back() + FVGap[0]);
#ifdef DEBUG

  std::cout << "}" << std::endl;
  std::cout << "\t\t: Y {" << YBins[0] << ", " << YBins[1] << ", " << YBins[2]
            << ", " << YBins[3] << " }" << std::endl;
  std::cout << "\t\t: Z {" << ZBins[0] << ", " << ZBins[1] << ", " << ZBins[2]
            << ", " << ZBins[3] << " }" << std::endl;
#endif

  if (NYSteps == 0) {
    NYSteps = 3;
  }
  if (NZSteps == 0) {
    NZSteps = 3;
  }

  if (NYSteps != 3) {
    YBins.clear();
    double step_y = (DetMax[1] - DetMin[1]) / double(NYSteps);
    YBins.push_back(DetMin[1]);
    for (int i = 0; i < NYSteps; ++i) {
      YBins.push_back(YBins.back() + step_y);
    }
  }

  if (NZSteps != 3) {
    ZBins.clear();
    double step_z = (DetMax[2] - DetMin[2]) / double(NZSteps);
    ZBins.push_back(DetMin[2]);
    for (int i = 0; i < NYSteps; ++i) {
      ZBins.push_back(ZBins.back() + step_z);
    }
  }

  TH3D* Dm =
      new TH3D("dm", "", (XBins.size() - 1), XBins.data(), (YBins.size() - 1),
               YBins.data(), (ZBins.size() - 1), ZBins.data());
  Dm->SetDirectory(nullptr);
  Dm->SetName("dm_c");

  return Dm;
}

TH3D* G4ArReader::GetCacheMap(size_t i) {
  if (i >= CacheDetectorMaps.size()) {
    size_t os = CacheDetectorMaps.size();
    CacheDetectorMaps.resize(i + 1);
    for (size_t it = os; it < CacheDetectorMaps.size(); ++it) {
      CacheDetectorMaps[it] =
          static_cast<TH3D*>(CacheDetectorMaps.front()->Clone());
      CacheDetectorMaps[it]->SetDirectory(nullptr);
      CacheDetectorMaps[it]->Reset();
    }

    return GetCacheMap(i);
  }
  return CacheDetectorMaps[i];
}

G4ArReader::G4ArReader(std::string inputG4ArFileName,
                       DetectorAndFVDimensions& detdims,
                       std::string inputGENIERooTrackerFileName,
                       double timesep_us, Long64_t MaxEntries)
    : InputG4ArFile(nullptr),
      InputG4ArTree(nullptr),
      InputGENIERooTrackerFile(nullptr),
      InputGENIERooTrackerTree(nullptr),
      timesep_us(timesep_us),
      NMaxTrackSteps(1000),
      RooTrackerInteractionCode(nullptr) {
  InputG4ArFile = new TFile(inputG4ArFileName.c_str(), "READ");
  if (!InputG4ArFile || !InputG4ArFile->IsOpen()) {
    std::cout << "[ERROR]: Could not open " << inputG4ArFileName
              << " for reading." << std::endl;
    exit(1);
  }

  InputG4ArTree = dynamic_cast<TTree*>(InputG4ArFile->Get("argon"));
  if (!InputG4ArTree) {
    std::cout << "[ERROR]: Could not open \"argon\" tree from  "
              << inputG4ArFileName << "." << std::endl;
    exit(1);
  }

  if (inputGENIERooTrackerFileName.size()) {
    InputGENIERooTrackerFile =
        new TFile(inputGENIERooTrackerFileName.c_str(), "READ");
    if (!InputGENIERooTrackerFile || !InputGENIERooTrackerFile->IsOpen()) {
      std::cout << "[ERROR]: Could not open " << inputGENIERooTrackerFileName
                << " for reading." << std::endl;
      exit(1);
    }

    InputGENIERooTrackerTree =
        dynamic_cast<TTree*>(InputGENIERooTrackerFile->Get("gRooTracker"));
    if (!InputGENIERooTrackerTree) {
      std::cout << "[ERROR]: Could not open \"gRooTracker\" tree from  "
                << inputGENIERooTrackerFileName << "." << std::endl;
      exit(1);
    }
  }

  NInputEntries = (MaxEntries == -1)
                      ? InputG4ArTree->GetEntries()
                      : std::min(InputG4ArTree->GetEntries(), MaxEntries);

  std::cout << "[INFO]: Reading " << NInputEntries
            << " entries from input files: (" << inputG4ArFileName << ", "
            << inputGENIERooTrackerFileName << ")." << std::endl;

  CacheDetectorMaps.push_back(detdims.BuildDetectorMap());

  SetBranchAddresses();
  ResetCurrentEntry();
}

bool G4ArReader::SetBranchAddresses() {
  InputG4ArTree->SetBranchAddress("ev", &ev_id);

  InputG4ArTree->SetBranchAddress("ekina", &nu_4mom[3]);
  InputG4ArTree->SetBranchAddress("xa", &nu_vtx_pos[0]);
  InputG4ArTree->SetBranchAddress("ya", &nu_vtx_pos[1]);
  InputG4ArTree->SetBranchAddress("za", &nu_vtx_pos[2]);
  InputG4ArTree->SetBranchAddress("pxa", &nu_4mom[0]);
  InputG4ArTree->SetBranchAddress("pya", &nu_4mom[1]);
  InputG4ArTree->SetBranchAddress("pza", &nu_4mom[2]);
  InputG4ArTree->SetBranchAddress("pida", &nu_PDG);

  InputG4ArTree->SetBranchAddress("ni", &n_prim);
  InputG4ArTree->SetBranchAddress("pidi", &prim_PDG);
  InputG4ArTree->SetBranchAddress("pxi", &prim_3mom_x);
  InputG4ArTree->SetBranchAddress("pyi", &prim_3mom_y);
  InputG4ArTree->SetBranchAddress("pzi", &prim_3mom_z);
  InputG4ArTree->SetBranchAddress("ekini", &prim_ekinetic);
  InputG4ArTree->SetBranchAddress("mi", &prim_mass);

  InputG4ArTree->SetBranchAddress("nstep", &n_steps);
  InputG4ArTree->SetBranchAddress("pid", &part_PDG);
  InputG4ArTree->SetBranchAddress("tid", &part_ID);
  InputG4ArTree->SetBranchAddress("parid", &part_parent_id);
  InputG4ArTree->SetBranchAddress("xs", &part_step_pos_start_x);
  InputG4ArTree->SetBranchAddress("xe", &part_step_pos_end_x);
  InputG4ArTree->SetBranchAddress("ys", &part_step_pos_start_y);
  InputG4ArTree->SetBranchAddress("ye", &part_step_pos_end_y);
  InputG4ArTree->SetBranchAddress("zs", &part_step_pos_start_z);
  InputG4ArTree->SetBranchAddress("ze", &part_step_pos_end_z);
  InputG4ArTree->SetBranchAddress("te", &part_step_time_end);
  InputG4ArTree->SetBranchAddress("pxs", &part_step_3mom_start_x);
  InputG4ArTree->SetBranchAddress("pxe", &part_step_3mom_end_x);
  InputG4ArTree->SetBranchAddress("pys", &part_step_3mom_start_y);
  InputG4ArTree->SetBranchAddress("pye", &part_step_3mom_end_y);
  InputG4ArTree->SetBranchAddress("pzs", &part_step_3mom_start_z);
  InputG4ArTree->SetBranchAddress("pze", &part_step_3mom_end_z);
  InputG4ArTree->SetBranchAddress("ekin", &part_step_ekinetic);
  InputG4ArTree->SetBranchAddress("edep", &part_step_edeposited);

  if (InputGENIERooTrackerTree) {
    InputGENIERooTrackerTree->SetBranchAddress("EvtCode",
                                               &RooTrackerInteractionCode);
    InputGENIERooTrackerTree->SetBranchAddress("StdHepN", &StdHepN);
    InputGENIERooTrackerTree->SetBranchAddress("StdHepStatus", StdHepStatus);
    InputGENIERooTrackerTree->SetBranchAddress("StdHepRescat", StdHepRescat);
    InputGENIERooTrackerTree->SetBranchAddress("StdHepPdg", StdHepPdg);
    InputGENIERooTrackerTree->SetBranchAddress("StdHepP4", StdHepP4);
  }
  return true;
}

void G4ArReader::ShoutRooTracker() {
  if (InputGENIERooTrackerTree) {
    double e = 0;
    double nmass = 0;
    std::cout << "[INFO]: GENIE RooTracker info: " << std::endl;
    for (Int_t i = 0; i < StdHepN; ++i) {
      std::cout << "[" << i << "] Status: " << StdHepStatus[i]
                << ", PDG: " << StdHepPdg[i] << ", Energy = " << StdHepP4[i][3]
                << ", Rescat = " << StdHepRescat[i] << std::endl;
      if (StdHepStatus[i] == 1) {
        e += StdHepP4[i][3];
        nmass += ((abs(StdHepPdg[i]) == 2212) || (abs(StdHepPdg[i]) == 2112))
                     ? StdHepP4[i][3] - 0.938
                     : StdHepP4[i][3];
      }
    }
    std::cout << "Total E = " << e << ", no nucl mass = " << nmass << std::endl;
  }
}

bool G4ArReader::GetNextEvent() {
  if (Entry == (NInputEntries - 1)) {
    return false;
  }
  Entry++;
  InputG4ArTree->GetEntry(Entry);
  if (InputGENIERooTrackerTree) {
    InputGENIERooTrackerTree->GetEntry(Entry);
  }
  return true;
}

bool G4ArReader::GetEvent(Long64_t ev_it) {
  if (ev_it >= NInputEntries) {
    return false;
  }
  Entry = ev_it;
  InputG4ArTree->GetEntry(Entry);
  if (InputGENIERooTrackerTree) {
    InputGENIERooTrackerTree->GetEntry(Entry);
  }
  return true;
}

void G4ArReader::ResetCurrentEntry() {
  Entry = -1;
  GetNextEvent();
}

G4ArReader::~G4ArReader() {
  InputG4ArFile->Close();
  InputGENIERooTrackerFile->Close();

  delete InputG4ArFile;
  delete InputGENIERooTrackerFile;

  for (TH3D* map : CacheDetectorMaps) {
    delete map;
  }
}

bool G4ArReader::IsBindino(int pdg) { return (pdg == 2000000101); }
bool G4ArReader::IsNuclearPDG(int pdg) { return (pdg > 1000000000); }

Event G4ArReader::BuildEvent() {
  Event ev;

  ev.ev_id = ev_id;
  ev.RooTrackerInteractionCode = RooTrackerInteractionCode;

  ev.VertexPosition = TVector3(nu_vtx_pos[0], nu_vtx_pos[1], nu_vtx_pos[2]);

  PrimaryParticle neutrino;
  neutrino.IsFinalState = false;
  neutrino.PDG = nu_PDG;
  neutrino.EKin = nu_4mom[3];
  neutrino.ThreeMom = TVector3(nu_4mom[0], nu_4mom[1], nu_4mom[2]);

  ev.PrimaryParticles.push_back(neutrino);

  // Get primary particles and seed deposit association map
  size_t cachemap_it = 0;
  for (Long64_t prim_part_it = 0; prim_part_it < n_prim; ++prim_part_it) {
    if (IsBindino(prim_PDG[prim_part_it])) {
      continue;
    }

    size_t prim_TID = prim_part_it + 1;

    PrimaryParticle pp;
    pp.IsFinalState = true;
    pp.PDG = prim_PDG[prim_part_it];
    pp.EKin = prim_ekinetic[prim_part_it];
    pp.EMass = prim_mass[prim_part_it];
    pp.ThreeMom = TVector3(prim_3mom_x[prim_part_it], prim_3mom_y[prim_part_it],
                           prim_3mom_z[prim_part_it]);

    ev.PrimaryParticles.push_back(pp);

    if (abs(pp.PDG) == 13) {
      ev.TrackedDeposits.emplace_back(prim_PDG[prim_part_it], prim_TID,
                                      NMaxTrackSteps, timesep_us);
#ifdef DEBUG
      std::cout << "[INFO]: Added new primary particle: PDG = "
                << ev.TrackedDeposits.back().PDG
                << ", TID = " << ev.TrackedDeposits.back().TrackID << std::endl;
#endif
      ev.TrackedDeposits.back().LendDepositMaps(GetCacheMap(cachemap_it),
                                                GetCacheMap(cachemap_it + 1));
      cachemap_it = cachemap_it + 2;

      if (std::find(TrackTimePDGs.begin(), TrackTimePDGs.end(), pp.PDG) !=
          TrackTimePDGs.end()) {
        ev.TrackedDeposits.back().SetTrackTime();
      }
    } else {
      ev.TotalDeposits.emplace_back(prim_PDG[prim_part_it], prim_TID,
                                    timesep_us);
#ifdef DEBUG
      std::cout << "[INFO]: Added new primary particle: PDG = "
                << ev.TotalDeposits.back().PDG
                << ", TID = " << ev.TotalDeposits.back().TrackID << std::endl;
#endif
      ev.TotalDeposits.back().LendDepositMaps(GetCacheMap(cachemap_it),
                                              GetCacheMap(cachemap_it + 1));
      cachemap_it = cachemap_it + 2;

      if (std::find(TrackTimePDGs.begin(), TrackTimePDGs.end(), pp.PDG) !=
          TrackTimePDGs.end()) {
        ev.TotalDeposits.back().SetTrackTime();
      }
    }
  }

  // Loop through the simulation steps
  for (Int_t step_it = 0; step_it < n_steps; ++step_it) {
    double pos[4] = {part_step_pos_end_x[step_it], part_step_pos_end_y[step_it],
                     part_step_pos_end_z[step_it], part_step_time_end[step_it]};

    DepoParticle* deposit;
    if ((deposit = ev.GetPrimaryParticle(part_parent_id[step_it],
                                         part_ID[step_it]))) {
      deposit->AddDeposit(pos, part_step_edeposited[step_it],
                          ev.GetIsPrimary(part_ID[step_it]));
#ifdef DEBUG
      std::cout << "[INFO]: Added deposit of " << part_step_edeposited[step_it]
                << " GeV at {" << pos[0] << ", " << pos[1] << ", " << pos[2]
                << " }." << std::endl;
#endif
    } else {
      std::string evcode = ev.RooTrackerInteractionCode->GetString().Data();
      std::cout << "[WARN]: Ev code of missed particle = " << evcode
                << std::endl;
      std::cout << "[WARN]: ParticleID =  " << part_ID[step_it]
                << ", ParentID = " << part_parent_id[step_it]
                << ", PDG = " << part_PDG[step_it] << std::endl;
    }

    DepoTracked* trackeddeposit;
    if ((trackeddeposit = ev.GetTrackedDeposit(part_ID[step_it]))) {
      double mom[3] = {part_step_3mom_end_x[step_it],
                       part_step_3mom_end_y[step_it],
                       part_step_3mom_end_z[step_it]};
      trackeddeposit->AddStep(pos, mom);
#ifdef DEBUG
      std::cout << "[INFO]: Added step to tracked particle at {" << pos[0]
                << ", " << pos[1] << ", " << pos[2] << " } with 3-mom {"
                << mom[0] << ", " << mom[1] << ", " << mom[2] << " }."
                << std::endl;
#endif
    }
  }

  return ev;
}
