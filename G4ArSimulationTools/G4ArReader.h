#include "DepoParticle.h"

#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"

#include <string>

struct DetectorAndFVDimensions {
  int NXSteps;
  double DetMin[3];
  double DetMax[3];
  double FVGap[3];
  int NYSteps;
  int NZSteps;
  DetectorAndFVDimensions()
      : NXSteps(400),
        DetMin{0, 0, 0},
        DetMax{0, 0, 0},
        FVGap{0, 0, 0},
        NYSteps(3),
        NZSteps(3) {}

  TH3D* BuildDetectorMap();
};

struct G4ArReader {
  TFile* InputG4ArFile;
  TTree* InputG4ArTree;

  TFile* InputGENIERooTrackerFile;
  TTree* InputGENIERooTrackerTree;

  std::vector<TH3D*> CacheDetectorMaps;

  double timesep_us;

  std::vector<int> TrackTimePDGs;
  void TrackTimeForPDG(int pdg) { TrackTimePDGs.push_back(pdg); }

  TH3D* GetCacheMap(size_t i);

  G4ArReader(std::string inputG4ArFileName, DetectorAndFVDimensions& detdims,
             std::string inputGENIERooTrackerFileName = "",
             double timesep_us = 0xdeadbeef, Long64_t MaxEntries = -1);

  bool SetBranchAddresses();

  bool GetNextEvent();
  bool GetEvent(Long64_t ev_it);

  void ResetCurrentEntry();

  Int_t NMaxTrackSteps;
  void SetNMaxTrackSteps(Int_t nmaxtrackssteps) {
    NMaxTrackSteps = nmaxtrackssteps;
  }

  static int const kMaxInit = 100;
  static int const kMaxTracks = 100000;
  static int const kMaxNQ = 10000000;

  Long64_t NInputEntries;
  Long64_t Entry;

  // ================== START Input branches ===================
  int ev_id;

  // Neutrino properties
  // Neutrino 4-momentum at vertex
  double nu_4mom[4];
  // Neutrino PDG
  int nu_PDG;
  double nu_vtx_pos[3];

  int n_prim;
  int prim_PDG[kMaxInit];
  double prim_3mom_x[kMaxInit];
  double prim_3mom_y[kMaxInit];
  double prim_3mom_z[kMaxInit];
  double prim_ekinetic[kMaxInit];
  double prim_mass[kMaxInit];

  int n_steps;

  int part_PDG[kMaxNQ];
  int part_ID[kMaxNQ];
  int part_parent_id[kMaxNQ];

  double part_step_pos_start_x[kMaxNQ];
  double part_step_pos_end_x[kMaxNQ];
  double part_step_pos_start_y[kMaxNQ];
  double part_step_pos_end_y[kMaxNQ];
  double part_step_pos_start_z[kMaxNQ];
  double part_step_pos_end_z[kMaxNQ];
  double part_step_time_end[kMaxNQ];
  double part_step_3mom_start_x[kMaxNQ];
  double part_step_3mom_end_x[kMaxNQ];
  double part_step_3mom_start_y[kMaxNQ];
  double part_step_3mom_end_y[kMaxNQ];
  double part_step_3mom_start_z[kMaxNQ];
  double part_step_3mom_end_z[kMaxNQ];
  double part_step_ekinetic[kMaxNQ];
  double part_step_edeposited[kMaxNQ];

  ~G4ArReader();

  TObjString* RooTrackerInteractionCode;
  Int_t StdHepN;
  Int_t StdHepStatus[1000];
  Int_t StdHepRescat[1000];
  Int_t StdHepPdg[1000];
  Double_t StdHepP4[1000][4];
  // ================== End Input branches ===================

  // ================== Start Util methods ===================

  static bool IsBindino(int pdg);
  static bool IsNuclearPDG(int pdg);

  void ShoutRooTracker();

  Event BuildEvent();
};
