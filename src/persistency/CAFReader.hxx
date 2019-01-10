#pragma once

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <functional>
#include <string>

class CAFReader {

public:
  TFile *file;
  TTree *caf;
  TTree *meta;
  TH1D *RunPOT;
  TH1D *StopFiles;
  size_t nfiles;

  // ND Reco info
  double Ev_reco;
  double Elep_reco;
  double theta_reco;
  double Ehad_veto;

  // ND Selection info
  int reco_q;
  int reco_numu;
  int reco_nue;
  int reco_nc;
  int muon_contained;
  int muon_tracker;
  int muon_ecal;
  int muon_exit;

  // FD Reco info
  double Ev_reco_nue;
  double Ev_reco_numu;

  // FD Selection info
  double cvnnumu;
  double cvnnue;

  // Truth info
  double Ev;
  int isCC;
  int nuPDG;
  int nuPDGunosc;
  int LepPDG;
  int mode;
  double Q2;
  double W;
  double Y;
  double X;

  // Near Detector off-axis position in meters.
  double det_x;

  // Interaction position in detector coordinates in cm;
  double vtx_x;
  double vtx_y;
  double vtx_z;

  // True energy of particles by species
  double eP;
  double eN;
  double ePip;
  double ePim;
  double ePi0;
  double eOther;

  double NuMomX;
  double NuMomY;
  double NuMomZ;
  double LepMomX;
  double LepMomY;
  double LepMomZ;
  double LepE;
  double LepNuAngle;

  // config
  int run;
  int isFD;
  int isFHC;

  double FilePOT;
  double POTWeight;

  bool HasRunPOTWeight;
  bool isNDFile;

  CAFReader()
      : file(nullptr), caf(nullptr), meta(nullptr), RunPOT(nullptr),
        StopFiles(nullptr), nfiles(0), HasRunPOTWeight(false), isNDFile(false) {
    Ev_reco = 0;
    Elep_reco = 0;
    theta_reco = 0;
    Ehad_veto = 0;
    reco_q = 0;
    reco_numu = 0;
    reco_nue = 0;
    reco_nc = 0;
    muon_contained = 0;
    muon_tracker = 0;
    muon_ecal = 0;
    muon_exit = 0;
    Ev_reco_nue = 0;
    Ev_reco_numu = 0;
    cvnnumu = 0;
    cvnnue = 0;
    Ev = 0;
    isCC = 0;
    nuPDG = 0;
    nuPDGunosc = 0;
    LepPDG = 0;
    mode = 0;
    Q2 = 0;
    W = 0;
    Y = 0;
    X = 0;
    det_x = 0;
    vtx_x = 0;
    vtx_y = 0;
    vtx_z = 0;
    eP = 0;
    eN = 0;
    ePip = 0;
    ePim = 0;
    ePi0 = 0;
    eOther = 0;
    NuMomX = 0;
    NuMomY = 0;
    NuMomZ = 0;
    LepMomX = 0;
    LepMomY = 0;
    LepMomZ = 0;
    LepE = 0;
    LepNuAngle = 0;
    run = 0;
    isFD = 0;
    isFHC = 0;
  }

  CAFReader(std::string const &filename);

  size_t GetEntries();
  void GetEntry(size_t i);

  CAFReader &operator=(CAFReader const &other);

  static CAFReader *MakeWriter(std::string const &filename, bool);

  void NewFile(std::function<bool(double const &)> const &IsSel);
  size_t GetNFiles();
  void Fill();

  ~CAFReader();
};
