#pragma once

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <string>

class CAFReader {

public:
  TFile *file;
  TTree *caf;
  TTree *meta;
  TH1D *RunPOT;
  TH1D *StopFiles;
  size_t nfiles;

  // Reco info
  double Ev_reco;
  double Elep_reco;
  double theta_reco;
  double Ehad_veto;

  // Selection info
  int reco_q;
  int reco_numu;
  int reco_nue;
  int reco_nc;
  int muon_contained;
  int muon_tracker;
  int muon_ecal;
  int muon_exit;

  // Truth info
  double Ev;
  int isCC;
  int nuPDG;
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

  CAFReader()
      : file(nullptr), caf(nullptr), meta(nullptr), RunPOT(nullptr),
        StopFiles(nullptr), nfiles(0), HasRunPOTWeight(false) {}

  CAFReader(std::string const &filename);

  size_t GetEntries();
  void GetEntry(size_t i);

  CAFReader &operator=(CAFReader const &other);

  static CAFReader *MakeWriter(std::string const &filename);

  void NewFile();
  size_t GetNFiles();
  void Fill();

  ~CAFReader();
};
