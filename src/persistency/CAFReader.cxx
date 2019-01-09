#include "CAFReader.hxx"
#include "DUNETDRNDHelper.hxx"

#include "ROOTUtility.hxx"

#include <iostream>

CAFReader::CAFReader(std::string const &filename) : CAFReader() {
  file = CheckOpenFile(filename.c_str(), "READ");

  caf = dynamic_cast<TTree *>(file->Get("caf"));
  if (!caf) {
    std::cout << "[ERROR]: Failed to open caf tree from TFile, " << filename
              << std::endl;
    return;
  }

  // If we are loading from a precombined file, the exposure of the run should
  // be in the file and can be used to build event weights.
  TH1D *RunPOTHist = GetHistogram<TH1D>(file, "RunPOT", false);
  HasRunPOTWeight = bool(RunPOTHist);

  if (!HasRunPOTWeight) {
    meta = dynamic_cast<TTree *>(file->Get("meta"));
    if (!meta) {
      std::cout << "[ERROR]: Failed to open meta tree from TFile " << filename
                << std::endl;
      throw;
    }

    double pot;
    meta->SetBranchAddress("pot", &pot);
    size_t nmeta_ents = meta->GetEntries();
    for (size_t meta_it = 0; meta_it < nmeta_ents; ++meta_it) {
      meta->GetEntry(meta_it);
      FilePOT += pot;
    }
    meta->ResetBranchAddresses();
  } else {
    RunPOT = RunPOTHist;
  }

  // Set caf branch addresses
  caf->SetBranchAddress("Ev_reco", &Ev_reco);
  caf->SetBranchAddress("Elep_reco", &Elep_reco);
  caf->SetBranchAddress("theta_reco", &theta_reco);
  caf->SetBranchAddress("Ehad_veto", &Ehad_veto);
  caf->SetBranchAddress("reco_q", &reco_q);
  caf->SetBranchAddress("reco_numu", &reco_numu);
  caf->SetBranchAddress("reco_nue", &reco_nue);
  caf->SetBranchAddress("reco_nc", &reco_nc);
  caf->SetBranchAddress("muon_contained", &muon_contained);
  caf->SetBranchAddress("muon_tracker", &muon_tracker);
  caf->SetBranchAddress("muon_ecal", &muon_ecal);
  caf->SetBranchAddress("muon_exit", &muon_exit);
  caf->SetBranchAddress("Ev", &Ev);
  caf->SetBranchAddress("isCC", &isCC);
  caf->SetBranchAddress("nuPDG", &nuPDG);
  caf->SetBranchAddress("LepPDG", &LepPDG);
  caf->SetBranchAddress("mode", &mode);
  caf->SetBranchAddress("Q2", &Q2);
  caf->SetBranchAddress("W", &W);
  caf->SetBranchAddress("Y", &Y);
  caf->SetBranchAddress("X", &X);
  caf->SetBranchAddress("det_x", &det_x);
  caf->SetBranchAddress("vtx_x", &vtx_x);
  caf->SetBranchAddress("vtx_y", &vtx_y);
  caf->SetBranchAddress("vtx_z", &vtx_z);
  caf->SetBranchAddress("eP", &eP);
  caf->SetBranchAddress("eN", &eN);
  caf->SetBranchAddress("ePip", &ePip);
  caf->SetBranchAddress("ePim", &ePim);
  caf->SetBranchAddress("ePi0", &ePi0);
  caf->SetBranchAddress("eOther", &eOther);
  caf->SetBranchAddress("NuMomX", &NuMomX);
  caf->SetBranchAddress("NuMomY", &NuMomY);
  caf->SetBranchAddress("NuMomZ", &NuMomZ);
  caf->SetBranchAddress("LepMomX", &LepMomX);
  caf->SetBranchAddress("LepMomY", &LepMomY);
  caf->SetBranchAddress("LepMomZ", &LepMomZ);
  caf->SetBranchAddress("LepE", &LepE);
  caf->SetBranchAddress("LepNuAngle", &LepNuAngle);
  caf->SetBranchAddress("run", &run);
  caf->SetBranchAddress("isFD", &isFD);
  caf->SetBranchAddress("isFHC", &isFHC);
}

size_t CAFReader::GetEntries() { return caf ? caf->GetEntries() : 0; }

void CAFReader::GetEntry(size_t i) {
  caf->GetEntry(i);
  double rpot = 0;
  if (HasRunPOTWeight) {
    double xpos = vtx_x * 1E-2 + det_x;
    rpot = RunPOT->GetBinContent(RunPOT->GetXaxis()->FindFixBin(xpos));
  }
  if (std::isnormal(rpot)) {
    POTWeight = 1.0 / rpot;
  } else {
    POTWeight = 0;
  }
}

CAFReader &CAFReader::operator=(CAFReader const &other) {

  Ev_reco = other.Ev_reco;
  Elep_reco = other.Elep_reco;
  theta_reco = other.theta_reco;
  Ehad_veto = other.Ehad_veto;
  reco_q = other.reco_q;
  reco_numu = other.reco_numu;
  reco_nue = other.reco_nue;
  reco_nc = other.reco_nc;
  muon_contained = other.muon_contained;
  muon_tracker = other.muon_tracker;
  muon_ecal = other.muon_ecal;
  muon_exit = other.muon_exit;
  Ev = other.Ev;
  isCC = other.isCC;
  nuPDG = other.nuPDG;
  LepPDG = other.LepPDG;
  mode = other.mode;
  Q2 = other.Q2;
  W = other.W;
  Y = other.Y;
  X = other.X;
  det_x = other.det_x;
  vtx_x = other.vtx_x;
  vtx_y = other.vtx_y;
  vtx_z = other.vtx_z;
  eP = other.eP;
  eN = other.eN;
  ePip = other.ePip;
  ePim = other.ePim;
  ePi0 = other.ePi0;
  eOther = other.eOther;
  NuMomX = other.NuMomX;
  NuMomY = other.NuMomY;
  NuMomZ = other.NuMomZ;
  LepMomX = other.LepMomX;
  LepMomY = other.LepMomY;
  LepMomZ = other.LepMomZ;
  LepE = other.LepE;
  LepNuAngle = other.LepNuAngle;
  run = other.run;
  isFD = other.isFD;
  isFHC = other.isFHC;

  FilePOT = other.FilePOT;

  return *this;
}

CAFReader *CAFReader::MakeWriter(std::string const &filename) {

  CAFReader *wrt = new CAFReader();

  wrt->file = CheckOpenFile(filename.c_str(), "RECREATE");

  wrt->caf = new TTree("caf", "");

  wrt->caf->Branch("Ev_reco", &wrt->Ev_reco, "Ev_reco/D");
  wrt->caf->Branch("Elep_reco", &wrt->Elep_reco, "Elep_reco/D");
  wrt->caf->Branch("theta_reco", &wrt->theta_reco, "theta_reco/D");
  wrt->caf->Branch("Ehad_veto", &wrt->Ehad_veto, "Ehad_veto/D");
  wrt->caf->Branch("reco_q", &wrt->reco_q, "reco_q/I");
  wrt->caf->Branch("reco_numu", &wrt->reco_numu, "reco_numu/I");
  wrt->caf->Branch("reco_nue", &wrt->reco_nue, "reco_nue/I");
  wrt->caf->Branch("reco_nc", &wrt->reco_nc, "reco_nc/I");
  wrt->caf->Branch("muon_contained", &wrt->muon_contained, "muon_contained/I");
  wrt->caf->Branch("muon_tracker", &wrt->muon_tracker, "muon_tracker/I");
  wrt->caf->Branch("muon_ecal", &wrt->muon_ecal, "muon_ecal/I");
  wrt->caf->Branch("muon_exit", &wrt->muon_exit, "muon_exit/I");
  wrt->caf->Branch("Ev", &wrt->Ev, "Ev/D");
  wrt->caf->Branch("isCC", &wrt->isCC, "isCC/I");
  wrt->caf->Branch("nuPDG", &wrt->nuPDG, "nuPDG/I");
  wrt->caf->Branch("LepPDG", &wrt->LepPDG, "LepPDG/I");
  wrt->caf->Branch("mode", &wrt->mode, "mode/I");
  wrt->caf->Branch("Q2", &wrt->Q2, "Q2/D");
  wrt->caf->Branch("W", &wrt->W, "W/D");
  wrt->caf->Branch("Y", &wrt->Y, "Y/D");
  wrt->caf->Branch("X", &wrt->X, "X/D");
  wrt->caf->Branch("det_x", &wrt->det_x, "det_x/D");
  wrt->caf->Branch("vtx_x", &wrt->vtx_x, "vtx_x/D");
  wrt->caf->Branch("vtx_y", &wrt->vtx_y, "vtx_y/D");
  wrt->caf->Branch("vtx_z", &wrt->vtx_z, "vtx_z/D");
  wrt->caf->Branch("eP", &wrt->eP, "eP/D");
  wrt->caf->Branch("eN", &wrt->eN, "eN/D");
  wrt->caf->Branch("ePip", &wrt->ePip, "ePip/D");
  wrt->caf->Branch("ePim", &wrt->ePim, "ePim/D");
  wrt->caf->Branch("ePi0", &wrt->ePi0, "ePi0/D");
  wrt->caf->Branch("eOther", &wrt->eOther, "eOther/D");
  wrt->caf->Branch("NuMomX", &wrt->NuMomX, "NuMomX/D");
  wrt->caf->Branch("NuMomY", &wrt->NuMomY, "NuMomY/D");
  wrt->caf->Branch("NuMomZ", &wrt->NuMomZ, "NuMomZ/D");
  wrt->caf->Branch("LepMomX", &wrt->LepMomX, "LepMomX/D");
  wrt->caf->Branch("LepMomY", &wrt->LepMomY, "LepMomY/D");
  wrt->caf->Branch("LepMomZ", &wrt->LepMomZ, "LepMomZ/D");
  wrt->caf->Branch("LepE", &wrt->LepE, "LepE/D");
  wrt->caf->Branch("LepNuAngle", &wrt->LepNuAngle, "LepNuAngle/D");
  wrt->caf->Branch("run", &wrt->run, "run/I");
  wrt->caf->Branch("isFD", &wrt->isFD, "isFD/I");
  wrt->caf->Branch("isFHC", &wrt->isFHC, "isFHC/I");

  wrt->RunPOT = new TH1D("RunPOT", "RunPOT;Position (m);POT", 4500, -5, 40);
  wrt->RunPOT->SetDirectory(wrt->file);
  wrt->StopFiles =
      new TH1D("StopFiles", "StopFiles;Position (m);NFiles", 4500, -5, 40);
  wrt->StopFiles->SetDirectory(wrt->file);

  return wrt;
}

void CAFReader::NewFile(std::function<bool(double const &)> const &IsSel) {
  if (!RunPOT) {
    std::cout << "[ERROR]: NewFile called on non-writer instance" << std::endl;
  } else {
    for (size_t i = 0; i < 700; ++i) {
      if (!IsSel(double(i) - 349.5)) {
        continue;
      }
      double det_pos = (double(i) - 349.5) / 100.0;
      RunPOT->Fill(det_pos + det_x, FilePOT);
      StopFiles->Fill(det_pos + det_x);
    }
    nfiles++;
  }
}
size_t CAFReader::GetNFiles() { return nfiles; }
void CAFReader::Fill() { caf->Fill(); }

CAFReader::~CAFReader() {
  if (file) {
    std::cout << "[INFO]: Closing file " << file->GetName() << std::endl;
    file->Close();
  }
}
