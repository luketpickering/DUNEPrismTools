#include "FileSystemUtility.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// Expected FHiCL like
// Inputs: [ {
//   InputDirectory: "/path/to/files"
//   FilePattern: "file_pattern*.root"
// }, ]
// OutputFile: "Output_CAF.root"

class CAFReader {

  TFile *file;
  TTree *caf;
  TTree *meta;
  TH1D *RunPOT;

public:
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
  double Elep;
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

  double POTWeight;
  double FilePOT;

  CAFReader() : file(nullptr), caf(nullptr), meta(nullptr), RunPOT(nullptr) {}

  CAFReader(std::string const &filename) : CAFReader() {
    file = new TFile(filename.c_str(), "READ");
    if (!file || !file->IsOpen()) {
      std::cout << "[ERROR]: Failed to open TFile, " << filename << std::endl;
      throw;
    }

    caf = dynamic_cast<TTree *>(file->Get("caf"));
    if (!caf) {
      std::cout << "[ERROR]: Failed to open caf tree from TFile, " << filename
                << std::endl;
      return;
    }

    meta = dynamic_cast<TTree *>(file->Get("meta"));
    if (!meta) {
      std::cout << "[ERROR]: Failed to open meta tree from TFile " << filename
                << std::endl;
      throw;
    }

    meta->SetBranchAddress("pot", &FilePOT);
    meta->GetEntry(0);

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
    caf->SetBranchAddress("Elep", &Elep);
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

  size_t GetEntries() { return caf ? caf->GetEntries() : 0; }

  void GetEntry(size_t i) { caf->GetEntry(i); }

  CAFReader &operator=(CAFReader const &other) {

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
    Elep = other.Elep;
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

    POTWeight = 1.0 / other.FilePOT;

    return *this;
  }

  static CAFReader *MakeWriter(std::string const &filename) {

    CAFReader *wrt = new CAFReader();

    wrt->file = new TFile(filename.c_str(), "RECREATE");
    if (!wrt->file || !wrt->file->IsOpen()) {
      std::cout << "[ERROR]: Failed to open TFile, " << filename << std::endl;
      throw;
    }

    wrt->caf = new TTree("caf","");

    wrt->caf->Branch("Ev_reco", &wrt->Ev_reco, "Ev_reco/D");
    wrt->caf->Branch("Elep_reco", &wrt->Elep_reco, "Elep_reco/D");
    wrt->caf->Branch("theta_reco", &wrt->theta_reco, "theta_reco/D");
    wrt->caf->Branch("Ehad_veto", &wrt->Ehad_veto, "Ehad_veto/D");
    wrt->caf->Branch("reco_q", &wrt->reco_q, "reco_q/D");
    wrt->caf->Branch("reco_numu", &wrt->reco_numu, "reco_numu/D");
    wrt->caf->Branch("reco_nue", &wrt->reco_nue, "reco_nue/D");
    wrt->caf->Branch("reco_nc", &wrt->reco_nc, "reco_nc/D");
    wrt->caf->Branch("muon_contained", &wrt->muon_contained,
                     "muon_contained/D");
    wrt->caf->Branch("muon_tracker", &wrt->muon_tracker, "muon_tracker/D");
    wrt->caf->Branch("muon_ecal", &wrt->muon_ecal, "muon_ecal/D");
    wrt->caf->Branch("muon_exit", &wrt->muon_exit, "muon_exit/D");
    wrt->caf->Branch("Ev", &wrt->Ev, "Ev/D");
    wrt->caf->Branch("Elep", &wrt->Elep, "Elep/D");
    wrt->caf->Branch("isCC", &wrt->isCC, "isCC/D");
    wrt->caf->Branch("nuPDG", &wrt->nuPDG, "nuPDG/D");
    wrt->caf->Branch("LepPDG", &wrt->LepPDG, "LepPDG/D");
    wrt->caf->Branch("mode", &wrt->mode, "mode/D");
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
    wrt->caf->Branch("run", &wrt->run, "run/D");
    wrt->caf->Branch("isFD", &wrt->isFD, "isFD/D");
    wrt->caf->Branch("isFHC", &wrt->isFHC, "isFHC/D");
    wrt->caf->Branch("POTWeight", &wrt->POTWeight, "POTWeight/D");

    wrt->RunPOT = new TH1D("RunPOT", "RunPOT;Position (m);POT", 4500, -5, 40);
    wrt->RunPOT->SetDirectory(wrt->file);

    return wrt;
  }

  void NewFile() {
    if (!RunPOT) {
      std::cout << "[ERROR]: NewFile called on non-writer instance"
                << std::endl;
    } else {
      for (size_t i = 0; i < 700; ++i) {
        double det_pos = (double(i) - 349.5) / 100.0;
        RunPOT->Fill(det_pos + det_x, 1.0 / POTWeight);
      }
    }
  }
  void Fill() { caf->Fill(); }

  ~CAFReader() {
    if (file) {
      file->Close();
    }
  }
};

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string oupf = ps.get<std::string>("OutputFile");

  CAFReader *wrtr = CAFReader::MakeWriter(oupf);

  for (fhicl::ParameterSet input :
       ps.get<std::vector<fhicl::ParameterSet>>("Inputs")) {
    std::string idir = input.get<std::string>("InputDirectory");
    std::string pat = input.get<std::string>("FilePattern");

    std::vector<std::string> MatchingFileList = GetMatchingFiles(idir, pat);

    for (std::string const &inpf : MatchingFileList) {
      CAFReader rdr(inpf);
      std::cout << "[INFO]: Reading " << rdr.GetEntries()
                << " CAF events from file " << inpf << std::endl;

      size_t fents = rdr.GetEntries();
      for(size_t ent = 0; ent < fents; ++ent){
        rdr.GetEntry(ent);

        (*wrtr) = rdr;

        if(!ent){
          wrtr->NewFile();
        }

        wrtr->Fill();
      }
    }
  }
}
