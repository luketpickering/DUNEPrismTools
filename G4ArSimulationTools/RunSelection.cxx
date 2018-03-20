#include "EDepTreeReader.h"
#include "TTree.h"
#include "Utils.hxx"

#include <string>
#include <vector>

std::string inpfile;
std::string oupfile;

double HadrVeto = 0;  // GeV
bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>    : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV.\n"
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
      SelMuExit = true;
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

  TChain *config_in = new TChain("configTree", "");
  config_in->Add(inpfile.c_str());
  Int_t NStops;
  config_in->SetBranchAddress("NStops", &NStops);
  double timesep_us = 0xdeadbeef;
  if (CheckTTreeHasBranch(config_in, "timesep_us")) {
    config_in->SetBranchAddress("timesep_us", &timesep_us);
  }
  config_in->GetEntry(0);

  TChain *stopConfig_in = new TChain("stopConfigTree", "");
  stopConfig_in->Add(inpfile.c_str());

  EDep edr("EDeps", inpfile);

  TFile *of = CheckOpenFile(oupfile, "RECREATE");

  TTree *config_copy = config_in->CloneTree();
  TTree *stopConfig_copy = stopConfig_in->CloneTree();
  config_copy->SetDirectory(of);
  stopConfig_copy->SetDirectory(of);

  TTree *OutputTree = new TTree("EDeps", "");

  TH1D *EventRates_Selected =
      new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
  TH1D *EventRates_True =
      new TH1D("EventRates_True", ";Offset (cm);Count", 400, -250, 3750);

  TH1D *ERec_Selected = new TH1D(
      "ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count", 150, 0, 15);
  TH1D *EHadr_Selected =
      new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count", 150, 0, 15);
  TH1D *EHadr_True =
      new TH1D("EHadr_True", ";E_{Hadr} (GeV);Count", 150, 0, 15);

  TH2D *EvRateERec_Selected = new TH2D(
      "EvRateERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
      100, 0, 10, 400, -250, 3750);

  TH2D *EvRateEHadr_Selected =
      new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);
  TH2D *EvRateEHadr_True =
      new TH2D("EvRateEHadr_True", ";E_{Hadr} (GeV);Offset (cm)", 100, 0, 10,
               400, -250, 3750);

  TH2D *EvRateEVeto = new TH2D("EvRateEVeto", ";E_{Hadr} (GeV);Offset (cm)",
                               100, 0, 2, 400, -250, 3750);

  TDirectory *oupD = of;

  TDirectory *wD = oupD->mkdir("StopEventRates");
  wD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEventRates_True;
  std::vector<TH1D *> StopEventRates_TrueUnWeighted;

  for (size_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_True.push_back(new TH1D(
        (std::string("StopEventRates_True_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_TrueUnWeighted.push_back(new TH1D(
        (std::string("StopEventRates_TrueUnWeighted_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
  }
  oupD->cd();

  EDep *OutputEDep = EDep::MakeTreeWriter(OutputTree, timesep_us, true);

  size_t loud_every = edr.GetEntries() / 10;

  Long64_t NEntries = edr.GetEntries();

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    if ((edr.stop < 0) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    // Fill truth histos
    EventRates_True->Fill(-1 * edr.vtx[0], edr.stop_weight);
    EHadr_True->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                     edr.stop_weight);
    EvRateEHadr_True->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                           -1 * edr.vtx[0], edr.stop_weight);
    StopEventRates_True[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);
    StopEventRates_TrueUnWeighted[edr.stop]->Fill(-1 * edr.vtx[0]);
    EvRateEVeto->Fill(-1 * edr.vtx[0], edr.TotalNonlep_Dep_veto);
    if ((edr.stop < 0) || (SelMuExitKE && (edr.LepExitKE < SelMuExitKE)) ||
        (edr.TotalNonlep_Dep_veto > HadrVeto) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    // Fill selected histos
    EventRates_Selected->Fill(-1 * edr.vtx[0], edr.stop_weight);
    ERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                            edr.TotalNonlep_Dep_veto,
                        edr.stop_weight);
    EHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                         edr.stop_weight);

    EvRateERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                  edr.TotalNonlep_Dep_veto,
                              -1 * edr.vtx[0], edr.stop_weight);
    EvRateEHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                               -1 * edr.vtx[0], edr.stop_weight);
    StopEventRates_Selected[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);

    OutputEDep->Copy(edr);
    OutputTree->Fill();
  }
  of->Write();
}
