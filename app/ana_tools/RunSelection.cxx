#include "DepositsSummaryTreeReader.hxx"
#include "SimConfigTreeReader.hxx"
#include "StopConfigTreeReader.hxx"
#include "SelectionSummaryTreeReader.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "BoundingBox.hxx"

#include "TTree.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

std::string InputDepostisSummaryFile;
std::string OutputFile;
bool BuildMissingProtonEFakeData = false;

double HadrVeto = 0;  // GeV
bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV
double HadrVisAcceptance = 0; //GeV
std::vector<double> ExtraVertexSelectionPadding;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <stopprocessor.root>     : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
         "\t-FV <fvx,y,z>               : Vertex selection fiducial volume "
         "padding inside of the \n"
         "\t                              non-veto active region "
         "{Default: 0,0,0}.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV.\n"
         "\t-A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower"
         " energy. Showers with\n"
         "\t                              more than this energy are cut and "
         "filled in by far detector MC.\n"
         "\t-FDproton                   : Build missing proton energy fake "
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
      InputDepostisSummaryFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
    } else if (std::string(argv[opt]) == "-A") {
      HadrVisAcceptance = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-FDproton") {
      BuildMissingProtonEFakeData = true;
    }
    else if (std::string(argv[opt]) == "-FV") {
      ExtraVertexSelectionPadding = ParseToVect<double>(argv[++opt], ",");
      if(ExtraVertexSelectionPadding.size() != 3){
        std::cout << "[ERROR]: -FV option contained " <<
          ExtraVertexSelectionPadding.size() << " entries, expected 3."
          << std::endl;
        SayUsage(argv);
        exit(1);
      }
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

  if(!ExtraVertexSelectionPadding.size()){
      int argc_dum = 3;
      char const * argv_dum[] = { argv[0], "-FV", "0,0,0"};
      handleOpts(argc_dum, argv_dum);
  }

  SimConfig simCRdr(InputDepostisSummaryFile);
  StopConfig csRdr(InputDepostisSummaryFile);

  std::vector<BoundingBox> FVs = csRdr.GetStopBoundingBoxes(true,
    {ExtraVertexSelectionPadding[0],
      ExtraVertexSelectionPadding[1], ExtraVertexSelectionPadding[2]});

  DepositsSummary DepSumRdr(InputDepostisSummaryFile);

  Double_t ProtonFakeDataWeight = 1;
  if (BuildMissingProtonEFakeData) {
    // Should probably check that this exists, but leaving hardcoded and
    // presumptive for now as CheckTTreeHasBranch doesn't seem to work on
    // pre-friended trees.
    DepSumRdr.tree->SetBranchAddress("EnuTp", &ProtonFakeDataWeight);
  }

  TFile *of = CheckOpenFile(OutputFile, "RECREATE");
  of->cd();

  SimConfig *simCWtr = SimConfig::MakeTreeWriter();
  StopConfig *csWtr = StopConfig::MakeTreeWriter();

  simCRdr.GetEntry(0);
  simCWtr->Copy(simCRdr);
  simCWtr->Fill();

  csRdr.DetermineNStops();
  for(Long64_t stop_it = 0; stop_it < csRdr.NStops; ++stop_it){
    csRdr.GetEntry(stop_it);
    csWtr->Copy(csRdr);
    csWtr->Fill();
  }

  //Selection summary tree
  SelectionSummary *ss = SelectionSummary::MakeTreeWriter();

  ss->SelectOnMuonExit = SelMuExit;
  ss->MuonExitKECut_MeV = SelMuExitKE*1E3;
  ss->HadronicShowerVetoCut_MeV = HadrVeto*1E3;
  std::copy_n(ExtraVertexSelectionPadding.begin(),3,ss->VertexSelectionFV);

  //Selected output tree
  DepositsSummary *OutputDepSum =
    DepositsSummary::MakeTreeWriter(simCRdr.timesep_us, true);

  size_t loud_every = DepSumRdr.GetEntries() / 10;
  Long64_t NEntries = DepSumRdr.GetEntries();

  for(Long64_t f_it = 0; f_it < simCRdr.GetEntries(); ++f_it){
    simCRdr.GetEntry(f_it);
    ss->TotalPOT += simCRdr.POTPerFile;
  }
  std::cout << "[INFO]: Total POT = " << ss->TotalPOT << std::endl;

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    DepSumRdr.GetEntry(e_it);

    ss->NTotal++;
    ss->NInStops += (DepSumRdr.stop > -1);
    ss->NNumuCC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == 13);
    ss->NNumuNC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == 14);
    ss->NNumubCC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == -13);
    ss->NNumubNC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == -14);
    ss->NNueCC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == 11);
    ss->NNueNC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == 12);
    ss->NNuebCC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == -11);
    ss->NNuebNC += (DepSumRdr.stop > -1) && (DepSumRdr.PrimaryLepPDG == -12);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << DepSumRdr.vtx[0] << ", " << DepSumRdr.vtx[1] << ", " << DepSumRdr.vtx[2]
                << "}, Enu: " << DepSumRdr.nu_4mom[3] << " )" << std::endl;
      if (BuildMissingProtonEFakeData) {
        std::cout << "\tProton fake data weight = " << ProtonFakeDataWeight
                  << std::endl;
      }
    }

    ///Add smearing here and then select on smeared quantities
    /// e.g. PID/charge mis-identification
    /// Muon momentum reconstruction.


    bool NumuInStop = ((DepSumRdr.stop >= 0) && (DepSumRdr.PrimaryLepPDG == 13));
    if(!NumuInStop){
      continue;
    }

    bool SelMu = ((SelMuExitKE == 0) || (DepSumRdr.LepExitKE > SelMuExitKE));
    bool SelHadr = (DepSumRdr.TotalNonlep_Dep_veto < HadrVeto);
    if(!(SelMu && SelHadr)){
      continue;
    }

    bool InFV = FVs[DepSumRdr.stop].Contains(
      {DepSumRdr.vtx[0], DepSumRdr.vtx[1], DepSumRdr.vtx[2]});
    if(!InFV){
      ss->NOOFV++;
      continue;
    }

    bool HadrAccepted = ((HadrVisAcceptance == 0) ||
      (DepSumRdr.GetProjection(DepositsSummary::kEHadr_vis) <= HadrVisAcceptance));
    if(!HadrAccepted){
      ss->NOOAcceptance++;
      continue;
    }

    if (BuildMissingProtonEFakeData) {
      DepSumRdr.stop_weight *= ProtonFakeDataWeight;
      DepSumRdr.TotalNonlep_Dep_veto -= DepSumRdr.ProtonDep_veto * 0.2;
      DepSumRdr.TotalNonlep_Dep_FV -= DepSumRdr.ProtonDep_FV * 0.2;
      DepSumRdr.ERecProxy_True -= DepSumRdr.EKinProton_True * 0.2;
    }

    ss->NSel++;

    OutputDepSum->Copy(DepSumRdr);
    OutputDepSum->Fill();
  }
  ss->Fill();

  std::cout << "[INFO]: Selection summary: " << std::endl;
  ss->tree->Show(0);

  of->Write();
}
