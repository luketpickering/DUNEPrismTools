#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"

void MakeSlimFluxSanityPlots(char const *add, bool append = false,
                             char const *subdir = "", bool maketree = true) {
  TChain *dk2nuTree_lite = new TChain("dk2nuTree_lite");
  dk2nuTree_lite->Add(add);

  TFile *oupF = new TFile("dk2nuTree_lite_SanityPlots.root",
                          append ? "UPDATE" : "RECREATE");

  TDirectory *oupD = oupF;
  std::string outputDir = subdir;
  if (outputDir.length()) {
    oupD = oupF->GetDirectory(outputDir.c_str());
    if (!oupD) {
      oupD = oupF->mkdir(outputDir.c_str());
    }
  }
  oupD->cd();

  TH2D *p_xy = new TH2D("p_xy", "", 100, -15, 15, 100, -30, 30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vx*1E-2 >> p_xy", "", "GOFF");
  TH2D *p_yz = new TH2D("p_yz", "", 100, 0, 300, 100, -30, 30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vz*1E-2 >> p_yz", "", "GOFF");

  TH2D *p_xy_imp = new TH2D("p_xy_imp", "", 100, -15, 15, 100, -30, 30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vx*1E-2 >> p_xy_imp",
                       "decay_nimpwt", "GOFF");
  TH2D *p_yz_imp = new TH2D("p_yz_imp", "", 100, 0, 300, 100, -30, 30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vz*1E-2 >> p_yz_imp",
                       "decay_nimpwt", "GOFF");

  TH2D *p_xy_imp_fine = new TH2D("p_xy_imp_fine", "", 40, -20, 20, 40, -20, 20);
  dk2nuTree_lite->Draw("decay_vy:decay_vx >> p_xy_imp_fine", "decay_nimpwt",
                       "GOFF");
  TH2D *p_yz_imp_fine =
      new TH2D("p_yz_imp_fine", "", 1000, 0, 10000, 60, -30, 30);
  dk2nuTree_lite->Draw("decay_vy:decay_vz >> p_yz_imp_fine", "decay_nimpwt",
                       "GOFF");

  if (maketree) {
    TTree *stree = new TTree("sum_tree", "");
    double parent_x, parent_y, parent_z, impweight;
    stree->Branch("parent_x", &parent_x);
    stree->Branch("parent_y", &parent_y);
    stree->Branch("parent_z", &parent_z);
    stree->Branch("impweight", &impweight);

    dk2nuTree_lite->SetBranchAddress("decay_vx", &parent_x);
    dk2nuTree_lite->SetBranchAddress("decay_vy", &parent_y);
    dk2nuTree_lite->SetBranchAddress("decay_vz", &parent_z);
    dk2nuTree_lite->SetBranchAddress("decay_nimpwt", &impweight);

    size_t npars = dk2nuTree_lite->GetEntries();
    for (size_t it = 0; it < npars: ++it) {
      dk2nuTree_lite->GetEntry(it);
      stree->Fill();
    }
    stree->SetDirectory(oupD);
  }

  p_xy->SetDirectory(oupD);
  p_yz->SetDirectory(oupD);
  p_xy_imp->SetDirectory(oupD);
  p_yz_imp->SetDirectory(oupD);
  p_xy_imp_fine->SetDirectory(oupD);
  p_yz_imp_fine->SetDirectory(oupD);

  oupF->Write();
  oupF->Close();
}
