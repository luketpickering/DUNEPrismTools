#include "TFile.h"
#include "TChain.h"
#include "TH2D.h"

void MakeSlimFluxSanityPlots(char const *add) {
  TChain *dk2nuTree_lite = new TChain("dk2nuTree_lite");
  dk2nuTree_lite->Add(add);

  TFile *outf = new TFile("dk2nuTree_lite_SanityPlots.root","RECREATE");

  TH2D * p_xy = new TH2D("p_xy","",100,-15,15,100,-30,30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vx*1E-2 >> p_xy","","GOFF");
  TH2D * p_yz = new TH2D("p_yz","",100,0,300,100,-30,30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vz*1E-2 >> p_yz","","GOFF");

  TH2D * p_xy_imp = new TH2D("p_xy_imp","",100,-15,15,100,-30,30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vx*1E-2 >> p_xy_imp","decay_nimpwt","GOFF");
  TH2D * p_yz_imp = new TH2D("p_yz_imp","",100,0,300,100,-30,30);
  dk2nuTree_lite->Draw("decay_vy*1E-2:decay_vz*1E-2 >> p_yz_imp","decay_nimpwt","GOFF");


  outf->Write();
  outf->Close();
}

