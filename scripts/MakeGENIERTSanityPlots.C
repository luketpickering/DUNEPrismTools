void MakeGENIERTSanityPlots(char const *add) {
  TChain *gRooTracker = TChain("gRooTracker");
  gRooTracker->Add(add);

  TFile *outf = new TFile("GENIE_SanityPlots.root","RECREATE")

  TH2D * p_xy = new TH2D("p_xy","",100,-1.2,0.2,100,-10,0);
  gRooTracker->Draw("StdHepP4[0][1]:StdHepP4[0][0] >> p_xy");
  TH2D * p_yz = new TH2D("p_yz","",100,0,100,100,-10,0);
  gRooTracker->Draw("StdHepP4[0][1]:StdHepP4[0][2] >> p_yz");
  TH2D * evt_xy = new TH2D("evt_xy","",180,-40,5,64,-2,2);
  gRooTracker->Draw("EvtVtx[1]:EvtVtx[0] >> evt_xy");
  TH2D * evt_yz = new TH2D("evt_yz","",100,-3,3,100,-3,3);
  gRooTracker->Draw("EvtVtx[1]:EvtVtx[2] >> evt_yz");

}

