void VALORXSecPlots(char const* inps, char const* wfile) {
  TH1D::SetDefaultSumw2(true);

  TChain* EDeps = new TChain("EDeps");
  EDeps->Add(inps);
  EDeps->AddFriend("XSecWeights", wfile);

  TCanvas* c1 = new TCanvas("c1", "");
  c1->Divide(1, 2);
  c1->cd(1);
  c1->SaveAs("VALORXSecPlots.pdf[");

  TH1D* q2_nom_QE = new TH1D("q2_nom_QE", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_nom_QE",
              "(stop >= 0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==1)",
              "GOFF");

  TH1D* q2_p1_QE = new TH1D("q2_p1_QE", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_p1_QE",
              "XSecWeights[0]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==1))",
              "GOFF");

  TH1D* q2_m1_QE = new TH1D("q2_m1_QE", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_m1_QE",
              "XSecWeights[1]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==1))",
              "GOFF");

  q2_p1_QE->SetLineStyle(2);
  q2_p1_QE->SetLineColor(kBlue);
  q2_m1_QE->SetLineStyle(2);
  q2_m1_QE->SetLineColor(kRed);

  q2_p1_QE->Draw("EHIST");
  q2_nom_QE->Draw("EHISTSAME");
  q2_m1_QE->Draw("EHISTSAME");

  c1->cd(2);

  TH1D* q2_m1_QE_rat = q2_m1_QE->Clone();
  TH1D* q2_p1_QE_rat = q2_p1_QE->Clone();

  q2_m1_QE_rat->Divide(q2_nom_QE);
  q2_p1_QE_rat->Divide(q2_nom_QE);

  q2_p1_QE_rat->GetYaxis()->SetRangeUser(0, 2);
  q2_p1_QE_rat->Draw("HIST");
  q2_m1_QE_rat->Draw("HISTSAME");

  c1->SaveAs("VALORXSecPlots.pdf");
  c1->Clear();


  c1->Divide(1, 2);
  c1->cd(1);

  TH1D* q2_nom_RES = new TH1D("q2_nom_RES", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_nom_RES",
              "(stop >= 0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==3)",
              "GOFF");

  TH1D* q2_p1_RES = new TH1D("q2_p1_RES", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_p1_RES",
              "XSecWeights[0]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==3))",
              "GOFF");

  TH1D* q2_m1_RES = new TH1D("q2_m1_RES", ";Q^{2} (GeV^{2});Count", 100, 0, 2);

  EDeps->Draw("Q2_True >> q2_m1_RES",
              "XSecWeights[1]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==3))",
              "GOFF");

  q2_p1_RES->SetLineStyle(2);
  q2_p1_RES->SetLineColor(kBlue);
  q2_m1_RES->SetLineStyle(2);
  q2_m1_RES->SetLineColor(kRed);

  q2_p1_RES->Draw("EHIST");
  q2_nom_RES->Draw("EHISTSAME");
  q2_m1_RES->Draw("EHISTSAME");

  c1->cd(2);

  TH1D* q2_m1_RES_rat = q2_m1_RES->Clone();
  TH1D* q2_p1_RES_rat = q2_p1_RES->Clone();

  q2_m1_RES_rat->Divide(q2_nom_RES);
  q2_p1_RES_rat->Divide(q2_nom_RES);

  q2_p1_RES_rat->GetYaxis()->SetRangeUser(0, 2);
  q2_p1_RES_rat->Draw("HIST");
  q2_m1_RES_rat->Draw("HISTSAME");

  c1->SaveAs("VALORXSecPlots.pdf");
  c1->Clear();



  c1->Divide(1, 2);
  c1->cd(1);

  TH1D* q2_nom_RES = new TH1D("q2_nom_DIS", ";Q^{2} (GeV^{2});Count", 100, 0, 20);

  EDeps->Draw("Q2_True >> q2_nom_DIS",
              "(stop >= 0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==4)",
              "GOFF");

  TH1D* q2_p1_DIS = new TH1D("q2_p1_DIS", ";Q^{2} (GeV^{2});Count", 100, 0, 20);

  EDeps->Draw("Q2_True >> q2_p1_DIS",
              "XSecWeights[0]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==4))",
              "GOFF");

  TH1D* q2_m1_DIS = new TH1D("q2_m1_DIS", ";Q^{2} (GeV^{2});Count", 100, 0, 20);

  EDeps->Draw("Q2_True >> q2_m1_DIS",
              "XSecWeights[1]*((stop >= "
              "0)&&(PrimaryLepPDG==13)&&(GENIEInteractionTopology==4))",
              "GOFF");

  q2_p1_DIS->SetLineStyle(2);
  q2_p1_DIS->SetLineColor(kBlue);
  q2_m1_DIS->SetLineStyle(2);
  q2_m1_DIS->SetLineColor(kRed);

  q2_p1_DIS->Draw("EHIST");
  q2_nom_DIS->Draw("EHISTSAME");
  q2_m1_DIS->Draw("EHISTSAME");

  c1->cd(2);

  TH1D* q2_m1_DIS_rat = q2_m1_DIS->Clone();
  TH1D* q2_p1_DIS_rat = q2_p1_DIS->Clone();

  q2_m1_DIS_rat->Divide(q2_nom_DIS);
  q2_p1_DIS_rat->Divide(q2_nom_DIS);

  q2_p1_DIS_rat->GetYaxis()->SetRangeUser(0, 2);
  q2_p1_DIS_rat->Draw("HIST");
  q2_m1_DIS_rat->Draw("HISTSAME");

  c1->SaveAs("VALORXSecPlots.pdf");
  c1->Clear();

  c1->SaveAs("VALORXSecPlots.pdf]");
}
