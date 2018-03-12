void EffCorrectionSanity(char const* inps, char const* wfile) {
  gStyle->SetOptStat(0);
  TH1D::SetDefaultSumw2(true);

  TChain* EDeps = new TChain("EDeps");
  EDeps->Add(inps);
  EDeps->AddFriend("EffWeights", wfile);

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SaveAs("EffCorrectionSanityPlots.pdf[");

  TH2D* true_enu_pos_fv_evrate =
      new TH2D("true_enu_pos_fv_evrate",
               "true_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195, -3750,
               150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> true_enu_pos_fv_evrate",
              "(stop>=0)&&(PrimaryLepPDG==13)", "GOFF");

  true_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* selmu_enu_pos_fv_evrate =
      new TH2D("selmu_enu_pos_fv_evrate",
               "selmu_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195, -3750,
               150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> selmu_enu_pos_fv_evrate",
              "((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh))",
              "GOFF");

  selmu_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corrmu_enu_pos_fv_evrate =
      new TH2D("corrmu_enu_pos_fv_evrate",
               "corrmu_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195,
               -3750, 150, 50, 0, 10);

  EDeps->Draw(
      "nu_4mom[3]:vtx[0] >> corrmu_enu_pos_fv_evrate",
      "MuEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh))",
      "GOFF");

  corrmu_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* selhadr_enu_pos_fv_evrate =
      new TH2D("selhadr_enu_pos_fv_evrate",
               "selhadr_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195,
               -3750, 150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> selhadr_enu_pos_fv_evrate",
              "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))",
              "GOFF");

  selhadr_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* sel_enu_pos_fv_evrate =
      new TH2D("sel_enu_pos_fv_evrate",
               "sel_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195, -3750,
               150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> sel_enu_pos_fv_evrate",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  sel_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  sel_enu_pos_fv_evrate->SetTitle(";X (cm);E_{#nu} (GeV);Eff");
  sel_enu_pos_fv_evrate->Divide(true_enu_pos_fv_evrate);
  sel_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corr_enu_pos_fv_evrate =
      new TH2D("corr_enu_pos_fv_evrate",
               "corr_enu_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195, -3750,
               150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> corr_enu_pos_fv_evrate",
              "EffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  corr_enu_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  corr_enu_pos_fv_evrate->SetTitle(";X (cm);E_{#nu} (GeV);Eff");
  corr_enu_pos_fv_evrate->Divide(true_enu_pos_fv_evrate);
  corr_enu_pos_fv_evrate->Draw("COLZ");
  corr_enu_pos_fv_evrate->GetZaxis()->SetRangeUser(0, 2);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH1D* ETrue = new TH1D("ETrue", "ETrue;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrue", "((stop>=0)&&(PrimaryLepPDG==13))",
              "GOFF");

  TH1D* ETrueSel_mu =
      new TH1D("ETrueSel_mu", "ETrueSel_mu;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_mu",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh))",
              "GOFF");

  TH1D* ERecSel_mu =
      new TH1D("ERecSel_mu", "ERecSel_mu;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecSel_mu",
      "((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh))",
      "GOFF");

  TH1D* ERecCorr_mu =
      new TH1D("ERecCorr_mu", "ERecCorr_mu;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecCorr_mu",
      "MuEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh))",
      "GOFF");

  ETrueSel_mu->SetLineColor(kRed);
  ERecSel_mu->SetLineColor(kBlack);
  ERecCorr_mu->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_mu->Draw("SAME");
  ERecSel_mu->Draw("SAME");
  ERecCorr_mu->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  TH1D* ETrueSel_hadr =
      new TH1D("ETrueSel_hadr", "ETrueSel_hadr;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_hadr",
              "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel_hadr =
      new TH1D("ERecSel_hadr", "ERecSel_hadr;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecSel_hadr",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))", "GOFF");

  TH1D* ERecCorr_hadr =
      new TH1D("ERecCorr_hadr", "ERecCorr_hadr;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecCorr_hadr",
      "HadrEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV)"
      ")",
      "GOFF");

  ETrueSel_hadr->SetLineColor(kRed);
  ERecSel_hadr->SetLineColor(kBlack);
  ERecCorr_hadr->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_hadr->Draw("SAME");
  ERecSel_hadr->Draw("SAME");
  ERecCorr_hadr->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  TH1D* ETrueSel_hadr_true = new TH1D(
      "ETrueSel_hadr_true", "ETrueSel_hadr_true;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_hadr_true",
              "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel_hadr_true = new TH1D(
      "ERecSel_hadr_true", "ERecSel_hadr_true;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecSel_hadr_true",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))", "GOFF");

  TH1D* ERecCorr_hadr_true = new TH1D(
      "ERecCorr_hadr_true", "ERecCorr_hadr_true;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecCorr_hadr_true",
      "HadrEffWeight_true*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "HadrShowerContainedInFV)"
      ")",
      "GOFF");

  ETrueSel_hadr_true->SetLineColor(kRed);
  ERecSel_hadr_true->SetLineColor(kBlack);
  ERecCorr_hadr_true->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_hadr_true->Draw("SAME");
  ERecSel_hadr_true->Draw("SAME");
  ERecCorr_hadr_true->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  TH1D* ETrueSel_hadr_true_mu =
      new TH1D("ETrueSel_hadr_true_mu", "ETrueSel_hadr_true_mu;E (GeV);Count",
               50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_hadr_true_mu",
              "((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh)&&("
              "HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel_hadr_true_mu = new TH1D(
      "ERecSel_hadr_true_mu", "ERecSel_hadr_true_mu;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecSel_hadr_true_mu",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh)&&("
      "HadrShowerContainedInFV))",
      "GOFF");

  TH1D* ERecCorr_hadr_true_mu =
      new TH1D("ERecCorr_hadr_true_mu", "ERecCorr_hadr_true_mu;E (GeV);Count",
               50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> "
      "ERecCorr_hadr_true_mu",
      "HadrEffWeight_true*MuEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  ETrueSel_hadr_true_mu->SetLineColor(kRed);
  ERecSel_hadr_true_mu->SetLineColor(kBlack);
  ERecCorr_hadr_true_mu->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_hadr_true_mu->Draw("SAME");
  ERecSel_hadr_true_mu->Draw("SAME");
  ERecCorr_hadr_true_mu->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  TH1D* ETrueSel = new TH1D("ETrueSel", "ETrueSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel = new TH1D("ERecSel", "ERecSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> ERecSel",
      "((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  TH1D* ERecCorr = new TH1D("ERecCorr", "ERecCorr;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+PrimaryLep_4mom[3] >> ERecCorr",
      "EffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  ETrueSel->SetLineColor(kRed);
  ERecSel->SetLineColor(kBlack);
  ERecCorr->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel->Draw("SAME");
  ERecSel->Draw("SAME");
  ERecCorr->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  c1->SaveAs("EffCorrectionSanityPlots.pdf]");
}
