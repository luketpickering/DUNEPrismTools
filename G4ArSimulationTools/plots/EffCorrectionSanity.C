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

  sel_enu_pos_fv_evrate->Divide(true_enu_pos_fv_evrate);
  sel_enu_pos_fv_evrate->SetTitle(";X (cm);E_{#nu} (GeV);Eff");
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

  corr_enu_pos_fv_evrate->Divide(true_enu_pos_fv_evrate);
  corr_enu_pos_fv_evrate->SetTitle(";X (cm);E_{#nu} (GeV);Eff");
  corr_enu_pos_fv_evrate->Draw("COLZ");
  corr_enu_pos_fv_evrate->GetZaxis()->SetRangeUser(0, 2);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corr_ehadr_pos_fv_evrate =
      new TH2D("corr_ehadr_pos_fv_evrate",
               "corr_ehadr_pos_fv_evrate;X (cm);E_{#nu} (GeV);Count", 195,
               -3750, 150, 50, 0, 10);

  EDeps->Draw("nu_4mom[3]:vtx[0] >> corr_ehadr_pos_fv_evrate",
              "DumbEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  corr_ehadr_pos_fv_evrate->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  corr_ehadr_pos_fv_evrate->Divide(true_enu_pos_fv_evrate);
  corr_ehadr_pos_fv_evrate->SetTitle(";X (cm);E_{#nu} (GeV);Eff");
  corr_ehadr_pos_fv_evrate->Draw("COLZ");
  corr_ehadr_pos_fv_evrate->GetZaxis()->SetRangeUser(0, 2);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  // 100, 0, 20, 80, -3750, 150
  TH2D* true_ehadr_pos_fv_evrate_samebin =
      new TH2D("true_ehadr_pos_fv_evrate_samebin",
               "true_ehadr_pos_fv_evrate_samebin;E_{#nu} (GeV);X (cm);Count",
               100, 0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:nu_4mom[3] "
      ">> true_ehadr_pos_fv_evrate_samebin",
      "((stop>=0)&&(PrimaryLepPDG==13))", "GOFF");

  true_ehadr_pos_fv_evrate_samebin->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* sel_ehadr_pos_fv_evrate_samebin =
      new TH2D("sel_ehadr_pos_fv_evrate_samebin",
               "sel_ehadr_pos_fv_evrate_samebin;E_{#nu} (GeV);X (cm);Count",
               100, 0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:nu_4mom[3] "
      ">> sel_ehadr_pos_fv_evrate_samebin",
      "((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  sel_ehadr_pos_fv_evrate_samebin->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corr_ehadr_pos_fv_evrate_samebin =
      new TH2D("corr_ehadr_pos_fv_evrate_samebin",
               "corr_ehadr_pos_fv_evrate_samebin;E_{#nu} (GeV);X (cm);Count",
               100, 0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:nu_4mom[3] "
      ">> corr_ehadr_pos_fv_evrate_samebin",
      "DumbEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  corr_ehadr_pos_fv_evrate_samebin->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  corr_ehadr_pos_fv_evrate_samebin->SetTitle(";E_{#nu} (GeV);X (cm);Eff");
  corr_ehadr_pos_fv_evrate_samebin->Divide(true_ehadr_pos_fv_evrate_samebin);
  corr_ehadr_pos_fv_evrate_samebin->Draw("COLZ");
  corr_ehadr_pos_fv_evrate_samebin->GetZaxis()->SetRangeUser(0, 2);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* true_vis_smear =
      new TH2D("true_vis_smear", "true_vis_smear;E_{#nu} (GeV);E_{Rec.};Count",
               50, 0, 20, 50, 0, 20);

  EDeps->Draw(
      "PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto:nu_4mom[3] "
      ">> true_vis_smear",
      "((stop>=0)&&(PrimaryLepPDG==13))", "GOFF");

  true_vis_smear->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* true_vis_smear_sel = new TH2D(
      "true_vis_smear_sel", "true_vis_smear_sel;E_{#nu} (GeV);E_{Rec.};Count",
      50, 0, 20, 50, 0, 20);

  EDeps->Draw(
      "PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto:nu_4mom[3] "
      ">> true_vis_smear_sel",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh)&&("
      "HadrShowerContainedInFV))",
      "GOFF");

  true_vis_smear_sel->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  // 100, 0, 20, 80, -3750, 150
  TH2D* true_ehadr_pos_fv_evrate_samebin_evis = new TH2D(
      "true_ehadr_pos_fv_evrate_samebin_evis",
      "true_ehadr_pos_fv_evrate_samebin_evis;E_{Rec.} (GeV);X (cm);Count", 100,
      0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto "
      ">> true_ehadr_pos_fv_evrate_samebin_evis",
      "((stop>=0)&&(PrimaryLepPDG==13))", "GOFF");

  true_ehadr_pos_fv_evrate_samebin_evis->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* sel_ehadr_pos_fv_evrate_samebin_evis = new TH2D(
      "sel_ehadr_pos_fv_evrate_samebin_evis",
      "sel_ehadr_pos_fv_evrate_samebin_evis;E_{Rec.} (GeV);X (cm);Count", 100,
      0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto "
      ">> sel_ehadr_pos_fv_evrate_samebin_evis",
      "((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  sel_ehadr_pos_fv_evrate_samebin_evis->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corr_ehadr_pos_fv_evrate_samebin_evis = new TH2D(
      "corr_ehadr_pos_fv_evrate_samebin_evis",
      "corr_ehadr_pos_fv_evrate_samebin_evis;E_{Rec.} (GeV);X (cm);Count", 100,
      0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto "
      ">> corr_ehadr_pos_fv_evrate_samebin_evis",
      "DumbEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  corr_ehadr_pos_fv_evrate_samebin_evis->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  corr_ehadr_pos_fv_evrate_samebin_evis->SetTitle(";E_{#nu} (GeV);X (cm);Eff");
  corr_ehadr_pos_fv_evrate_samebin_evis->Divide(
      true_ehadr_pos_fv_evrate_samebin_evis);
  corr_ehadr_pos_fv_evrate_samebin_evis->Draw("COLZ");
  corr_ehadr_pos_fv_evrate_samebin_evis->GetZaxis()->SetRangeUser(0, 2);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* sel_ehadr_pos_fvonly_evrate_samebin_evis = new TH2D(
      "sel_ehadr_pos_fvonly_evrate_samebin_evis",
      "sel_ehadr_pos_fvonly_evrate_samebin_evis;E_{Rec.} (GeV);X (cm);Count",
      100, 0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto "
      ">> sel_ehadr_pos_fvonly_evrate_samebin_evis",
      "((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  sel_ehadr_pos_fvonly_evrate_samebin_evis->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  TH2D* corr_ehadr_pos_fvonly_evrate_samebin_evis = new TH2D(
      "corr_ehadr_pos_fvonly_evrate_samebin_evis",
      "corr_ehadr_pos_fvonly_evrate_samebin_evis;E_{Rec.} (GeV);X (cm);Count",
      100, 0, 20, 80, -3750, 150);

  EDeps->Draw(
      "vtx[0]:PrimaryLep_4mom[3]+TotalNonlep_Dep_FV+TotalNonlep_Dep_veto "
      ">> corr_ehadr_pos_fvonly_evrate_samebin_evis",
      "DumbEffWeight_FV*((stop>=0)&&(PrimaryLepPDG==13)&&("
      "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
      "GOFF");

  corr_ehadr_pos_fvonly_evrate_samebin_evis->Draw("COLZ");
  c1->SaveAs("EffCorrectionSanityPlots.pdf");
  c1->Clear();

  corr_ehadr_pos_fvonly_evrate_samebin_evis->SetTitle(
      ";E_{#nu} (GeV);X (cm);Eff");
  corr_ehadr_pos_fvonly_evrate_samebin_evis->Divide(
      true_ehadr_pos_fv_evrate_samebin_evis);
  corr_ehadr_pos_fvonly_evrate_samebin_evis->Draw("COLZ");
  corr_ehadr_pos_fvonly_evrate_samebin_evis->GetZaxis()->SetRangeUser(0, 2);
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
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> "
      "ERecSel_hadr_true",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(HadrShowerContainedInFV))", "GOFF");

  TH1D* ERecCorr_hadr_true = new TH1D(
      "ERecCorr_hadr_true", "ERecCorr_hadr_true;E (GeV);Count", 50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> "
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
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> "
      "ERecSel_hadr_true_mu",
      "((stop>=0)&&(PrimaryLepPDG==13)&&(LepExit_AboveThresh)&&("
      "HadrShowerContainedInFV))",
      "GOFF");

  TH1D* ERecCorr_hadr_true_mu =
      new TH1D("ERecCorr_hadr_true_mu", "ERecCorr_hadr_true_mu;E (GeV);Count",
               50, 0, 10);
  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> "
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
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecSel",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecCorr = new TH1D("ERecCorr", "ERecCorr;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecCorr",
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

  TH1D* ETrueSel_offsetcorr =
      new TH1D("ETrueSel_offsetcorr", "ETrueSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_offsetcorr",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel_offsetcorr =
      new TH1D("ERecSel_offsetcorr", "ERecSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecSel_offsetcorr",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecCorr_offsetcorr = new TH1D(
      "ERecCorr_offsetcorr", "ERecCorr (edep,x);E (GeV);Count", 50, 0, 10);
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecCorr_offsetcorr",
              "DumbEffWeight*((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  ETrueSel_offsetcorr->SetLineColor(kRed);
  ERecSel_offsetcorr->SetLineColor(kBlack);
  ERecCorr_offsetcorr->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_offsetcorr->Draw("SAME");
  ERecSel_offsetcorr->Draw("SAME");
  ERecCorr_offsetcorr->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  TH1D* ETrueSel_offsetcorr_fv =
      new TH1D("ETrueSel_offsetcorr_fv", "ETrueSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("nu_4mom[3] >> ETrueSel_offsetcorr_fv",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecSel_offsetcorr_fv =
      new TH1D("ERecSel_offsetcorr_fv", "ERecSel;E (GeV);Count", 50, 0, 10);
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecSel_offsetcorr_fv",
              "((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  TH1D* ERecCorr_offsetcorr_fv = new TH1D(
      "ERecCorr_offsetcorr_fv", "ERecCorr (edep,x);E (GeV);Count", 50, 0, 10);
  EDeps->Draw("TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3] >> ERecCorr_offsetcorr_fv",
              "DumbEffWeight_FV*((stop>=0)&&(PrimaryLepPDG==13)&&("
              "LepExit_AboveThresh)&&(HadrShowerContainedInFV))",
              "GOFF");

  ETrueSel_offsetcorr_fv->SetLineColor(kRed);
  ERecSel_offsetcorr_fv->SetLineColor(kBlack);
  ERecCorr_offsetcorr_fv->SetLineColor(kMagenta);

  ETrue->Draw();
  ETrueSel_offsetcorr_fv->Draw("SAME");
  ERecSel_offsetcorr_fv->Draw("SAME");
  ERecCorr_offsetcorr_fv->Draw("SAME");

  c1->BuildLegend(0.7, 0.7, 0.95, 0.95);
  c1->SaveAs("EffCorrectionSanityPlots.pdf");

  c1->SaveAs("EffCorrectionSanityPlots.pdf]");
}
