void SanityPlots_EDeps_Tree(const char* inps) {
  TChain* EDeps = new TChain("EDeps");
  EDeps->Add(inps);

  EDeps->SetAlias("nu_3mag",
                  "sqrt(nu_4mom[0]*nu_4mom[0]+nu_4mom[1]*nu_4mom[1]+nu_4mom[2]*"
                  "nu_4mom[2])");
  EDeps->SetAlias("totalfs_3mag",
                  "sqrt(TotalFS_3mom[0]*TotalFS_3mom[0]+TotalFS_3mom[1]*"
                  "TotalFS_3mom[1]+TotalFS_3mom[2]*TotalFS_3mom[2])");
  EDeps->SetAlias(
      "PrimaryLep_3mag",
      "sqrt(PrimaryLep_4mom[0]*PrimaryLep_4mom[0]+PrimaryLep_4mom[1]*"
      "PrimaryLep_4mom[1]+PrimaryLep_4mom[2]*PrimaryLep_4mom[2])");
  EDeps->SetAlias("diff_nu_tfs_x", "(TotalFS_3mom[0]-nu_4mom[0])");
  EDeps->SetAlias("diff_nu_tfs_y", "(TotalFS_3mom[1]-nu_4mom[1])");
  EDeps->SetAlias("diff_nu_tfs_z", "(TotalFS_3mom[2]-nu_4mom[2])");
  EDeps->SetAlias("diff_nu_tfs_3mag",
                  "sqrt(diff_nu_tfs_x*diff_nu_tfs_x+diff_nu_tfs_y*diff_nu_tfs_"
                  "y+diff_nu_tfs_z*diff_nu_tfs_z)");
  EDeps->SetAlias("ENuQERec",
                  "( (2 * 0.938 * PrimaryLep_4mom[3] - 0.106*0.106 + "
                  "0.938*0.938 - 0.938*0.938)/ (2*(0.938 - "
                  "PrimaryLep_4mom[3] + PrimaryLep_3mag * "
                  "(PrimaryLep_4mom[2]/PrimaryLep_3mag))) )");

  TCanvas* c1 = new TCanvas("c1", "");

  std::cout << "Making enu" << std::endl;
  TH1D* enu = new TH1D("enu", ";E_{#nu} (GeV);Count", 250, 0, 25);

  EDeps->Draw("nu_4mom[3] >> enu", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  enu->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf[");
  c1->Clear();

  std::cout << "Making nu_pdg" << std::endl;
  TH1D* nu_pdg = new TH1D("nu_pdg", "PDG_{#nu};PDG_{#nu};Count", 6, 11, 17);

  EDeps->Draw("nu_PDG >> nu_pdg", "(stop >= 0)", "GOFF");

  std::cout << "Making nu_pdg_neg" << std::endl;
  TH1D* nu_pdg_neg =
      new TH1D("nu_pdg_neg", "abs(PDG_{#bar{#nu}});PDG_{#nu};Count", 6, 11, 17);

  EDeps->Draw("abs(nu_PDG) >> nu_pdg_neg", "nu_PDG<0", "GOFF");

  nu_pdg->Draw();
  nu_pdg_neg->SetLineColor(kRed);
  nu_pdg_neg->Draw("SAME");

  c1->BuildLegend(0.2, 0.8, 0.35, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making vtx_pos_xy" << std::endl;
  TH2D* vtx_pos_xy = new TH2D("vtx_pos_xy", ";VtxPos_{x} (cm);VtxPos_{y} (cm)",
                              880, -4000, 400, 80, -200, 200);

  EDeps->Draw("vtx[1]:vtx[0] >> vtx_pos_xy", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  vtx_pos_xy->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making vtx_pos_yz" << std::endl;
  TH2D* vtx_pos_yz = new TH2D("vtx_pos_yz", ";VtxPos_{z} (cm);VtxPos_{y} (cm)",
                              120, -300, 300, 80, -200, 200);

  EDeps->Draw("vtx[1]:vtx[2] >> vtx_pos_yz", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  vtx_pos_yz->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making q2" << std::endl;
  TH1D* q2 = new TH1D("q2", ";Q^{2} (GeV^{2});Count", 100, 0, 10);

  EDeps->Draw("Q2_True >> q2", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  q2->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making y" << std::endl;
  TH1D* y = new TH1D("y", "All;Elasticity;Count", 100, 0, 1);

  EDeps->Draw("(1-y_True) >> y", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  y->Draw();

  TH1D* y_hadr_cont =
      new TH1D("y_hadr_cont", "HadrContained;Elasticity;Count", 100, 0, 1);

  EDeps->Draw("(1-y_True) >> y_hadr_cont",
              "(stop >= 0)&&(HadrShowerContainedInFV==1)&&(PrimaryLepPDG==13)",
              "GOFF");

  y_hadr_cont->SetLineColor(kRed);
  y_hadr_cont->Draw("SAME");

  TH1D* y_mu_cont =
      new TH1D("y_mu_cont", "MuContained;Elasticity;Count", 100, 0, 1);

  EDeps->Draw(
      "(1-y_True) >> y_mu_cont",
      "(stop >= 0)&&(PrimaryLeptonContainedInFV==1)&&(PrimaryLepPDG==13)",
      "GOFF");

  y_mu_cont->SetLineColor(kBlue);
  y_mu_cont->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making wrest" << std::endl;
  TH1D* wrest = new TH1D("wrest", ";W_{rest} (GeV);Count", 100, 0.8, 3);

  EDeps->Draw("W_Rest >> wrest", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  wrest->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NLep" << std::endl;
  TH1D* NLep = new TH1D("NLep", ";N_{l};Count", 5, 0, 5);

  EDeps->Draw("NLep >> NLep", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NLep->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NPi0" << std::endl;
  TH1D* NPi0 = new TH1D("NPi0", ";N_{#pi^{0}};Count", 10, 0, 10);

  EDeps->Draw("NPi0 >> NPi0", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NPi0->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NPiC" << std::endl;
  TH1D* NPiC = new TH1D("NPiC", ";N_{#pi^{#pm}};Count", 10, 0, 10);

  EDeps->Draw("NPiC >> NPiC", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NPiC->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NProton" << std::endl;
  TH1D* NProton = new TH1D("NProton", ";N_{proton};Count", 20, 0, 20);

  EDeps->Draw("NProton >> NProton", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NProton->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NNeutron" << std::endl;
  TH1D* NNeutron = new TH1D("NNeutron", ";N_{neutron};Count", 20, 0, 20);

  EDeps->Draw("NNeutron >> NNeutron", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  NNeutron->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NGamma" << std::endl;
  TH1D* NGamma = new TH1D("NGamma", ";N_{#gamma};Count", 5, 0, 5);

  EDeps->Draw("NGamma >> NGamma", "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NGamma->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NBaryonicRes" << std::endl;
  TH1D* NBaryonicRes =
      new TH1D("NBaryonicRes", ";N_{baryonres};Count", 5, 0, 5);

  EDeps->Draw("NBaryonicRes >> NBaryonicRes",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  NBaryonicRes->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making primlep_e" << std::endl;
  TH2D* primlep_e = new TH2D(
      "primlep_e", ";T_{l} + N_{l}#timesM_{l} (GeV);EDep_{l,prim.} (GeV)", 100,
      0, 10, 100, 0, 10);

  EDeps->Draw("LepDep_FV+LepDep_veto:PrimaryLep_4mom[3] >> primlep_e",
              "(PrimaryLepPDG==13)&&(NLep==1)", "GOFF");

  primlep_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making primlep_e_and_desc" << std::endl;
  TH2D* primlep_e_and_desc =
      new TH2D("primlep_e_and_desc",
               ";T_{l} + N_{l}#timesM_{l} (GeV);EDep_{l+desc.} (GeV)", 100, 0,
               10, 100, 0, 10);

  EDeps->Draw(
      "LepDep_FV+LepDep_veto+LepDepDescendent_FV+LepDepDescendent_veto:"
      "PrimaryLep_4mom[3] >> primlep_e_and_desc",
      "(PrimaryLepPDG==13)&&(NLep==1)", "GOFF");

  primlep_e_and_desc->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making pi0_e" << std::endl;
  TH2D* pi0_e = new TH2D(
      "pi0_e",
      ";T_{#pi^{0}} + N_{#pi^{0}}#timesM_{#pi^{0}} (GeV);EDep_{#pi^{0}} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw("Pi0Dep_FV+Pi0Dep_veto:EKinPi0_True+EMassPi0_True >> pi0_e", "",
              "GOFF");

  pi0_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making pic_e" << std::endl;
  TH2D* pic_e =
      new TH2D("pic_e",
               ";T_{#pi^{#pm}} + "
               "N_{#pi^{#pm}}#timesM_{#pi^{#pm}} (GeV);EDep_{#pi^{#pm}} (GeV)",
               100, 0, 10, 100, 0, 10);

  EDeps->Draw("PiCDep_FV+PiCDep_veto:EKinPiC_True+EMassPiC_True >> pic_e", "",
              "GOFF");

  pic_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making neutron_e" << std::endl;
  TH2D* neutron_e = new TH2D(
      "neutron_e",
      ";T_{neutron} + N_{neutron}#timesM_{neutron} (GeV);EDep_{neutron} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "NeutronDep_FV+NeutronDep_veto:EKinNeutron_True+EMassNeutron_True >> "
      "neutron_e",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  neutron_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making proton_e" << std::endl;
  TH2D* proton_e = new TH2D(
      "proton_e",
      ";T_{proton} + N_{proton}#timesM_{proton} (GeV);EDep_{proton} (GeV)", 100,
      0, 10, 100, 0, 10);

  EDeps->Draw(
      "ProtonDep_FV+ProtonDep_veto:EKinProton_True+EMassProton_True >> "
      "proton_e",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  proton_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making neutron_e_nom" << std::endl;
  TH2D* neutron_e_nom = new TH2D(
      "neutron_e_nom",
      ";T_{neutron} + N_{neutron}#timesM_{neutron} (GeV);EDep_{neutron} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "NeutronDep_FV+NeutronDep_veto:EKinNeutron_True >> "
      "neutron_e_nom",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  neutron_e_nom->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making proton_e_nom" << std::endl;
  TH2D* proton_e_nom = new TH2D(
      "proton_e_nom",
      ";T_{proton} + N_{proton}#timesM_{proton} (GeV);EDep_{proton} (GeV)", 100,
      0, 10, 100, 0, 10);

  EDeps->Draw(
      "ProtonDep_FV+ProtonDep_veto:EKinProton_True >> "
      "proton_e_nom",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  proton_e_nom->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_enu_2D" << std::endl;
  TH2D* efs_vs_enu_2D = new TH2D(
      "efs_vs_enu_2D",
      ";#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + #sumT_{nucl.} + E_{#mu} "
      "(GeV);E_{#nu} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "nu_4mom[3]:ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3]"
      " >> efs_vs_enu_2D",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_enu_2D->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_enu_2D_nobarres" << std::endl;
  TH2D* efs_vs_enu_2D_nobarres =
      new TH2D("efs_vs_enu_2D_nobarres",
               ";#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + #sumT_{nucl.} "
               "- N_{baryonres}#timesm_{nuc} + E_{#mu} "
               "(GeV);E_{#nu} (GeV)",
               100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "nu_4mom[3]:ERecProxy_True"
      " >> efs_vs_enu_2D_nobarres",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_enu_2D_nobarres->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_enu_1D_0pi" << std::endl;
  TH1D* efs_vs_enu_1D_0pi =
      new TH1D("efs_vs_enu_1D_0pi",
               "0Pi;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_0pi",
      "(stop >= 0)&&(Is0Pi)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_enu_1D_0pi->SetLineColor(kRed);
  efs_vs_enu_1D_0pi->Draw();

  std::cout << "Making efs_vs_enu_1D_1Pi" << std::endl;
  TH1D* efs_vs_enu_1D_1Pi =
      new TH1D("efs_vs_enu_1D_1Pi",
               "1Pi;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_1Pi",
      "(stop >= 0)&&(Is1Pi)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_enu_1D_1Pi->SetLineColor(kBlue);
  efs_vs_enu_1D_1Pi->Draw("SAME");

  std::cout << "Making efs_vs_enu_1D_Other" << std::endl;
  TH1D* efs_vs_enu_1D_Other =
      new TH1D("efs_vs_enu_1D_Other",
               "Other;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_Other",
      "(stop >= 0)&&(IsNPi||IsOther)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_enu_1D_Other->SetLineColor(kMagenta);
  efs_vs_enu_1D_Other->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);
  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_enu_1D_QE" << std::endl;
  TH1D* efs_vs_enu_1D_QE =
      new TH1D("efs_vs_enu_1D_QE",
               "QE;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_QE",
      "(stop >= 0)&&(GENIEInteractionTopology==1)&&(PrimaryLepPDG==13)",
      "GOFF");

  efs_vs_enu_1D_QE->SetLineColor(kRed);

  std::cout << "Making efs_vs_enu_1D_RES" << std::endl;
  TH1D* efs_vs_enu_1D_RES =
      new TH1D("efs_vs_enu_1D_RES",
               "RES;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_RES",
      "(stop >= 0)&&(GENIEInteractionTopology==3)&&(PrimaryLepPDG==13)",
      "GOFF");

  efs_vs_enu_1D_RES->SetLineColor(kBlue);

  std::cout << "Making efs_vs_enu_1D_DIS" << std::endl;
  TH1D* efs_vs_enu_1D_DIS =
      new TH1D("efs_vs_enu_1D_DIS",
               "DIS;#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
               "#sumT_{nucl.} + E_{#mu} - E_{#nu} (GeV); Count",
               100, -0.1, 0.1);

  EDeps->Draw(
      "nu_4mom[3]-(ENonPrimaryLep_KinNucleonTotalOther_True+PrimaryLep_4mom[3])"
      " >> efs_vs_enu_1D_DIS",
      "(stop >= 0)&&(GENIEInteractionTopology==4)&&(PrimaryLepPDG==13)",
      "GOFF");

  efs_vs_enu_1D_DIS->SetLineColor(kMagenta);

  efs_vs_enu_1D_QE->Draw();
  efs_vs_enu_1D_DIS->Draw("SAME");
  efs_vs_enu_1D_RES->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);
  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enuQERec_1D_QE" << std::endl;
  TH1D* enuQERec_1D_QE = new TH1D(
      "enuQERec_1D_QE", "QE;ERec_{QE}/E_{#nu} - 1; Count", 100, -1.5, 1);

  EDeps->Draw(
      "(ENuQERec/nu_4mom[3]) - 1"
      " >> enuQERec_1D_QE",
      "(stop >= 0)&&(GENIEInteractionTopology==1)&&(PrimaryLepPDG==13)",
      "GOFF");

  enuQERec_1D_QE->SetLineColor(kRed);

  std::cout << "Making enuQERec_1D_RES" << std::endl;
  TH1D* enuQERec_1D_RES = new TH1D(
      "enuQERec_1D_RES", "RES;ERec_{QE}/E_{#nu} - 1; Count", 100, -1.5, 1);

  EDeps->Draw(
      "(ENuQERec/nu_4mom[3]) - 1"
      " >> enuQERec_1D_RES",
      "(stop >= 0)&&(GENIEInteractionTopology==3)&&(PrimaryLepPDG==13)",
      "GOFF");

  enuQERec_1D_RES->SetLineColor(kBlue);

  std::cout << "Making enuQERec_1D_DIS" << std::endl;
  TH1D* enuQERec_1D_DIS = new TH1D(
      "enuQERec_1D_DIS", "DIS;ERec_{QE}/E_{#nu} - 1; Count", 100, -1.5, 1);

  EDeps->Draw(
      "(ENuQERec/nu_4mom[3]) - 1"
      " >> enuQERec_1D_DIS",
      "(stop >= 0)&&(GENIEInteractionTopology==4)&&(PrimaryLepPDG==13)",
      "GOFF");

  enuQERec_1D_DIS->SetLineColor(kMagenta);

  enuQERec_1D_QE->Draw();
  enuQERec_1D_DIS->Draw("SAME");
  enuQERec_1D_RES->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);
  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making miss_mom_p_ct" << std::endl;
  TH2D* miss_mom_p_ct = new TH2D(
      "miss_mom_p_ct",
      ";p_{N,i} = |p_{#sumFS} - p_{#nu}| (GeV/#it{c});cos(#theta_{p_{N,i}})",
      50, 0, 1, 50, -1, 1);

  EDeps->Draw(
      "(diff_nu_tfs_z/"
      "diff_nu_tfs_3mag):sqrt(diff_nu_tfs_x*diff_nu_tfs_x+diff_nu_tfs_z*diff_"
      "nu_tfs_z)"
      " >> miss_mom_p_ct",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  miss_mom_p_ct->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_dep" << std::endl;
  TH2D* efs_vs_dep = new TH2D("efs_vs_dep",
                              ";#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + "
                              "#sumT_{nucl.} + E_{#mu} (GeV);EDep_{det} (GeV)",
                              100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto+"
      "LepDepDescendent_FV+LepDepDescendent_veto:ENonPrimaryLep_"
      "KinNucleonTotalOther_True+PrimaryLep_4mom[3]"
      " >> efs_vs_dep",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_dep->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making efs_vs_dep_FV" << std::endl;
  TH2D* efs_vs_dep_FV = new TH2D(
      "efs_vs_dep_FV",
      ";#sumE_{#pi} + #sumE_{#gamma} + #sumE_{Other} + #sumT_{nucl.} + E_{#mu} "
      "(GeV);EDep_{FV} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+LepDep_FV+LepDepDescendent_FV:ENonPrimaryLep_"
      "KinNucleonTotalOther_True+PrimaryLep_4mom[3] >> efs_vs_dep_FV",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  efs_vs_dep_FV->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_dep" << std::endl;
  TH2D* enu_vs_dep = new TH2D("enu_vs_dep", ";E_{#nu} (GeV);EDep_{det} (GeV)",
                              100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto+"
      "LepDepDescendent_FV+LepDepDescendent_veto:nu_4mom[3]"
      " >> enu_vs_dep",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  enu_vs_dep->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_dep_FV" << std::endl;
  TH2D* enu_vs_dep_FV =
      new TH2D("enu_vs_dep_FV", ";E_{#nu} (GeV);EDep_{FV} (GeV)", 100, 0, 10,
               100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+LepDep_FV+LepDepDescendent_FV:nu_4mom[3] >> "
      "enu_vs_dep_FV",
      "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  enu_vs_dep_FV->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_erec_mu_cont" << std::endl;
  TH2D* enu_vs_erec_mu_cont =
      new TH2D("enu_vs_erec_mu_cont",
               ";E_{#nu} (GeV);E_{rec} = EDep_{#mu} + EDep_{hadr} (GeV)", 100,
               0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto+"
      "LepDepDescendent_FV+LepDepDescendent_veto:nu_4mom[3]"
      " >> enu_vs_erec_mu_cont",
      "(stop >= 0)&&(LepExit==0)&&(PrimaryLepPDG==13)", "GOFF");

  enu_vs_erec_mu_cont->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_erec_mu_exit" << std::endl;
  TH2D* enu_vs_erec_mu_exit =
      new TH2D("enu_vs_erec_mu_exit",
               "; EDep_{FV} (GeV);E_{rec} = E_{#mu} + EDep_{hadr} (GeV)", 100,
               0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3]:nu_4mom[3] "
      ">> enu_vs_erec_mu_exit",
      "(stop >= 0)&&(LepExit_AboveThresh==1)&&(PrimaryLepPDG==13)", "GOFF");

  enu_vs_erec_mu_exit->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_erec_mu_cont_hadr_cont" << std::endl;
  TH2D* enu_vs_erec_mu_cont_hadr_cont =
      new TH2D("enu_vs_erec_mu_cont_hadr_cont",
               ";E_{#nu} (GeV);E_{rec} = EDep_{#mu} + EDep_{hadr} (GeV)", 100,
               0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+LepDep_FV+LepDep_veto+LepDepDescendent_FV+"
      "LepDepDescendent_veto:nu_4mom[3]"
      " >> enu_vs_erec_mu_cont_hadr_cont",
      "(stop >= "
      "0)&&(LepExit==0)&&(HadrShowerContainedInFV==1)&&(PrimaryLepPDG==13)",
      "GOFF");

  enu_vs_erec_mu_cont_hadr_cont->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_erec_mu_exit_hadr_cont" << std::endl;
  TH2D* enu_vs_erec_mu_exit_hadr_cont =
      new TH2D("enu_vs_erec_mu_exit_hadr_cont",
               "; EDep_{FV} (GeV);E_{rec} = E_{#mu} + EDep_{hadr} (GeV)", 100,
               0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3]:nu_4mom[3] "
      ">> enu_vs_erec_mu_exit_hadr_cont",
      "(stop >= "
      "0)&&(LepExit_AboveThresh==1)&&(HadrShowerContainedInFV==1)&&("
      "PrimaryLepPDG==13)",
      "GOFF");

  enu_vs_erec_mu_exit_hadr_cont->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making erec_bias_mu_exit" << std::endl;
  TH1D* erec_bias_mu_exit =
      new TH1D("erec_bias_mu_exit", "LepExit_AllHadr;E_{rec}/E_{#nu} - 1;Count",
               100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3])/"
      "nu_4mom[3]) - 1 >> erec_bias_mu_exit",
      "(stop >= "
      "0)&&(LepExit_AboveThresh==1)&&(PrimaryLepPDG==13)",
      "GOFF");

  erec_bias_mu_exit->Draw();

  TH1D* erec_bias_mu_cont =
      new TH1D("erec_bias_mu_cont", "LepCont_AllHadr;E_{rec}/E_{#nu} - 1;Count",
               100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto+"
      "LepDepDescendent_FV+LepDepDescendent_veto)/"
      "nu_4mom[3]) - 1 >> erec_bias_mu_cont",
      "(stop >= "
      "0)&&(LepExit==0)&&(PrimaryLepPDG==13)",
      "GOFF");

  erec_bias_mu_cont->SetLineColor(kRed);
  erec_bias_mu_cont->Draw("SAME");

  TH1D* erec_bias_mu_exit_hadrcont =
      new TH1D("erec_bias_mu_exit_hadrcont",
               "LepExit_HadrCont;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3])/nu_4mom[3]) - 1 >> "
      "erec_bias_mu_exit_hadrcont",
      "(stop >= "
      "0)&&(LepExit_AboveThresh==1)&&(HadrShowerContainedInFV==1)&&("
      "PrimaryLepPDG==13)",
      "GOFF");

  erec_bias_mu_exit_hadrcont->SetLineColor(kBlue);
  erec_bias_mu_exit_hadrcont->Draw("SAME");

  TH1D* erec_bias_mu_cont_hadrcont =
      new TH1D("erec_bias_mu_cont_hadrcont",
               "LepCont_HadrCont;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+LepDep_FV+LepDep_veto+LepDepDescendent_FV+"
      "LepDepDescendent_veto)/nu_4mom[3]) - 1 >> "
      "erec_bias_mu_cont_hadrcont",
      "(stop >= "
      "0)&&(LepExit==0)&&(HadrShowerContainedInFV==1)&&(PrimaryLepPDG=="
      "13)",
      "GOFF");

  erec_bias_mu_cont_hadrcont->SetLineColor(kMagenta);
  erec_bias_mu_cont_hadrcont->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making erec_bias_mu_exit" << std::endl;
  TH1D* erec_bias_0pi_mu_exit =
      new TH1D("erec_bias_0pi_mu_exit",
               "LepExit_AllHadr;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3])/"
      "nu_4mom[3]) - 1 >> erec_bias_0pi_mu_exit",
      "(stop >= "
      "0)&&(LepExit_AboveThresh==1)&&(PrimaryLepPDG==13)&&IsCC&&Is0Pi",
      "GOFF");

  erec_bias_0pi_mu_exit->Draw();

  TH1D* erec_bias_0pi_mu_cont =
      new TH1D("erec_bias_0pi_mu_cont",
               "LepCont_AllHadr;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto+"
      "LepDepDescendent_FV+LepDepDescendent_veto)/"
      "nu_4mom[3]) - 1 >> erec_bias_0pi_mu_cont",
      "(stop >= "
      "0)&&(LepExit==0)&&(PrimaryLepPDG==13)&&IsCC&&Is0Pi",
      "GOFF");

  erec_bias_0pi_mu_cont->SetLineColor(kRed);
  erec_bias_0pi_mu_cont->Draw("SAME");

  TH1D* erec_bias_0pi_mu_exit_hadrcont =
      new TH1D("erec_bias_0pi_mu_exit_hadrcont",
               "LepExit_HadrCont;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+PrimaryLep_4mom[3])/nu_4mom[3]) - 1 >> "
      "erec_bias_0pi_mu_exit_hadrcont",
      "(stop >= "
      "0)&&(LepExit_AboveThresh==1)&&(HadrShowerContainedInFV==1)&&("
      "PrimaryLepPDG==13)&&IsCC&&Is0Pi",
      "GOFF");

  erec_bias_0pi_mu_exit_hadrcont->SetLineColor(kBlue);
  erec_bias_0pi_mu_exit_hadrcont->Draw("SAME");

  TH1D* erec_bias_0pi_mu_cont_hadrcont =
      new TH1D("erec_bias_0pi_mu_cont_hadrcont",
               "LepCont_HadrCont;E_{rec}/E_{#nu} - 1;Count", 100, -1, 1);

  EDeps->Draw(
      "((TotalNonlep_Dep_FV+LepDep_FV+LepDep_veto+LepDepDescendent_FV+"
      "LepDepDescendent_veto)/nu_4mom[3]) - 1 >> "
      "erec_bias_0pi_mu_cont_hadrcont",
      "(stop >= "
      "0)&&(LepExit==0)&&(HadrShowerContainedInFV==1)&&(PrimaryLepPDG=="
      "13)&&IsCC&&Is0Pi",
      "GOFF");

  erec_bias_0pi_mu_cont_hadrcont->SetLineColor(kMagenta);
  erec_bias_0pi_mu_cont_hadrcont->Draw("SAME");

  c1->BuildLegend(0.2, 0.7, 0.4, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit" << std::endl;
  TH1D* lep_exit = new TH1D("lep_exit", ";LepExitMode;Count", 8, 0, 8);

  EDeps->Draw("LepExitTopology >> lep_exit", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  lep_exit->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_x" << std::endl;
  TH1D* lep_exit_x =
      new TH1D("lep_exit_x", ";LepExit_{x} (cm);Count", 4400, -4000, 400);

  EDeps->Draw("LepExitingPos[0] >> lep_exit_x",
              "(stop >= 0)&&(LepExit_AboveThresh==1)", "GOFF");

  lep_exit_x->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_y" << std::endl;
  TH1D* lep_exit_y =
      new TH1D("lep_exit_y", ";LepExit_{y} (cm);Count", 400, -200, 200);

  EDeps->Draw("LepExitingPos[1] >> lep_exit_y",
              "(stop >= 0)&&(LepExit_AboveThresh==1)", "GOFF");

  lep_exit_y->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_z" << std::endl;
  TH1D* lep_exit_z =
      new TH1D("lep_exit_z", ";LepExit_{z} (cm);Count", 600, -300, 300);

  EDeps->Draw("LepExitingPos[2] >> lep_exit_z",
              "(stop >= 0)&&(LepExit_AboveThresh==1)", "GOFF");

  lep_exit_z->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making true_topo" << std::endl;
  TH1D* true_topo = new TH1D("true_topo", ";True topology;Count", 16, -8, 8);

  EDeps->Draw("Topology >> true_topo", "(stop >= 0)", "GOFF");

  true_topo->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lepdep_FV" << std::endl;
  TH1D* lepdep_FV = new TH1D("lepdep_FV", "FV;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDep_FV >> lepdep_FV", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  std::cout << "Making lepdep_veto" << std::endl;
  TH1D* lepdep_veto =
      new TH1D("lepdep_veto", "Veto;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDep_veto >> lepdep_veto", "(stop >= 0)&&(PrimaryLepPDG==13)",
              "GOFF");

  std::cout << "Making lepdep_timesep_FV" << std::endl;
  TH1D* lepdep_timesep_FV = new TH1D(
      "lepdep_timesep_FV", "FV T > 250 #mus;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDep_timesep_FV >> lepdep_timesep_FV",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making lepdep_timesep_veto" << std::endl;
  TH1D* lepdep_timesep_veto = new TH1D(
      "lepdep_timesep_veto", "Veto T > 250 #mus;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDep_timesep_veto >> lepdep_timesep_veto",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  lepdep_veto->SetLineColor(kRed);
  lepdep_timesep_FV->SetLineColor(kBlue);
  lepdep_timesep_veto->SetLineColor(kMagenta);

  lepdep_FV->Draw();
  lepdep_veto->Draw("SAME");
  lepdep_timesep_FV->Draw("SAME");
  lepdep_timesep_veto->Draw("SAME");
  c1->BuildLegend(0.7, 0.8, 0.95, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lepdepdesc_FV" << std::endl;
  TH1D* lepdepdesc_FV =
      new TH1D("lepdepdesc_FV", "FV;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDepDescendent_FV >> lepdepdesc_FV",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making lepdepdesc_veto" << std::endl;
  TH1D* lepdepdesc_veto =
      new TH1D("lepdepdesc_veto", "Veto;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDepDescendent_veto >> lepdepdesc_veto",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making lepdepdesc_timesep_FV" << std::endl;
  TH1D* lepdepdesc_timesep_FV = new TH1D(
      "lepdepdesc_timesep_FV", "FV T > 250 #mus;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("LepDepDescendent_timesep_FV >> lepdepdesc_timesep_FV",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making lepdepdesc_timesep_veto" << std::endl;
  TH1D* lepdepdesc_timesep_veto =
      new TH1D("lepdepdesc_timesep_veto", "Veto T > 250 #mus;EDep_{#mu};Count",
               100, 0, 5);

  EDeps->Draw("LepDepDescendent_timesep_veto >> lepdepdesc_timesep_veto",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  lepdepdesc_veto->SetLineColor(kRed);
  lepdepdesc_timesep_FV->SetLineColor(kBlue);
  lepdepdesc_timesep_veto->SetLineColor(kMagenta);

  lepdepdesc_FV->Draw();
  lepdepdesc_veto->Draw("SAME");
  lepdepdesc_timesep_FV->Draw("SAME");
  lepdepdesc_timesep_veto->Draw("SAME");
  c1->BuildLegend(0.7, 0.8, 0.95, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making nonlepdep_FV" << std::endl;
  TH1D* nonlepdep_FV =
      new TH1D("nonlepdep_FV", "FV;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("TotalNonlep_Dep_FV >> nonlepdep_FV",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making nonlepdep_veto" << std::endl;
  TH1D* nonlepdep_veto =
      new TH1D("nonlepdep_veto", "Veto;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("TotalNonlep_Dep_veto >> nonlepdep_veto",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making nonlepdep_timesep_FV" << std::endl;
  TH1D* nonlepdep_timesep_FV = new TH1D(
      "nonlepdep_timesep_FV", "FV T > 250 #mus;EDep_{#mu};Count", 100, 0, 5);

  EDeps->Draw("TotalNonlep_Dep_timesep_FV >> nonlepdep_timesep_FV",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  std::cout << "Making nonlepdep_timesep_veto" << std::endl;
  TH1D* nonlepdep_timesep_veto =
      new TH1D("nonlepdep_timesep_veto", "Veto T > 250 #mus;EDep_{#mu};Count",
               100, 0, 5);

  EDeps->Draw("TotalNonlep_Dep_timesep_veto >> nonlepdep_timesep_veto",
              "(stop >= 0)&&(PrimaryLepPDG==13)", "GOFF");

  nonlepdep_veto->SetLineColor(kRed);
  nonlepdep_timesep_FV->SetLineColor(kBlue);
  nonlepdep_timesep_veto->SetLineColor(kMagenta);

  nonlepdep_FV->Draw();
  nonlepdep_veto->Draw("SAME");
  nonlepdep_timesep_FV->Draw("SAME");
  nonlepdep_timesep_veto->Draw("SAME");
  c1->BuildLegend(0.7, 0.8, 0.95, 0.95);

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf]");
  c1->Clear();
}
