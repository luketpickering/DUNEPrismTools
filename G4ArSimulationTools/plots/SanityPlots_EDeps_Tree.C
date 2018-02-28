{
  TCanvas* c1 = new TCanvas("c1", "");

  std::cout << "Making enu" << std::endl;
  TH1D* enu = new TH1D("enu", ";E_{#nu} (GeV);Count", 250, 0, 25);

  EDeps->Draw("nu_4mom[3] >> enu", "", "GOFF");

  enu->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf[");
  c1->Clear();

  std::cout << "Making nu_pdg" << std::endl;
  TH1D* nu_pdg = new TH1D("nu_pdg", ";PDG_{#nu};Count", 6, 11, 17);

  EDeps->Draw("nu_PDG >> nu_pdg", "", "GOFF");

  std::cout << "Making nu_pdg_neg" << std::endl;
  TH1D* nu_pdg_neg = new TH1D("nu_pdg_neg", ";PDG_{#nu};Count", 6, 11, 17);

  EDeps->Draw("abs(nu_PDG) >> nu_pdg_neg", "nu_PDG<0", "GOFF");

  nu_pdg->Draw();
  nu_pdg_neg->SetLineColor(kRed);
  nu_pdg_neg->Draw("SAME");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making vtx_pos_xy" << std::endl;
  TH2D* vtx_pos_xy = new TH2D("vtx_pos_xy", ";VtxPos_{x} (cm);VtxPos_{y} (cm)",
                              4400, -4000, 400, 400, -200, 200);

  EDeps->Draw("vtx[1]:vtx[0] >> vtx_pos_xy", "", "GOFF");

  vtx_pos_xy->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making vtx_pos_yz" << std::endl;
  TH2D* vtx_pos_yz = new TH2D("vtx_pos_yz", ";VtxPos_{z} (cm);VtxPos_{y} (cm)",
                              600, -300, 300, 400, -200, 200);

  EDeps->Draw("vtx[1]:vtx[2] >> vtx_pos_yz", "", "GOFF");

  vtx_pos_yz->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making q2" << std::endl;
  TH1D* q2 = new TH1D("q2", ";Q^{2} (GeV^{2});Count", 100, 0, 10);

  EDeps->Draw("Q2_True >> q2", "", "GOFF");

  q2->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making y" << std::endl;
  TH1D* y = new TH1D("y", ";Elasticity;Count", 100, 0, 1);

  EDeps->Draw("y_True >> y", "", "GOFF");

  y->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making wrest" << std::endl;
  TH1D* wrest = new TH1D("wrest", ";W_{rest} (GeV);Count", 100, 0.8, 3);

  EDeps->Draw("W_Rest >> wrest", "", "GOFF");

  wrest->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NLep" << std::endl;
  TH1D* NLep = new TH1D("NLep", ";N_{#ell};Count", 20, 0, 20);

  EDeps->Draw("NLep >> NLep", "", "GOFF");

  NLep->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NPi0" << std::endl;
  TH1D* NPi0 = new TH1D("NPi0", ";N_{#pi^{0}};Count", 20, 0, 20);

  EDeps->Draw("NPi0 >> NPi0", "", "GOFF");

  NPi0->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NPiC" << std::endl;
  TH1D* NPiC = new TH1D("NPiC", ";N_{#pi^{#pm}};Count", 20, 0, 20);

  EDeps->Draw("NPiC >> NPiC", "", "GOFF");

  NPiC->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NProton" << std::endl;
  TH1D* NProton = new TH1D("NProton", ";N_{proton};Count", 20, 0, 20);

  EDeps->Draw("NProton >> NProton", "", "GOFF");

  NProton->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NNeutron" << std::endl;
  TH1D* NNeutron = new TH1D("NNeutron", ";N_{neutron};Count", 20, 0, 20);

  EDeps->Draw("NNeutron >> NNeutron", "", "GOFF");

  NNeutron->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making NGamma" << std::endl;
  TH1D* NGamma = new TH1D("NGamma", ";N_{#gamma};Count", 20, 0, 20);

  EDeps->Draw("NGamma >> NGamma", "", "GOFF");

  NGamma->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making primlep_e" << std::endl;
  TH2D* primlep_e = new TH2D(
      "primlep_e", ";EDep_{#ell} (GeV);T_{#ell} + N_{#ell}*M_{#ell} (GeV)",
      100, 0, 10, 100, 0, 10);

  EDeps->Draw("LepDep_FV+LepDep_veto:PrimaryLep_4mom[3] >> primlep_e", "",
              "GOFF");

  primlep_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making pi0_e" << std::endl;
  TH2D* pi0_e = new TH2D(
      "pi0_e",
      ";EDep_{#pi^{0}} (GeV);T_{#pi^{0}} + N_{#pi^{0}}*M_{#pi^{0}} (GeV)", 100,
      0, 10, 100, 0, 10);

  EDeps->Draw("Pi0Dep_FV+Pi0Dep_veto:EKinPi0_True+EMassPi0_True >> pi0_e", "",
              "GOFF");

  pi0_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making pic_e" << std::endl;
  TH2D* pic_e = new TH2D("pic_e",
                         ";EDep_{#pi^{#pm}} (GeV);T_{#pi^{#pm}} + N_{#pi^{#pm}}*M_{#pi^{#pm}} (GeV)",
                         100, 0, 10, 100, 0, 10);

  EDeps->Draw("PiCDep_FV+PiCDep_veto:EKinPiC_True+EMassPiC_True >> pic_e", "",
              "GOFF");

  pic_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making neutron_e" << std::endl;
  TH2D* neutron_e = new TH2D(
      "neutron_e",
      ";EDep_{neutron} (GeV);T_{neutron} + N_{neutron}*M_{neutron} (GeV)", 100,
      0, 10, 100, 0, 10);

  EDeps->Draw(
      "NeutronDep_FV+NeutronDep_veto:EKinNeutron_True+EMassNeutron_True >> "
      "neutron_e",
      "", "GOFF");

  neutron_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making proton_e" << std::endl;
  TH2D* proton_e =
      new TH2D("proton_e",
               ";EDep_{proton} (GeV);T_{proton} + N_{proton}*M_{proton} (GeV)",
               100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "ProtonDep_FV+ProtonDep_veto:EKinProton_True+EMassProton_True >> "
      "proton_e",
      "", "GOFF");

  proton_e->Draw("COLZ");

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_dep" << std::endl;
  TH2D* enu_vs_dep = new TH2D("enu_vs_dep", ";E_{#nu} (GeV); EDep_{det} (GeV)",
                              100, 0, 10, 100, 0, 10);

  EDeps->Draw(
      "TotalNonlep_Dep_FV+TotalNonlep_Dep_veto+LepDep_FV+LepDep_veto:nu_4mom[3] >> enu_vs_dep",
      "", "GOFF");

  enu_vs_dep->Draw("COLZ");


  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making enu_vs_dep_FV" << std::endl;
  TH2D* enu_vs_dep_FV =
      new TH2D("enu_vs_dep_FV", ";E_{#nu} (GeV); EDep_{FV} (GeV)", 100, 0, 10,
               100, 0, 10);

  EDeps->Draw("TotalNonlep_Dep_FV+LepDep_FV:nu_4mom[3] >> enu_vs_dep_FV", "",
              "GOFF");

  enu_vs_dep_FV->Draw("COLZ");


  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit" << std::endl;
  TH1D* lep_exit = new TH1D("lep_exit", ";LepExitMode;Count", 8, 0, 8);

  EDeps->Draw(
      "LepExitTopology >> lep_exit",
      "", "GOFF");

  lep_exit->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_x" << std::endl;
  TH1D* lep_exit_x = new TH1D("lep_exit_x", ";LepExit_{x} (cm);Count", 4400,-4000,400);

  EDeps->Draw(
      "LepExitingPos[0] >> lep_exit_x",
      "", "GOFF");

  lep_exit_x->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_y" << std::endl;
  TH1D* lep_exit_y = new TH1D("lep_exit_y", ";LepExit_{y} (cm);Count", 400,-200,200);

  EDeps->Draw(
      "LepExitingPos[1] >> lep_exit_y",
      "", "GOFF");

  lep_exit_y->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making lep_exit_z" << std::endl;
  TH1D* lep_exit_z = new TH1D("lep_exit_z", ";LepExit_{z} (cm);Count", 600,-300,300);

  EDeps->Draw(
      "LepExitingPos[2] >> lep_exit_z",
      "", "GOFF");

  lep_exit_z->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  std::cout << "Making true_topo" << std::endl;
  TH1D* true_topo = new TH1D("true_topo", ";True topology;Count", 16, -8, 8);

  EDeps->Draw(
      "Topology >> true_topo",
      "", "GOFF");

  true_topo->Draw();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf");
  c1->Clear();

  c1->SaveAs("SanityPlots_EDeps_Tree.pdf]");
  c1->Clear();
}
