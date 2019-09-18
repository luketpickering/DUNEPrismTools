#pragma cling load("../../build/Linux/lib/libTH2Jagged.so");

void Dump1D() {

  TFile *f = new TFile(
      "DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_uncert_jagged_opt.root");

  TH2 *nom_ND = static_cast<TH2 *>(f->Get("ND_nu_ppfx/LBNF_numu_flux_Nom"));
  TH2 *WL_ND = static_cast<TH2 *>(f->Get("ND_nu_WL_p1/LBNF_numu_flux"));
  TH2 *HCp1_ND = static_cast<TH2 *>(f->Get("ND_nu_HC_m1/LBNF_numu_flux"));
  TH2 *HCm1_ND = static_cast<TH2 *>(f->Get("ND_nu_HC_p1/LBNF_numu_flux"));
  TH2 *TargetDensity_m1_ND =
      static_cast<TH2 *>(f->Get("ND_nu_TargetDensity_m1/LBNF_numu_flux"));

  TH1D *nom_ND_1D = nom_ND->ProjectionX("nom_ND", 1, 1);
  TH1D *WL_ND_1D = WL_ND->ProjectionX("WL_ND", 1, 1);
  TH1D *HCp1_ND_1D = HCp1_ND->ProjectionX("HCp1_ND", 1, 1);
  TH1D *HCm1_ND_1D = HCm1_ND->ProjectionX("HCm1_ND", 1, 1);
  TH1D *TargetDensity_m1_ND_1D =
      TargetDensity_m1_ND->ProjectionX("TargetDensity_m1_ND", 1, 1);

  TH1 *nom_FD = static_cast<TH1 *>(f->Get("FD_nu_ppfx/LBNF_numu_flux_Nom"));
  TH1 *WL_FD = static_cast<TH1 *>(f->Get("FD_nu_WL_p1/LBNF_numu_flux"));
  TH1 *HCp1_FD = static_cast<TH1 *>(f->Get("FD_nu_HC_m1/LBNF_numu_flux"));
  TH1 *HCm1_FD = static_cast<TH1 *>(f->Get("FD_nu_HC_p1/LBNF_numu_flux"));
  TH1 *TargetDensity_m1_FD =
      static_cast<TH1 *>(f->Get("FD_nu_TargetDensity_m1/LBNF_numu_flux"));

  nom_ND_1D->SetDirectory(nullptr);
  WL_ND_1D->SetDirectory(nullptr);
  HCp1_ND_1D->SetDirectory(nullptr);
  HCm1_ND_1D->SetDirectory(nullptr);
  TargetDensity_m1_ND_1D->SetDirectory(nullptr);

  nom_FD->SetDirectory(nullptr);
  WL_FD->SetDirectory(nullptr);
  HCp1_FD->SetDirectory(nullptr);
  HCm1_FD->SetDirectory(nullptr);
  TargetDensity_m1_FD->SetDirectory(nullptr);

  TFile *f1 = new TFile("1DFluxes.root", "RECREATE");

  f1->WriteTObject(nom_ND_1D, "nom_ND_1D");
  f1->WriteTObject(WL_ND_1D, "WL_ND_1D");
  f1->WriteTObject(HCp1_ND_1D, "HCp1_ND_1D");
  f1->WriteTObject(HCm1_ND_1D, "HCm1_ND_1D");
  f1->WriteTObject(TargetDensity_m1_ND_1D, "TargetDensity_m1_ND_1D");

  WL_ND_1D->Divide(nom_ND_1D);
  HCp1_ND_1D->Divide(nom_ND_1D);
  HCm1_ND_1D->Divide(nom_ND_1D);
  TargetDensity_m1_ND_1D->Divide(nom_ND_1D);

  f1->WriteTObject(WL_ND_1D, "WL_ND_1D_twk");
  f1->WriteTObject(HCp1_ND_1D, "HCp1_ND_1D_twk");
  f1->WriteTObject(HCm1_ND_1D, "HCm1_ND_1D_twk");
  f1->WriteTObject(TargetDensity_m1_ND_1D, "TargetDensity_m1_ND_1D_twk");

  f1->WriteTObject(nom_FD, "nom_FD");
  f1->WriteTObject(WL_FD, "WL_FD");
  f1->WriteTObject(HCp1_FD, "HCp1_FD");
  f1->WriteTObject(HCm1_FD, "HCm1_FD");
  f1->WriteTObject(TargetDensity_m1_FD, "TargetDensity_m1_FD");

  WL_FD->Divide(nom_FD);
  HCp1_FD->Divide(nom_FD);
  HCm1_FD->Divide(nom_FD);
  TargetDensity_m1_FD->Divide(nom_FD);

  f1->WriteTObject(nom_FD, "nom_FD_twk");
  f1->WriteTObject(WL_FD, "WL_FD_twk");
  f1->WriteTObject(HCp1_FD, "HCp1_FD_twk");
  f1->WriteTObject(HCm1_FD, "HCm1_FD_twk");
  f1->WriteTObject(TargetDensity_m1_FD, "TargetDensity_m1_FD_twk");


  f1->Close();
}
