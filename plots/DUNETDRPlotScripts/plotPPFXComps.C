void plotPPFXComps() {
  TFile f("FluxErrors_OnAxis_PPFX_AllWeights.root");

  TH1D *Diagnostics_PPFX =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_abs =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_abs/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_att =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_att/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_ttpCpi =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_ttpCpi/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_ttpCk =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_ttpCk/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_ttpCnu =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_ttpCnu/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_ttnua =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_ttnua/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_ttmesinc =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_ttmesinc/FD_nu_numu_stddev"));
  TH1D *Diagnostics_PPFX_oth =
      static_cast<TH1D *>(f.Get("Diagnostics_PPFX_oth/FD_nu_numu_stddev"));

  Diagnostics_PPFX->SetLineColor(kBlack);
  Diagnostics_PPFX_abs->SetLineColor(kAzure + 6);
  Diagnostics_PPFX_att->SetLineColor(kGreen);
  Diagnostics_PPFX_ttpCpi->SetLineColor(kOrange + 3);
  Diagnostics_PPFX_ttpCk->SetLineColor(kMagenta);
  Diagnostics_PPFX_ttnua->SetLineColor(kGreen+3);
  Diagnostics_PPFX_ttmesinc->SetLineColor(kOrange);
  Diagnostics_PPFX_oth->SetLineColor(kRed);

  Diagnostics_PPFX->SetLineWidth(3);
  Diagnostics_PPFX_abs->SetLineWidth(3);
  Diagnostics_PPFX_att->SetLineWidth(3);
  Diagnostics_PPFX_ttpCpi->SetLineWidth(3);
  Diagnostics_PPFX_ttpCk->SetLineWidth(3);
  Diagnostics_PPFX_ttpCnu->SetLineWidth(3);
  Diagnostics_PPFX_ttnua->SetLineWidth(3);
  Diagnostics_PPFX_ttmesinc->SetLineWidth(3);
  Diagnostics_PPFX_oth->SetLineWidth(3);

  Diagnostics_PPFX->SetTitle("Total");
  Diagnostics_PPFX_abs->SetTitle("abs");
  Diagnostics_PPFX_att->SetTitle("att");
  Diagnostics_PPFX_ttpCpi->SetTitle("ttpCpi");
  Diagnostics_PPFX_ttpCk->SetTitle("ttpCk");
  Diagnostics_PPFX_ttpCnu->SetTitle("ttpCnu");
  Diagnostics_PPFX_ttnua->SetTitle("ttnua");
  Diagnostics_PPFX_ttmesinc->SetTitle("ttmesinc");
  Diagnostics_PPFX_oth->SetTitle("oth");

  TCanvas c1("c1","c1");
  Diagnostics_PPFX->Draw("HIST");
  Diagnostics_PPFX->GetYaxis()->SetRangeUser(0,0.3);
  Diagnostics_PPFX_oth->Draw("HISTSAME");
  Diagnostics_PPFX_ttpCpi->Draw("HISTSAME");
  Diagnostics_PPFX_ttpCk->Draw("HISTSAME");
  Diagnostics_PPFX_ttmesinc->Draw("HISTSAME");
  Diagnostics_PPFX_ttnua->Draw("HISTSAME");
  Diagnostics_PPFX_att->Draw("HISTSAME");
  Diagnostics_PPFX_abs->Draw("HISTSAME");

  c1.BuildLegend();

  c1.Print("FD_nu_numu_PPFXAll.pdf");
}
