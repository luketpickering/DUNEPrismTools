#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "TH2Jagged/build/Linux/include/TH2Jagged.h"
#include "TPad.h"

enum spec { kNuMu, kNuE, kNuMuBar, kNuEBar };
enum beam { kNu, kNuBar };

void DumpBinningPlots() {

  gStyle->SetOptStat(false);

  std::map<beam, std::string> filenames = {
      {kNu, "DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_uncert_jagged_opt"},
  };
  std::vector<std::pair<spec, std::string>> specnames = {
      {kNuMu, "ND_nu_ppfx/LBNF_numu_flux_Nom"},
      {kNuMuBar, "ND_nu_ppfx/LBNF_numubar_flux_Nom"},
      {kNuE, "ND_nu_ppfx/LBNF_nue_flux_Nom"}};

  std::map<spec, TH2JaggedD *> hs;

  for (auto bf : filenames) {
    TFile f((bf.second + ".root").c_str());
    for (auto s : specnames) {
      hs.emplace(s.first, dynamic_cast<TH2JaggedD *>(f.Get(s.second.c_str())));
      hs[s.first]->SetDirectory(nullptr);
      hs[s.first]->Scale(0);
      hs[s.first]->fOTitle = "";
    }
  }

  TCanvas c1("c1", "", 400, 600);
  c1.Divide(1, 3);

  TLatex l;
  l.SetTextSize(0.14);
  l.SetTextAlign(32);

  TVirtualPad *p_kNuMu = c1.cd(1);
  p_kNuMu->SetTopMargin(0.15);
  p_kNuMu->SetBottomMargin(0.2);
  p_kNuMu->SetLeftMargin(0.15);
  p_kNuMu->SetRightMargin(0.03);
  TH2Poly *h_kNuMu = hs[kNuMu]->ToTH2Poly();
  h_kNuMu->SetYTitle("Off axis position (m)");
  h_kNuMu->SetYTitle("Off axis position (m)");
  h_kNuMu->SetXTitle("E_{\\nu}\\textrm{ (GeV)}");
  h_kNuMu->GetXaxis()->SetNdivisions(505);
  h_kNuMu->GetYaxis()->SetLabelSize(0.1);
  h_kNuMu->GetYaxis()->SetTitleSize(0.1);
  h_kNuMu->GetYaxis()->SetTitleOffset(0.75);
  h_kNuMu->GetXaxis()->SetLabelSize(0.1);
  h_kNuMu->GetXaxis()->SetTitleSize(0.1);
  h_kNuMu->GetXaxis()->SetTitleOffset(1);
  h_kNuMu->GetXaxis()->SetNdivisions(505);
  h_kNuMu->GetYaxis()->SetNdivisions(505);
  h_kNuMu->Draw();
  l.DrawLatexNDC(0.95, 0.925,
                 "\\textrm{Right sign }\\nu_{\\mu}\\textrm{ binning}");
  TVirtualPad *p_kNuMuBar = c1.cd(2);
  p_kNuMuBar->SetTopMargin(0.15);
  p_kNuMuBar->SetBottomMargin(0.2);
  p_kNuMuBar->SetLeftMargin(0.15);
  p_kNuMuBar->SetRightMargin(0.03);
  TH2Poly *h_kNuMuBar = hs[kNuMuBar]->ToTH2Poly();
  h_kNuMuBar->SetYTitle("Off axis position (m)");
  h_kNuMuBar->SetXTitle("E_{\\nu}\\textrm{ (GeV)}");
  h_kNuMuBar->GetXaxis()->SetNdivisions(505);
  h_kNuMuBar->GetYaxis()->SetLabelSize(0.1);
  h_kNuMuBar->GetYaxis()->SetTitleSize(0.1);
  h_kNuMuBar->GetYaxis()->SetTitleOffset(0.75);
  h_kNuMuBar->GetXaxis()->SetLabelSize(0.1);
  h_kNuMuBar->GetXaxis()->SetTitleSize(0.1);
  h_kNuMuBar->GetXaxis()->SetTitleOffset(1);
  h_kNuMuBar->GetXaxis()->SetNdivisions(505);
  h_kNuMuBar->GetYaxis()->SetNdivisions(505);
  h_kNuMuBar->Draw();
  l.DrawLatexNDC(0.95, 0.925,
                 "\\textrm{Wrong sign }\\nu_{\\mu}\\textrm{ binning}");
  TVirtualPad *p_kNuE = c1.cd(3);
  p_kNuE->SetTopMargin(0.15);
  p_kNuE->SetBottomMargin(0.2);
  p_kNuE->SetLeftMargin(0.15);
  p_kNuE->SetRightMargin(0.03);
  TH2Poly *h_kNuE = hs[kNuE]->ToTH2Poly();
  h_kNuE->SetYTitle("Off axis position (m)");
  h_kNuE->SetYTitle("Off axis position (m)");
  h_kNuE->SetXTitle("E_{\\nu}\\textrm{ (GeV)}");
  h_kNuE->GetXaxis()->SetNdivisions(505);
  h_kNuE->GetYaxis()->SetLabelSize(0.1);
  h_kNuE->GetYaxis()->SetTitleSize(0.1);
  h_kNuE->GetYaxis()->SetTitleOffset(0.75);
  h_kNuE->GetXaxis()->SetLabelSize(0.1);
  h_kNuE->GetXaxis()->SetTitleSize(0.1);
  h_kNuE->GetXaxis()->SetTitleOffset(1);
  h_kNuE->GetXaxis()->SetNdivisions(505);
  h_kNuE->GetYaxis()->SetNdivisions(505);
  h_kNuE->Draw();
  l.DrawLatexNDC(0.95, 0.925, "\\nu_{e}\\textrm{ binning}");

  c1.SaveAs("tikz/Binning.tex");
}
