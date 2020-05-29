#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "colordef.h"

void PlotEVDistrib() {

  DeclareColors();

  gStyle->SetOptStat(false);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);

  TFile f("FluxErrors_UncertPoints_Total_onaxis_allpca.root");

  TH1 *ev = nullptr;
  f.GetObject("pca_eigenvalues", ev);

  TCanvas c1("c1", "", 600, 600);

  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.15);
  c1.SetTopMargin(0.15);
  c1.SetBottomMargin(0.15);
  c1.SetLogy();
  ev->GetXaxis()->SetNdivisions(505);
  ev->GetYaxis()->SetNdivisions(505);

  ev->SetLineWidth(2);
  ev->GetYaxis()->SetTitle("\\lambda_{\\it{i}}\\ \\textrm{Magnitude}");
  ev->GetXaxis()->SetTitle("\\textrm{Spectral component number, }\\it{i}");
  ev->Draw();

  c1.SaveAs("tikz/allpca_evs.tex");
}