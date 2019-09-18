TH2 *MatToHist(char const *n, TMatrixD *m) {
  TH2 *matrixHist = new TH2D(n, "", m->GetNcols(), 0, m->GetNcols(),
                             m->GetNrows(), 0, m->GetNrows());

  for (int i = 0; i < m->GetNcols(); ++i) {
    for (int j = 0; j < m->GetNrows(); ++j) {
      int gbin = matrixHist->GetBin(i + 1, j + 1);
      matrixHist->SetBinContent(gbin, (*m)[i][j]);
    }
  }
  return matrixHist;
}

void DumpMatrix() {
  gStyle->SetOptStat(false);

  TFile f("FluxErrors_UncertPoints_Total_onaxis.root");
  TH2 *covmat =
      MatToHist("covmat", static_cast<TMatrixD *>(f.Get("Matrices/covmat")));
  TH2 *corrmat =
      MatToHist("corrmat", static_cast<TMatrixD *>(f.Get("Matrices/corrmat")));

  TCanvas c1("c1", "", 600, 600);

  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.15);
  c1.SetTopMargin(0.15);
  c1.SetBottomMargin(0.15);

  gStyle->SetPalette(kTemperatureMap);

  covmat->GetZaxis()->SetRangeUser(0, 0.0225);
  covmat->GetZaxis()->SetNdivisions(505);
  covmat->GetZaxis()->SetTitleOffset(1.6);
  covmat->GetXaxis()->SetTitle("Bin number");
  covmat->GetYaxis()->SetTitle("Bin number");
  covmat->Draw("COLZ");

  TLatex ltx;
  ltx.SetTextAlign(22);

  ltx.SetTextSize(0.04);
  ltx.SetTextAngle(90);
  ltx.DrawLatexNDC(0.95, 0.9, "C_{ij}");
  ltx.SetTextAngle(0);

  ltx.SetTextSize(0.04);
  ltx.DrawLatexNDC(0.35, 0.95, "Near detector");

  ltx.SetTextSize(0.03);
  ltx.DrawLatexNDC(0.25, 0.91, "\\nu-mode");

  ltx.DrawLatexNDC(0.2, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32, 0.874, "\\bar{\\nu}_{e}");

  ltx.DrawLatexNDC(0.425, 0.91, "\\bar{\\nu}-mode");

  ltx.DrawLatexNDC(0.2 + 0.15, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.15, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.175, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.175, 0.874, "\\bar{\\nu}_{e}");

  ltx.SetTextSize(0.04);
  ltx.DrawLatexNDC(0.35 + 0.35, 0.95, "Far detector");

  ltx.SetTextSize(0.03);
  ltx.DrawLatexNDC(0.25 + 0.35, 0.91, "\\nu-mode");

  ltx.DrawLatexNDC(0.2 + 0.35, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.35, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.35, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.35, 0.874, "\\bar{\\nu}_{e}");

  ltx.DrawLatexNDC(0.425 + 0.35, 0.91, "\\bar{\\nu}-mode");

  ltx.DrawLatexNDC(0.2 + 0.15 + 0.35, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.15 + 0.35, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.175 + 0.35, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.175 + 0.35, 0.874, "\\bar{\\nu}_{e}");

  c1.SaveAs("ErrorMatrices_covmat.pdf");
  corrmat->GetZaxis()->SetRangeUser(-1,1);
  corrmat->GetZaxis()->SetNdivisions(505);
  corrmat->GetZaxis()->SetTitleOffset(1.6);
  corrmat->GetXaxis()->SetTitle("Bin number");
  corrmat->GetYaxis()->SetTitle("Bin number");
  corrmat->Draw("COLZ");

  ltx.SetTextSize(0.04);
  ltx.SetTextAngle(90);
  ltx.DrawLatexNDC(0.95, 0.875, "C_{ij}/(C_{ii}\\times C_{jj})^{1/2}");
  ltx.SetTextAngle(0);

  ltx.SetTextSize(0.04);
  ltx.DrawLatexNDC(0.35, 0.95, "Near detector");

  ltx.SetTextSize(0.03);
  ltx.DrawLatexNDC(0.25, 0.91, "\\nu-mode");

  ltx.DrawLatexNDC(0.2, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32, 0.874, "\\bar{\\nu}_{e}");

  ltx.DrawLatexNDC(0.425, 0.91, "\\bar{\\nu}-mode");

  ltx.DrawLatexNDC(0.2 + 0.15, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.15, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.175, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.175, 0.874, "\\bar{\\nu}_{e}");

  ltx.SetTextSize(0.04);
  ltx.DrawLatexNDC(0.35 + 0.35, 0.95, "Far detector");

  ltx.SetTextSize(0.03);
  ltx.DrawLatexNDC(0.25 + 0.35, 0.91, "\\nu-mode");

  ltx.DrawLatexNDC(0.2 + 0.35, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.35, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.35, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.35, 0.874, "\\bar{\\nu}_{e}");

  ltx.DrawLatexNDC(0.425 + 0.35, 0.91, "\\bar{\\nu}-mode");

  ltx.DrawLatexNDC(0.2 + 0.15 + 0.35, 0.87, "\\nu_{\\mu}");
  ltx.DrawLatexNDC(0.265 + 0.15 + 0.35, 0.874, "\\bar{\\nu}_{\\mu}");
  ltx.DrawLatexNDC(0.295 + 0.175 + 0.35, 0.87, "\\nu_{e}");
  ltx.DrawLatexNDC(0.32 + 0.175 + 0.35, 0.874, "\\bar{\\nu}_{e}");

  c1.SaveAs("ErrorMatrices_corrmat.pdf");
}
