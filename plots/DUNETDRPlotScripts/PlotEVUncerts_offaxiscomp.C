#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "TH2Jagged/build/Linux/include/TH2Jagged.h"

#include "colordef.h"

#include "TPad.h"

enum spec { kNuMu, kNuE, kNuMuBar, kNuEBar };
enum beam { kNu, kNuBar };

void PlotEVUncerts_offaxiscomp() {

  DeclareColors();

  gStyle->SetOptStat(false);
  gStyle->SetPadGridX(true);

  std::map<int, int> Colors;
  Colors[0] = kBlack;
  Colors[1] = kDUNEOrange;
  Colors[2] = kDUNEBlue;
  Colors[3] = kMSUGreen;
  Colors[4] = kMSUPurple;
  Colors[5] = kGray + 2;

  std::string filename_oa = "FluxErrors_UncertPoints_Total.root";
  std::string filename = "FluxErrors_UncertPoints_Total_onaxis_onebin.root";
  std::vector<std::pair<spec, std::string>> specnames = {{kNuMu, "numu"},
                                                         {kNuMuBar, "numubar"},
                                                         {kNuE, "nue"},
                                                         {kNuEBar, "nuebar"}};

  std::vector<std::pair<beam, std::string>> beamnames = {{kNu, "ND_nu"},
                                                         {kNuBar, "ND_nubar"}};

  std::map<spec, double> specscale = {
      {kNuMu, 0.26}, {kNuE, 0.19}, {kNuMuBar, 0.26}, {kNuEBar, 0.19}};

  std::map<spec, std::string> specyscalestr = {{kNuMu, "10^{-9}"},
                                               {kNuE, "10^{-11}"},
                                               {kNuMuBar, "10^{-9}"},
                                               {kNuEBar, "10^{-11}"}};

  std::map<spec, std::string> speclatexnames = {
      {kNuMu, "\\nu_{\\mu}"},
      {kNuE, "\\nu_{e}"},
      {kNuMuBar, "\\bar{\\nu}_{\\mu}"},
      {kNuEBar, "\\bar{\\nu}_{e}"}};

  std::map<int, std::map<beam, std::map<spec, TH1 *>>> Hists;
  std::map<int, std::map<beam, std::map<spec, TH1 *>>> Hists_oa;

  TFile f(filename.c_str());
  TFile f_oa(filename_oa.c_str());
  for (int i = 0; i < 6; ++i) {
    for (auto bf : beamnames) {
      for (auto s : specnames) {

        std::string pname;
        if (i == 0) {
          pname = "TotalError1D_total/" + bf.second + "_" + s.second;
        }

        else {
          pname = "FluxParameters/param_pca_" + std::to_string(i - 1) + "/" +
                  bf.second + "_" + s.second;
        }
        TH1F *hist = nullptr;
        f.GetObject(pname.c_str(), hist);
        if (!hist) {
          std::cout << "[ERROR] Failed to read: " << pname << std::endl;
          abort();
        }
        Hists[i][bf.first][s.first] = hist;

        TH2JaggedF *hist_oa = nullptr;
        f_oa.GetObject(pname.c_str(), hist_oa);
        if (!hist) {
          std::cout << "[ERROR] Failed to read: " << pname << std::endl;
          abort();
        }
        Hists_oa[i][bf.first][s.first] = hist_oa->NonUniformSlice(3);
        Hists_oa[i][bf.first][s.first]->SetName((pname + "_s1").c_str());
        ;
      }
    }
  }

  TCanvas c1("", "", 500, 700);
  c1.Clear();
  c1.cd();

  std::vector<TVirtualPad *> pads;

  double XLL = 0;
  double XLR = 0.5;
  double XRL = 0.5;
  double XRR = 1;
  double YMin = 0;
  double NomPadHeight = 0.23;
  double PadMax = (4 * NomPadHeight);

  double BottomMarg = 0.1;
  double NotBottomMarg = 0.045;
  double TopMargin = 0.02;
  double RightMargin = 0.04;

  // 3*PHB*BottomMarg + PHNB*NotBottomMarg + 3*PHB*(1-BottomMarg) + PHNB*(1
  // -NotBottomMarg) = 4*NomPadHeight
  // -- Want same draw area
  // PHB*(1-BottomMarg) = PHNB*(1-NotBottomMarg)
  // PHB = PHNB*(1-NotBottomMarg)/(1-BottomMarg)
  double OMMargRat = (1.0 - NotBottomMarg) / (1.0 - BottomMarg);
  // PHNB * (3*OMMargRat*BottomMarg + NotBottomMarg +
  // 3*OMMargRat*(1-BottomMarg)
  // + (1 - NotBottomMarg)) = 4*NomPadHeight
  double PadHeightNotBottom =
      (4 * NomPadHeight) /
      (3 * OMMargRat * BottomMarg + NotBottomMarg +
       3 * OMMargRat * (1 - BottomMarg) + (1 - NotBottomMarg));
  double PadHeightBottom = PadHeightNotBottom * OMMargRat;

  TLegend l(0.1, PadMax + 0.04, 1, 1);
  l.SetBorderSize(0);
  l.SetNColumns(3);

  size_t p_id = 0;
  for (auto s : specnames) {
    for (auto bf : beamnames) {
      c1.cd();
      bool isLeftCol = !(p_id & 1);
      int NRow = (p_id >> 1);
      std::string n = std::to_string(isLeftCol) + s.second;
      double YMax;
      double YMin;
      if (p_id > 5) {
        YMax = PadMax - 3 * PadHeightNotBottom;
        YMin = PadMax - 3 * PadHeightNotBottom - PadHeightBottom;
      } else {
        YMax = PadMax - (NRow + 1) * PadHeightNotBottom;
        YMin = PadMax - NRow * PadHeightNotBottom;
      }

      TVirtualPad *p = new TPad(n.c_str(), "", isLeftCol ? XLL : XRL, YMin,
                                isLeftCol ? XLR : XRR, YMax);
      p->AppendPad();
      p->cd();
      pads.push_back(p);
      p->SetTopMargin(TopMargin);
      p->SetLeftMargin(0.2);
      p->SetRightMargin(RightMargin);
      if (p_id > 5) {
        p->SetBottomMargin(BottomMarg);
      } else {
        p->SetBottomMargin(NotBottomMarg);
      }

      bool drawn = false;
      for (int i = 0; i < 6; ++i) {

        TH1 *h = Hists[i][bf.first][s.first];
        TH1 *h_oa = Hists_oa[i][bf.first][s.first];

        if (!p_id) {
          if (i == 0) {
            l.AddEntry(h, "Total Error", "l");
          } else {
            l.AddEntry(h, ("Component " + std::to_string(i - 1)).c_str(), "l");
          }
        }

        h->SetLineColor(Colors[i]);
        h->SetLineWidth(2);

        h_oa->SetLineColor(Colors[i]);
        h_oa->SetLineWidth(2);
        h_oa->SetLineStyle(2);

        if (!drawn) {
          std::string ytitle = "Relative variation";
          std::string xtitle = "E_{\\nu}\\textrm{ (GeV)}";

          h->GetYaxis()->SetNdivisions(505);
          h->GetXaxis()->SetNdivisions(508);
          h->GetYaxis()->SetRangeUser(0, specscale[s.first]);
          h->GetXaxis()->SetRangeUser(0, 10);

          double ts = 0.125;
          double ls = 0.125;
          if (p_id > 5) {
            h->GetYaxis()->SetLabelSize(ls *
                                        (PadHeightNotBottom / PadHeightBottom));
            h->GetYaxis()->SetLabelOffset(
                h->GetYaxis()->GetLabelOffset() /
                (PadHeightNotBottom / PadHeightBottom));

            h->GetXaxis()->SetLabelSize(ls);
          } else {
            h->GetYaxis()->SetLabelSize(ls);
            h->GetXaxis()->SetLabelSize(0);
          }
          h->Draw("HIST");

          // Manually doing axis titles because ROOT auto-movement is a
          // massive pain
          TLatex ltx;
          ltx.SetTextAngle(90);
          double YTS = 0.1;
          if (p_id > 5) {
            ltx.SetTextSize(YTS / (PadHeightNotBottom / PadHeightBottom));
          } else {
            ltx.SetTextSize(YTS);
          }
          ltx.DrawLatexNDC(0.05, 0.2, ytitle.c_str());
          ltx.SetTextAngle(0);

          if (p_id > 5) {
            ltx.SetTextSize(0.125);
            ltx.DrawLatexNDC(0.675, -0.1, xtitle.c_str());
          }
          ltx.SetTextSize(0.3);
          ltx.SetTextAlign(22);
          ltx.DrawLatexNDC(0.885, 0.85, speclatexnames[s.first].c_str());
        } else {
          h->Draw("HISTSAME");
        }
        h_oa->Draw("HISTSAME");
        drawn = true;
      }
      p_id++;
    }
  }

  for (auto p : pads) {
    p->RedrawAxis("g");
  }

  c1.cd();
  l.SetTextSize(0.03);
  l.Draw();
  TLatex ltx;
  ltx.SetTextSize(0.075);
  ltx.SetTextAlign(22);
  ltx.DrawLatexNDC(0.2785, 0.935, "\\nu\\textrm{-mode}");
  ltx.DrawLatexNDC(0.7785, 0.935, "\\bar{\\nu}\\textrm{-mode}");
  ltx.SetTextSize(0.04);
  c1.SaveAs("tikz/EvUncerts_oacomp.tex");
}
