#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "TH2Jagged/build/Linux/include/TH2Jagged.h"

#include "colordef.h"

#include "TPad.h"

enum spec { kNuMu, kNuE, kNuMuBar, kNuEBar };
enum beam { kNu, kNuBar, kNu280, kNuBar280 };

void HornCurrentFluxDiffPred() {

  DeclareColors();

  gStyle->SetOptStat(false);
  gStyle->SetPadGridX(true);

  std::map<int, int> Colors;
  Colors[293] = kDUNEOrange;
  Colors[280] = kDUNEBlue;

  std::string filename = "DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_"
                         "uncertbin_offaxis_280kAOnAxis.root";
  std::vector<std::pair<spec, std::string>> specnames = {
      {kNuMu, "LBNF_numu_flux"},
      {kNuMuBar, "LBNF_numubar_flux"},
      {kNuE, "LBNF_nue_flux"},
      {kNuEBar, "LBNF_nuebar_flux"}};
  std::vector<std::pair<beam, std::string>> beamnames = {
      {kNu, "ND_nu"},
      {kNuBar, "ND_nubar"},
      {kNu280, "ND_nu_280kA"},
      {kNuBar280, "ND_nubar_280kA"}};

  std::map<spec, double> specscale = {
      {kNuMu, 1E9}, {kNuE, 1E11}, {kNuMuBar, 1E9}, {kNuEBar, 1E11}};

  std::map<spec, std::string> specyscalestr = {{kNuMu, "10^{-9}"},
                                               {kNuE, "10^{-11}"},
                                               {kNuMuBar, "10^{-9}"},
                                               {kNuEBar, "10^{-11}"}};

  std::map<spec, std::string> speclatexnames = {
      {kNuMu, "\\nu_{\\mu}"},
      {kNuE, "\\nu_{e}"},
      {kNuMuBar, "\\bar{\\nu}_{\\mu}"},
      {kNuEBar, "\\bar{\\nu}_{e}"}};

        std::map<int, std::string> currlatexnames = {
      {293, "\\textrm{Nominal}, \\pm293\\textrm{ kA}"},
      {280, "\\textrm{Varied}, \\pm280\\textrm{ kA}"},};

  std::map<beam, std::map<spec, TH1 *>> Hists;

  TFile f(filename.c_str());
  for (auto bf : beamnames) {
    for (auto s : specnames) {

      std::string pname = bf.second + "_ppfx/" + s.second + "_Nom";
      TH2JaggedD *hist = nullptr;
      f.GetObject(pname.c_str(), hist);
      if (!hist) {
        std::cout << "[ERROR] Failed to read: " << pname << std::endl;
        abort();
      }
      TH1 *proj =
          hist->UniformRange((pname + "_s1").c_str(), 2, 6, true)->NonUniformSlice(1);
      proj->Scale(specscale[s.first]);
      proj->SetDirectory(nullptr);
      Hists[bf.first][s.first] = proj;
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

  TLegend l(0.1, PadMax + 0.02, 1, 1);
  l.SetBorderSize(0);
  l.SetNColumns(3);
  l.SetTextSize(0.03);

  size_t p_id = 0;
  for (auto s : specnames) {
    for (bool isLeftCol : {true, false}) {
      c1.cd();
      int NRow = (p_id >> 1);
      std::string n =  std::to_string(isLeftCol) + s.second;
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
      for (int curr : {293, 280}) {
        auto bf = kNu;
        if (isLeftCol) {
          bf = (curr == 293) ? kNu : kNu280;
        } else {
          bf = (curr == 293) ? kNuBar : kNuBar280;
        }

        TH1 *h = Hists[bf][s.first];

        if (!p_id) {
          l.AddEntry(h, currlatexnames[curr].c_str(), "l");
        }

        h->SetLineColor(Colors[curr]);
        h->SetLineWidth(2);
        if (!drawn) {
          std::string ytitle =
              speclatexnames[s.first] +
              " /\\textrm{cm}^{2}/\\textrm{POT per 1 GeV }#times " +
              specyscalestr[s.first];
          std::string xtitle = "E_{\\nu}\\textrm{ (GeV)}";

          h->GetYaxis()->SetNdivisions(505);
          h->GetXaxis()->SetNdivisions(508);
          h->GetXaxis()->SetRangeUser(0,10);

          double ts = 0.08;
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
          double YTS = 0.08;
          if (p_id > 5) {
            ltx.SetTextSize(YTS / (PadHeightNotBottom / PadHeightBottom));
          } else {
            ltx.SetTextSize(YTS);
          }
          ltx.DrawLatexNDC(0.05, 0.05, ytitle.c_str());
          ltx.SetTextAngle(0);

          if (p_id > 5) {
            ltx.SetTextSize(0.125);
            ltx.DrawLatexNDC(0.675, -0.1, xtitle.c_str());
          }
          ltx.SetTextSize(0.3);
          ltx.SetTextAlign(22);
          ltx.DrawLatexNDC(0.86, 0.85, speclatexnames[s.first].c_str());
        } else {
          h->Draw("HISTSAME");
        }
        drawn = true;
      }
      p_id++;
    }
  }

  for (auto p : pads) {
    p->RedrawAxis("g");
  }

  c1.cd();
  l.SetTextSize(0.05);
  l.Draw();
  TLatex ltx;
  ltx.SetTextSize(0.05);
  ltx.SetTextAlign(22);
  ltx.DrawLatexNDC(0.2785, 0.935, "\\nu\\textrm{-mode}");
  ltx.DrawLatexNDC(0.7785, 0.935, "\\bar{\\nu}\\textrm{-mode}");
  ltx.SetTextSize(0.04);
  c1.SaveAs((std::string("tikz/") + "HornCurrentFluxPredictions.tex")
                .c_str());
}
