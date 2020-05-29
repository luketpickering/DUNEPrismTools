#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "TH2Jagged/build/Linux/include/TH2Jagged.h"

#include "colordef.h"

#include "TPad.h"

enum spec { kNuMu, kNuE, kNuMuBar, kNuEBar };
enum beam { kNu, kNuBar };

void PlotEVUncertsOffAxis() {

  DeclareColors();

  gStyle->SetOptStat(false);
  gStyle->SetPadGridX(true);
  gStyle->SetPalette(kTemperatureMap);

  std::string filename = "FluxErrors_UncertPoints_Total_allpca.root";
  std::vector<std::pair<spec, std::string>> specnames = {{kNuMu, "numu"},
                                                         {kNuMuBar, "numubar"},
                                                         {kNuE, "nue"},
                                                         {kNuEBar, "nuebar"}};

  std::vector<std::pair<beam, std::string>> beamnames = {{kNu, "ND_nu"},
                                                         {kNuBar, "ND_nubar"}};

  std::map<spec, double> specscale = {
      {kNuMu, 0.22}, {kNuE, 0.15}, {kNuMuBar, 0.22}, {kNuEBar, 0.15}};

  std::map<spec, std::string> speclatexnames = {
      {kNuMu, "\\nu_{\\mu}"},
      {kNuE, "\\nu_{e}"},
      {kNuMuBar, "\\bar{\\nu}_{\\mu}"},
      {kNuEBar, "\\bar{\\nu}_{e}"}};

  std::map<int, std::map<beam, std::map<spec, TH2JaggedF *>>> Hists;


  TFile f(filename.c_str());
  for (int i = 0; i < 6; ++i) {
    for (auto bf : beamnames) {
      for (auto s : specnames) {

        std::string pname;
        if (i == 0) {
          pname = "TotalError1D_pcaonly/" + bf.second + "_" + s.second;
        }

        else {
          pname = "FluxParameters/param_pca_" + std::to_string(i - 1) + "/" +
                  bf.second + "_" + s.second;
        }
        TH2JaggedF *hist = nullptr;
        f.GetObject(pname.c_str(), hist);
        if (!hist) {
          std::cout << "[ERROR] Failed to read: " << pname << std::endl;
          abort();
        }
        Hists[i][bf.first][s.first] = hist;
      }
    }
  }

  TCanvas c1("", "", 500, 700);
  c1.Clear();
  c1.cd();

  for (int i = 0; i < 6; ++i) {
    std::vector<TVirtualPad *> pads;

    double XLL = 0;
    double XLR = 0.5;
    double XRL = 0.5;
    double XRR = 1;
    double YMin = 0;
    double NomPadHeight = 0.25;
    double PadMax = (4 * NomPadHeight);

    double BottomMarg = 0.1;
    double NotBottomMarg = 0.045;
    double TopMargin = 0.02;
    double RightMargin = 0.25;

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

        TH1 *h = Hists[i][bf.first][s.first]->ToUniformTH2("width");

        std::string ytitle = "Off axis position (m)";
        std::string xtitle = "E_{#nu} (GeV)";

        h->GetYaxis()->SetNdivisions(505);
        h->GetXaxis()->SetNdivisions(508);
        h->GetXaxis()->SetRangeUser(0, 10);

        double ts = 0.125;
        double ls = 0.125;
        if (p_id > 5) {
          h->GetYaxis()->SetLabelSize(ls *
                                      (PadHeightNotBottom / PadHeightBottom));
          h->GetYaxis()->SetLabelOffset(h->GetYaxis()->GetLabelOffset() /
                                        (PadHeightNotBottom / PadHeightBottom));

          h->GetXaxis()->SetLabelSize(ls);
          h->GetZaxis()->SetLabelSize(ls*.6);
          h->GetZaxis()->SetNdivisions(505);
        } else {
          h->GetZaxis()->SetLabelSize(ls*.6);
          h->GetZaxis()->SetNdivisions(505);
          h->GetYaxis()->SetLabelSize(ls);
          h->GetXaxis()->SetLabelSize(0);
        }
        h->SetMaximum(specscale[s.first] * (i == 0 ? 1 : 0.8));
        h->SetMinimum(-specscale[s.first] * (i == 0 ? 1 : 0.8));
        h->Draw("COLZ");

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
        ltx.DrawLatexNDC(0.075, 0.2, ytitle.c_str());
        ltx.SetTextAngle(270);
        if(i != 0){
                ltx.DrawLatexNDC(0.925, 0.9, "Relative variation");} else {
                  ltx.DrawLatexNDC(0.925, 0.9, "Fracional error");
                }
        ltx.SetTextAngle(0);

        ltx.SetTextSize(0.18);
        ltx.SetTextAlign(22);
        ltx.DrawLatexNDC(0.7, 0.85, speclatexnames[s.first].c_str());

        p_id++;

        if (p_id > 5) {
          c1.cd();
          ltx.SetTextSize(0.03);
          ltx.DrawLatexNDC(0.35, 0.015, xtitle.c_str());
          ltx.DrawLatexNDC(0.85, 0.015, xtitle.c_str());
        }
      }
    }

    for (auto p : pads) {
      p->RedrawAxis("g");
    }

    if (i == 0) {
      c1.SaveAs("EvUncerts_offaxis_total.pdf");
    } else {
      c1.SaveAs(
          ("EvUncerts_offaxis_component_" + std::to_string(i - 1) + ".pdf")
              .c_str());
    }
  }
}
