#pragma cling load("TH2Jagged/build/Linux/lib/libTH2Jagged.so");

#include "colordef.h"

#include "TPad.h"

enum spec { kNuMu, kNuE, kNuMuBar, kNuEBar };
std::vector<std::pair<spec, std::string>> specnames = {
    {kNuMu, "numu"}, {kNuMuBar, "numubar"}, {kNuE, "nue"}, {kNuEBar, "nuebar"}};
std::map<spec, std::string> specyscalestr = {{kNuMu, "10^{-9}"},
                                             {kNuE, "10^{-11}"},
                                             {kNuMuBar, "10^{-9}"},
                                             {kNuEBar, "10^{-11}"}};

std::map<spec, std::string> speclatexnames = {{kNuMu, "\\nu_{\\mu}"},
                                              {kNuE, "\\nu_{e}"},
                                              {kNuMuBar, "\\bar{\\nu}_{\\mu}"},
                                              {kNuEBar, "\\bar{\\nu}_{e}"}};

enum error {
  kDecayPipeRadius,
  kWaterLayer,
  kHornCurrent,
  kTargetDensity,
  kHorn1XShift,
  kHorn2XShift,
  kHorn1YShift,
  kHorn2YShift,
  kBeamTheta,
  kBeamThetaPhi,
  kBeamSigma,
  kBeamOffsetX,
  kPPFX,
  kPOTCounting,
  kFocussingTweaks,
  kAlignmentTweaks,
  kBeamAlignmentTweaks,
  kPPFX_abs,
  kPPFX_att,
  kPPFX_ttpCpi,
  kPPFX_ttpCk,
  kPPFX_ttnCpi,
  kPPFX_ttpCnu,
  kPPFX_ttnua,
  kPPFX_ttmesinc,
  kPPFX_oth,
  kTotal,
};
std::map<error, std::string> errstrings = {
    {kDecayPipeRadius, "DecayPipeRadius"},
    {kWaterLayer, "WaterLayer"},
    {kHornCurrent, "HornCurrent"},
    {kTargetDensity, "TargetDensity"},
    {kHorn1XShift, "Horn1XShift"},
    {kHorn2XShift, "Horn2XShift"},
    {kHorn1YShift, "Horn1YShift"},
    {kHorn2YShift, "Horn2YShift"},
    {kBeamTheta, "BeamTheta"},
    {kBeamThetaPhi, "BeamThetaPhi"},
    {kBeamSigma, "BeamSigma"},
    {kBeamOffsetX, "BeamOffsetX"},
    {kPPFX, "PPFX"},
    {kPOTCounting, "POTCounting"},
    {kFocussingTweaks, "Focussing"},
    {kAlignmentTweaks, "Alignment"},
    {kBeamAlignmentTweaks, "BeamAlignment"},
    {kPPFX_abs, "PPFX_abs"},
    {kPPFX_att, "PPFX_att"},
    {kPPFX_ttpCpi, "PPFX_ttpCpi"},
    {kPPFX_ttpCk, "PPFX_ttpCk"},
    {kPPFX_ttnCpi, "PPFX_ttnCpi"},
    {kPPFX_ttpCnu, "PPFX_ttpCnu"},
    {kPPFX_ttnua, "PPFX_ttnua"},
    {kPPFX_ttmesinc, "PPFX_ttmesinc"},
    {kPPFX_oth, "PPFX_oth"},
    {kTotal, "Total"},
};

std::map<error, std::string> errslegtrings = {
    {kDecayPipeRadius, "Decay pipe radius (10 cm)"},
    {kWaterLayer, "Water Layer Thickness (0.5 mm)"},
    {kHornCurrent, "Horn Current (3 kA)"},
    {kTargetDensity, "\\textrm{Target Density (0.4 g cm}^{-3})"},
    {kHorn1XShift, "Horn1 X Shift (1 mm)"},
    {kHorn2XShift, "Horn2 X Shift (1 mm)"},
    {kHorn1YShift, "Horn1 Y Shift (1 mm)"},
    {kHorn2YShift, "Horn2 Y Shift (1 mm)"},
    {kBeamTheta, "Beam Tilt (0.07 mrad)"},
    {kBeamThetaPhi,
     "\\textrm{Beam Tilt + Rotate }(0.07 \\textrm{mrad} + 90^{#circ})"},
    {kBeamSigma, "Beam Width (0.1 mm)"},
    {kBeamOffsetX, "Beam Offset X (0.45 mm)"},
    {kPPFX, "Hadron production (total)"},
    {kPOTCounting, "POT counting"},
    {kFocussingTweaks, "Focussing (total)"},
    {kAlignmentTweaks, "Alignment (total)"},
    {kBeamAlignmentTweaks, "Beam Alignment (total)"},
    {kPPFX_abs, "Other absorption"},
    {kPPFX_att, "Target absorption"},
    {kPPFX_ttpCpi, "\\textrm{pC}\\rightarrow\\pi"},
    {kPPFX_ttpCk, "\\textrm{pC}\\rightarrow \\textrm{K}"},
    {kPPFX_ttnCpi, "PPFX_ttnCpi"},
    {kPPFX_ttpCnu, "PPFX_ttpCnu"},
    {kPPFX_ttnua, "NucleonA"},
    {kPPFX_ttmesinc, "Meson inclusive"},
    {kPPFX_oth, "Other"},
    {kTotal, "Total"},
};

std::map<std::string, std::vector<error>> errfiles = {
    {"FluxRatios.root",
     {kDecayPipeRadius, kWaterLayer, kHornCurrent, kTargetDensity, kHorn1XShift,
      kHorn2XShift, kHorn1YShift, kHorn2YShift, kBeamTheta, kBeamThetaPhi,
      kBeamSigma, kBeamOffsetX, kPPFX, kPOTCounting, kFocussingTweaks,
      kAlignmentTweaks, kBeamAlignmentTweaks, kTotal}},
    {"FluxRatios_PPFX_nu.root",
     {kPPFX_abs, kPPFX_att, kPPFX_ttpCpi, kPPFX_ttpCk, kPPFX_ttnCpi,
      kPPFX_ttpCnu, kPPFX_ttnua, kPPFX_ttmesinc, kPPFX_oth}}};

enum location { kND, kNDPos6m, kNDPos12m, kNDPos18m, kNDPos24m, kNDPos30m };
std::map<location, std::string> locstrings = {
    {kND, "ND_relative"},
    {kNDPos6m, "ND_OA6.000000_relative"},
    {kNDPos12m, "ND_OA12.000000_relative"},
    {kNDPos18m, "ND_OA18.000000_relative"},
    {kNDPos24m, "ND_OA24.000000_relative"},
    {kNDPos30m, "ND_OA30.000000_relative"},
};
std::map<location, std::string> loclatexstr = {
    {kND, "\\textrm{ND, on axis, }"},
    {kNDPos6m, "\\textrm{ND, 6 m off axis, }"},
    {kNDPos12m, "\\textrm{ND, 12 m off axis, }"},
    {kNDPos18m, "\\textrm{ND, 18 m off axis, }"},
    {kNDPos24m, "\\textrm{ND, 24 m off axis, }"},
    {kNDPos30m, "\\textrm{ND, 30 m off axis, }"},
};

enum beam { kNu, kNuBar };
std::map<beam, std::string> beamstr = {
    {kNu, "nu"},
    {kNuBar, "nubar"},
};
std::map<beam, std::string> beamlatexstr = {
    {kNu, "\\nu\\textrm{-mode}"}, {kNuBar, "\\bar{\\nu}\\textrm{-mode}"}};

void PlotNearOffAxisErrorSources() {

  DeclareColors();

  gStyle->SetOptStat(false);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);

  // std::map<parent, int> Colors;
  // Colors[kTotal] = kBlack;
  // Colors[kPi] = kDUNEOrange;
  // Colors[kCK] = kDUNEBlue;
  // Colors[kK] = kMSUGreen;
  // Colors[kMu] = kMSUPurple;

  std::map<beam, std::map<spec, std::map<error, std::map<location, TH1 *>>>>
      Hists;

  bool doppfxall = false;

  for (auto errf : errfiles) {
    if (!doppfxall && (errf.first == "FluxRatios_PPFX_nu.root")) {
      continue;
    }
    TFile f(errf.first.c_str());
    for (auto e : errf.second) {
      for (auto b : beamstr) {
        for (auto s : specnames) {
          for (auto l : locstrings) {
            std::stringstream hname("");
            hname << b.second << "/" << s.second << "/" << errstrings[e]
                  << "/LBNF_" << s.second << "_flux_" << l.second;

            Hists[b.first][s.first][e][l.first] =
                dynamic_cast<TH1 *>(f.Get(hname.str().c_str()));

            if (!Hists[b.first][s.first][e][l.first]) {
              continue;
            }
            std::cout << "[INFO]: Read " << hname.str() << " from "
                      << errf.first << std::endl;

            Hists[b.first][s.first][e][l.first]->SetDirectory(nullptr);
          }
        }
      }
    }
  }

  std::map<std::pair<std::string, std::string>,
           std::vector<std::pair<error, int>>>
      PlotGroups = {
          {{"ErrType", ""},
           {{kTotal, kBlack},
            {kPPFX, kRed},
            {kPOTCounting, kGray + 2},
            {kFocussingTweaks, kDUNEBlue},
            {kAlignmentTweaks, kDUNEOrange},
            {kBeamAlignmentTweaks, kMSUGreen}}},
      };

  TCanvas c1("", "", 600, 700);
  // c1.SaveAs("dum.pdf[");
  for (auto pg : PlotGroups) {
    for (auto b : beamstr) {
      for (auto s : specnames) {

        c1.Clear();
        c1.cd();

        TLatex ltx;
        ltx.SetTextSize(0.09);
        ltx.SetTextAlign(22);

        TLegend leg(0, 0.825, 1, 1);
        leg.SetHeader(pg.first.second.c_str());
        leg.SetBorderSize(0);
        leg.SetNColumns(2);
        leg.SetTextSize(0.04);

        size_t it = 0;
        for (auto l : locstrings) {
          c1.cd();
          bool isLeft = !(it & 1);

          double h = 0.275;
          double bot = (2 - (it / 2)) * h;
          double top = (3 - (it / 2)) * h;

          TVirtualPad *pnd =
              new TPad("nd", "", isLeft ? 0 : 0.5, bot, isLeft ? 0.5 : 1, top);
          pnd->SetLeftMargin(0.2);
          pnd->SetRightMargin(0.06);
          pnd->SetTopMargin(0.1);
          pnd->SetBottomMargin(0.2);
          pnd->AppendPad();
          pnd->cd();

          bool drawn = false;
          for (auto e : pg.second) {
            std::cout << "Plotting " << s.second << " for " << b.second
                      << " @ error: " << errstrings[e.first] << std::endl;
            if (!Hists[b.first][s.first][e.first][l.first]) {
              std::cout << "Couldn't get!" << std::endl;
              abort();
            }
            Hists[b.first][s.first][e.first][l.first]->SetLineWidth(2);
            Hists[b.first][s.first][e.first][l.first]->SetLineColor(e.second);

            if (!drawn) {
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetRangeUser(0, 0.15);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetNdivisions(505);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetTitleSize(0.09);
              Hists[b.first][s.first][e.first][l.first]->GetYaxis()->SetTitle(
                  "Fractional uncertainty");
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetTitleOffset(1.1);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetLabelOffset(0);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetYaxis()
                  ->SetLabelSize(0.09);

              Hists[b.first][s.first][e.first][l.first]->GetXaxis()->SetTitle(
                  "E_{\\nu}\\textrm{ (GeV)}");
              Hists[b.first][s.first][e.first][l.first]
                  ->GetXaxis()
                  ->SetNdivisions(505);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetXaxis()
                  ->SetTitleSize(0.09);
              Hists[b.first][s.first][e.first][l.first]
                  ->GetXaxis()
                  ->SetLabelSize(0.09);
            }
            Hists[b.first][s.first][e.first][l.first]->Draw(drawn ? "HISTSAME"
                                                                  : "HIST");
            if (l.first == kND) {
              leg.AddEntry(Hists[b.first][s.first][e.first][l.first],
                           errslegtrings[e.first].c_str(), "l");
            }
            drawn = true;
          }
          std::stringstream phdr("");
          phdr.str("");
          phdr << loclatexstr[l.first] << beamlatexstr[b.first] << ", "
               << speclatexnames[s.first];
          ltx.DrawLatexNDC(0.6, 0.95, phdr.str().c_str());
          it++;
        }

        c1.cd();
        leg.Draw();
        std::stringstream oname("");
        oname << "tikz/" << b.second << "mode_" << s.second << "_"
              << pg.first.first << "_OffAxis.tex";

        c1.SaveAs(oname.str().c_str());
        // c1.SaveAs("dum.pdf");
      }
    }
  }
  // c1.SaveAs("dum.pdf]");
}
