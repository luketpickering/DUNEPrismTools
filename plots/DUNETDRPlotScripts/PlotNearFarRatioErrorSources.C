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
  kTotalNonHP
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
    {kTotalNonHP, "Total"}};

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
    {kTotalNonHP, "Total"}};

std::map<std::string, std::vector<error>> errfiles = {
    {"FluxRatios.root",
     {kDecayPipeRadius, kWaterLayer, kHornCurrent, kTargetDensity, kHorn1XShift,
      kHorn2XShift, kHorn1YShift, kHorn2YShift, kBeamTheta, kBeamThetaPhi,
      kBeamSigma, kBeamOffsetX, kPPFX, kPOTCounting, kFocussingTweaks,
      kAlignmentTweaks, kBeamAlignmentTweaks, kTotal}},
    {"FluxRatios_NonHP.root", {kTotalNonHP}},
    {"FluxRatios_PPFX_nu.root",
     {kPPFX_abs, kPPFX_att, kPPFX_ttpCpi, kPPFX_ttpCk, kPPFX_ttnCpi,
      kPPFX_ttpCnu, kPPFX_ttnua, kPPFX_ttmesinc, kPPFX_oth}}};

enum location { kND, kFD, kNDFDRatio };
std::map<location, std::string> locstrings = {
    {kND, "ND_relative"},
    {kFD, "FD_relative"},
    {kNDFDRatio, "ND_FD_Ratio"},
};
enum beam { kNu, kNuBar };
std::map<beam, std::string> beamstr = {
    {kNu, "nu"},
    {kNuBar, "nubar"},
};
std::map<beam, std::string> beamlatexstr = {
    {kNu, "\\nu\\textrm{-mode}"}, {kNuBar, "\\bar{\\nu}\\textrm{-mode}"}};

void PlotNearFarRatioErrorSources() {

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

  bool doother = false;

  for (auto errf : errfiles) {
    if (!doother && (errf.first != "FluxRatios.root")) {
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

  std::map<std::tuple<std::string, std::string, double>,
           std::vector<std::pair<error, int>>>
      PlotGroups = {
          {{"ErrType", "", 0.15},
           {{kTotal, kBlack},
            {kPPFX, kRed},
            {kPOTCounting, kGray + 2},
            {kFocussingTweaks, kDUNEBlue},
            {kAlignmentTweaks, kDUNEOrange},
            {kBeamAlignmentTweaks, kMSUGreen}}},
          {{"Focussing", "Focussing errors", 0.08},
           {{kDecayPipeRadius, kRed},
            {kWaterLayer, kDUNEBlue},
            {kHornCurrent, kDUNEOrange},
            {kTargetDensity, kMSUGreen}}},
          {{"HornAlignment", "Horn Alignment errors", 0.08},
           {{kHorn1XShift, kRed},
            {kHorn2XShift, kDUNEBlue},
            {kHorn1YShift, kDUNEOrange},
            {kHorn2YShift, kMSUGreen}}},
          {{"BeamAlignment", "Beam Alignment errors", 0.08},
           {{kBeamTheta, kRed},
            {kBeamThetaPhi, kDUNEBlue},
            {kBeamSigma, kDUNEOrange},
            {kBeamOffsetX, kMSUGreen}}},
          {{"PPFX", "Hadron production errors", 0.15},
           {
               {kPPFX, kBlack},
               {kPPFX_oth, kDUNEBlue},
               {kPPFX_ttpCpi, kDUNEOrange},
               {kPPFX_ttpCk, kMSUGreen},
               {kPPFX_ttmesinc, kMSUPurple},
               {kPPFX_ttnua, kRed},
               {kPPFX_att, kGray + 2},
               {kPPFX_abs, kLilac},
           }},
          {{"NonHPErrs", "", 0.08},
           {
               {kTotalNonHP, kBlack},
               {kHornCurrent, kRed},
               {kWaterLayer, kDUNEBlue},
               {kDecayPipeRadius, kMSUPurple},
               {kAlignmentTweaks, kDUNEOrange},
               {kBeamSigma, kAzure - 4},
               {kBeamOffsetX, kLilac},
               {kTargetDensity, kMSUGreen},
               {kPOTCounting, kGray + 2},
           }},
      };

  TCanvas c1("", "", 600, 600);
  for (auto pg : PlotGroups) {
    if (!doother && ((std::get<0>(pg.first) == "PPFX") ||
                     (std::get<0>(pg.first) == "NonHPErrs"))) {
      continue;
    }
    for (auto b : beamstr) {
      for (auto s : specnames) {

        if (std::get<0>(pg.first) == "PPFX") {
          if (b.first != kNu) {
            continue;
          }
          if (s.first != kNuMu) {
            continue;
          }
        }

        c1.Clear();
        c1.cd();

        TLatex ltx;
        ltx.SetTextSize(0.09);
        ltx.SetTextAlign(22);

        TLegend l(0.5, 0.05, 0.95, 0.45);
        l.SetHeader(std::get<1>(pg.first).c_str());
        l.SetBorderSize(0);
        l.SetNColumns(1);
        l.SetTextSize(0.045);

        TVirtualPad *pnd = new TPad("nd", "", 0, 0.5, 0.5, 1);
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
          if (!Hists[b.first][s.first][e.first][kND]) {
            std::cout << "couldn't!" << std::endl;
            abort();
          }
          Hists[b.first][s.first][e.first][kND]->SetLineWidth(2);
          Hists[b.first][s.first][e.first][kND]->SetLineColor(e.second);

          if (!drawn) {
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetRangeUser(
                0, std::get<2>(pg.first));
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetNdivisions(
                505);
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetTitleSize(
                0.09);
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetTitle(
                "Fractional uncertainty");
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetTitleOffset(
                1.3);
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetLabelOffset(
                0);
            Hists[b.first][s.first][e.first][kND]->GetYaxis()->SetLabelSize(
                0.09);

            Hists[b.first][s.first][e.first][kND]->GetXaxis()->SetTitle(
                "E_{\\nu}\\textrm{ (GeV)}");
            Hists[b.first][s.first][e.first][kND]->GetXaxis()->SetNdivisions(
                505);
            Hists[b.first][s.first][e.first][kND]->GetXaxis()->SetTitleSize(
                0.09);
            Hists[b.first][s.first][e.first][kND]->GetXaxis()->SetLabelSize(
                0.09);
          }
          Hists[b.first][s.first][e.first][kND]->Draw(drawn ? "HISTSAME"
                                                            : "HIST");
          l.AddEntry(Hists[b.first][s.first][e.first][kND],
                     errslegtrings[e.first].c_str(), "l");
          drawn = true;
        }
        std::stringstream phdr("");
        phdr.str("");
        phdr << "\\textrm{Near detector, }" << beamlatexstr[b.first] << ", "
             << speclatexnames[s.first];
        ltx.DrawLatexNDC(0.6, 0.95, phdr.str().c_str());

        c1.cd();
        TVirtualPad *pfd = new TPad("fd", "", 0.5, 0.5, 1, 1);
        pfd->SetLeftMargin(0.2);
        pfd->SetRightMargin(0.05);
        pfd->SetTopMargin(0.1);
        pfd->SetBottomMargin(0.2);
        pfd->AppendPad();
        pfd->cd();

        drawn = false;
        for (auto e : pg.second) {
          std::cout << "Plotting " << s.second << " for " << b.second
                    << " @ error: " << errstrings[e.first] << std::endl;
          Hists[b.first][s.first][e.first][kFD]->SetLineWidth(2);
          Hists[b.first][s.first][e.first][kFD]->SetLineColor(e.second);
          if (!drawn) {

            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetRangeUser(
                0, std::get<2>(pg.first));
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetNdivisions(
                505);
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetTitleSize(
                0.09);
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetTitle(
                "Fractional uncertainty");
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetTitleOffset(
                1.3);
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetLabelOffset(
                0);
            Hists[b.first][s.first][e.first][kFD]->GetYaxis()->SetLabelSize(
                0.09);

            Hists[b.first][s.first][e.first][kFD]->GetXaxis()->SetTitle(
                "E_{\\nu}\\textrm{ (GeV)}");
            Hists[b.first][s.first][e.first][kFD]->GetXaxis()->SetNdivisions(
                505);
            Hists[b.first][s.first][e.first][kFD]->GetXaxis()->SetTitleSize(
                0.09);
            Hists[b.first][s.first][e.first][kFD]->GetXaxis()->SetLabelSize(
                0.09);
          }
          Hists[b.first][s.first][e.first][kFD]->Draw(drawn ? "HISTSAME"
                                                            : "HIST");
          drawn = true;
        }
        phdr.str("");
        phdr << "\\textrm{Far detector, }" << beamlatexstr[b.first] << ", "
             << speclatexnames[s.first];
        ltx.DrawLatexNDC(0.6, 0.95, phdr.str().c_str());

        c1.cd();
        TVirtualPad *pndfd = new TPad("ndfd", "", 0, 0, 0.5, 0.5);
        pndfd->SetLeftMargin(0.2);
        pndfd->SetRightMargin(0.05);
        pndfd->SetTopMargin(0.1);
        pndfd->SetBottomMargin(0.2);
        pndfd->AppendPad();
        pndfd->cd();

        drawn = false;
        for (auto e : pg.second) {
          std::cout << "Plotting " << s.second << " for " << b.second
                    << " @ error: " << errstrings[e.first] << std::endl;
          Hists[b.first][s.first][e.first][kNDFDRatio]->SetLineWidth(2);
          Hists[b.first][s.first][e.first][kNDFDRatio]->SetLineColor(e.second);
          if (!drawn) {

            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetRangeUser(0, std::get<2>(pg.first));
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetNdivisions(505);
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetTitleSize(0.09);
            Hists[b.first][s.first][e.first][kNDFDRatio]->GetYaxis()->SetTitle(
                "Fractional uncertainty");
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetTitleOffset(1.3);
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetLabelOffset(0);
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetYaxis()
                ->SetLabelSize(0.09);

            Hists[b.first][s.first][e.first][kNDFDRatio]->GetXaxis()->SetTitle(
                "E_{\\nu}\\textrm{ (GeV)}");
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetXaxis()
                ->SetNdivisions(505);
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetXaxis()
                ->SetTitleSize(0.09);
            Hists[b.first][s.first][e.first][kNDFDRatio]
                ->GetXaxis()
                ->SetLabelSize(0.09);
          }
          Hists[b.first][s.first][e.first][kNDFDRatio]->Draw(drawn ? "HISTSAME"
                                                                   : "HIST");
          drawn = true;
        }
        phdr.str("");
        phdr << "\\textrm{Far/Near ratio, }" << beamlatexstr[b.first] << ", "
             << speclatexnames[s.first];
        ltx.DrawLatexNDC(0.6, 0.95, phdr.str().c_str());

        c1.cd();
        l.Draw();
        std::stringstream oname("");
        oname << "tikz/" << b.second << "mode_" << s.second << "_"
              << std::get<0>(pg.first) << ".tex";

        c1.SaveAs(oname.str().c_str());
      }
    }
  }
}
