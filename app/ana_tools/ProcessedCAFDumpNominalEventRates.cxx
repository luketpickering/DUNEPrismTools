#include "CAFReader.hxx"
#include "PlotProcessor.hxx"
#include "ROOTUtility.hxx"

#include "DUNETDRDetHelper.hxx"
#include "OscillationHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

// Expect a fhicl like
// Input_ND: ND_CAF.root
// Input_FD: FD_CAF.root
// Oscillation:
int main(int argc, char const *argv[]) {

  TH1::SetDefaultSumw2();

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  std::cout << FD_ND_FVRatio(50) << std::endl;

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string inp_ND = ps.get<std::string>("Input_ND");
  CAFReader ND_rdr(inp_ND);

  std::string inp_FD = ps.get<std::string>("Input_FD");
  CAFReader FD_rdr(inp_FD);

  FD_rdr.GetEntry(0);

  bool isFHC = FD_rdr.isFHC;

  std::vector<std::unique_ptr<CAFEvProcessor>> FDPlots;

  Plot1DCAF FD_ETrue(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev, ev.POTWeight};
      },
      "Enu_true", ";E_{#nu} (GeV); Count/POT", 60, 0, 6);

  size_t kFDTrueAll = 0;
  FDPlots.emplace_back(new Plot1DCAF(FD_ETrue.Clone("Enu_true")));

  size_t kFDTrueSel = 1;
  FDPlots.emplace_back(new Plot1DCAF(FD_ETrue.Clone("Enu_sel")));
  size_t kFDTrueSelBad = 2;
  FDPlots.emplace_back(new Plot1DCAF(FD_ETrue.Clone("Enu_sel_bad")));

  Plot1DCAF FD_ERec(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev_reco_numu, ev.POTWeight};
      },
      "Erec_true", ";E_{#nu}^{Reco} (GeV); Count/POT", 60, 0, 6);

  size_t kFDRecSel = 3;
  FDPlots.emplace_back(new Plot1DCAF(FD_ERec.Clone("Erec_sel")));
  size_t kFDRecSelBad = 4;
  FDPlots.emplace_back(new Plot1DCAF(FD_ERec.Clone("Erec_sel_bad")));

  size_t kFDTrueRecSel = 5;
  FDPlots.emplace_back(new Plot2DCAF(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev_reco_numu, ev.Ev, ev.POTWeight};
      },
      "Enu_true_rec", ";E_{#nu}^{Reco} (GeV);E_{#nu} (GeV); Count/POT", 200, 0,
      10, 200, 0, 10));

  size_t kFDTruefluxbin = 6;
  FDPlots.emplace_back(new Plot1DCAF(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev, ev.POTWeight};
      },
      "Enu_true_fluxbin", ";E_{#nu} (GeV); Count/POT", 400, 0, 10));

  OscillationHelper oh;
  oh.Setup(ps.get<fhicl::ParameterSet>("Oscillation"));

  size_t FDEvents = FD_rdr.GetEntries();
  for (size_t fd_it = 0; fd_it < FDEvents; ++fd_it) {
    FD_rdr.GetEntry(fd_it);
    if (!FV_Select(FD_rdr)) {
      continue;
    }

    oh.SetOscillationChannel(FD_rdr.nuPDGunosc, FD_rdr.nuPDG);
    double ow = 1;
    oh.GetWeight(FD_rdr.Ev);

    bool true_sel = FHC_Numu_True_Select(FD_rdr);
    if (true_sel) {
      FDPlots[kFDTrueAll]->Process(FD_rdr, ow);
      FDPlots[kFDTruefluxbin]->Process(FD_rdr, ow);
    }

    if (FD_FHC_Numu_Select(FD_rdr)) {
      FDPlots[kFDTrueSel]->Process(FD_rdr, ow);
      FDPlots[kFDRecSel]->Process(FD_rdr, ow);
      if (!true_sel) {
        FDPlots[kFDTrueSelBad]->Process(FD_rdr, ow);
        FDPlots[kFDRecSelBad]->Process(FD_rdr, ow);
      }
      FDPlots[kFDTrueRecSel]->Process(FD_rdr, ow);
    }
  }

  std::vector<std::unique_ptr<CAFEvProcessor>> NDPlots;

  Plot2DCAF NDETrue(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateEnu", "", 60, 0, 6, 90, -5.25, 39.75);
  Plot2DCAF NDERec(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev_reco, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateErec", "", 60, 0, 6, 90, -5.25, 39.75);

  size_t kNDTrueAll = 0;
  NDPlots.emplace_back(new Plot2DCAF(NDETrue.Clone("OffAxisEvRate")));
  size_t kNDTrueSel = 1;
  NDPlots.emplace_back(new Plot2DCAF(NDETrue.Clone("OffAxisEvRateEnu_sel")));
  size_t kNDTrueSelBad = 2;
  NDPlots.emplace_back(
      new Plot2DCAF(NDETrue.Clone("OffAxisEvRateEnu_sel_bad")));

  size_t kNDRecSel = 3;
  NDPlots.emplace_back(new Plot2DCAF(NDERec.Clone("OffAxisEvRateERec_sel")));
  size_t kNDRecSelBad = 4;
  NDPlots.emplace_back(
      new Plot2DCAF(NDERec.Clone("OffAxisEvRateERec_sel_bad")));

  Plot1DCAF NDETrue1D(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev, ev.POTWeight * (ev.det_x == 0) * (ev.vtx_x < 25)};
      },
      "OffAxisEvRateErec", "", 60, 0, 6);

  size_t kND1DTrueAll = 5;
  NDPlots.emplace_back(new Plot1DCAF(NDETrue1D.Clone("NDETrue1_all")));
  size_t kND1DTrueSel = 6;
  NDPlots.emplace_back(new Plot1DCAF(NDETrue1D.Clone("NDETrue1_sel")));

  size_t kND1DRecSel = 7;
  NDPlots.emplace_back(new Plot1DCAF(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev_reco, ev.POTWeight * (ev.det_x == 0) * (ev.vtx_x < 25)};
      },
      "NDERec1_sel", "", 60, 0, 6));

  size_t kNDTrueAll_fluxbin = 8;
  NDPlots.emplace_back(new Plot2DCAF(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateEnu_fluxbin", "", 400, 0, 10, 91, -0.25, 45.25));

  size_t NDEvents = ND_rdr.GetEntries();
  for (size_t nd_it = 0; nd_it < NDEvents; ++nd_it) {
    ND_rdr.GetEntry(nd_it);
    if (!FV_Select(ND_rdr)) {
      continue;
    }

    bool true_sel = FHC_Numu_True_Select(ND_rdr);
    if (true_sel) {
      NDPlots[kNDTrueAll]->Process(ND_rdr);
      NDPlots[kND1DTrueAll]->Process(ND_rdr);
      NDPlots[kNDTrueAll_fluxbin]->Process(ND_rdr);
    }

    if (ND_FHC_Numu_Select(ND_rdr)) {
      NDPlots[kNDTrueSel]->Process(ND_rdr);
      NDPlots[kNDRecSel]->Process(ND_rdr);
      NDPlots[kND1DTrueSel]->Process(ND_rdr);
      NDPlots[kND1DRecSel]->Process(ND_rdr);
      if (!true_sel) {
        NDPlots[kNDTrueSelBad]->Process(ND_rdr);
        NDPlots[kNDRecSelBad]->Process(ND_rdr);
      }
    }
  }

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  TDirectory *fdplotsdir = oupFile->mkdir("fdplots");
  for (auto &p : FDPlots) {
    p->Write(fdplotsdir);
  }

  TDirectory *ndplotsdir = oupFile->mkdir("ndplots");
  for (auto &p : NDPlots) {
    p->Write(ndplotsdir);
  }
  oupFile->Write();
  oupFile->Close();
}
