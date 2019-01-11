#include "CAFReader.hxx"
#include "PlotProcessor.hxx"
#include "ROOTUtility.hxx"

#include "DUNETDRDetHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

size_t const kAll = 0;
size_t const kTrueSel = 1;
size_t const kNumuSel = 2;
size_t const kNSel = 3;

std::vector<std::string> StageSelectionString = {"_All", "_TrueSel",
                                                 "_NumuSel"};

// Expect a fhicl like
// InputFile: ND_FHC_CAF.root
int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string inpf = ps.get<std::string>("InputFile");
  CAFReader rdr(inpf);

  size_t NMaxEvents =
      ps.get<size_t>("NMaxEvents", std::numeric_limits<size_t>::max());
  size_t fents = std::min(NMaxEvents, rdr.GetEntries());

  std::vector<Plot1DCAF> Plot1DDefinitions;

  Plot1DDefinitions.emplace_back(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.det_x + ev.vtx_x * 1E-2, 1};
      },
      "OffAxisEvRate_noweight", "", 900, -5, 40);

  Plot1DDefinitions.emplace_back(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRate", "", 900, -5, 40);

  Plot1DDefinitions.emplace_back(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.vtx_x * 1E-2, (ev.det_x == 0)};
      },
      "vtx_x_noweight", "", 200, -3, 3);

  std::vector<std::vector<Plot1DCAF>> All1DPlots;
  All1DPlots.resize(kNSel);

  std::vector<Plot2DCAF> Plot2DDefinitions;

  Plot2DDefinitions.emplace_back(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateEnu", "", 20, 0, 6, 90, -5, 40);

  std::vector<std::vector<Plot2DCAF>> All2DPlots;
  All2DPlots.resize(kNSel);

  for (size_t i = kAll; i < kNSel; ++i) {
    for (Plot1DCAF const &p : Plot1DDefinitions) {
      All1DPlots[i].emplace_back(
          p.Clone(p.GetName() + StageSelectionString[i]));
    }
    for (Plot2DCAF const &p : Plot2DDefinitions) {
      All2DPlots[i].emplace_back(
          p.Clone(p.GetName() + StageSelectionString[i]));
    }
  }
  Plot1DDefinitions.clear();
  Plot2DDefinitions.clear();

  for (size_t ev_it = 0; ev_it < fents; ++ev_it) {
    rdr.GetEntry(ev_it);

    if (ev_it && (fents / 20) && !(ev_it % (fents / 20))) {
      std::cout << "\r[INFO]: Processing event " << ev_it << "/" << fents
                << " (" << std::round((double(ev_it) / double(fents)) * 100)
                << " %)." << std::flush;
    }

    if (!FV_Select(rdr)) {
      continue;
    }

    for (Plot1DCAF &p : All1DPlots[kAll]) {
      p.Process(rdr);
    }
    for (Plot2DCAF &p : All2DPlots[kAll]) {
      p.Process(rdr);
    }

    if (FHC_Numu_True_Select(rdr)) {
      for (Plot1DCAF &p : All1DPlots[kTrueSel]) {
        p.Process(rdr);
      }
      for (Plot2DCAF &p : All2DPlots[kTrueSel]) {
        p.Process(rdr);
      }
    }
    if (ND_FHC_Numu_Select(rdr)) {
      for (Plot1DCAF &p : All1DPlots[kNumuSel]) {
        p.Process(rdr);
      }
      for (Plot2DCAF &p : All2DPlots[kNumuSel]) {
        p.Process(rdr);
      }
    }
  }
  std::cout << "\n[INFO]: Processed " << fents << " events." << std::endl;

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  for (auto &SelectionPlots : All1DPlots) {
    for (Plot1DCAF &p : SelectionPlots) {
      p.Write(oupFile);
    }
  }
  for (auto &SelectionPlots : All2DPlots) {
    for (Plot2DCAF &p : SelectionPlots) {
      p.Write(oupFile);
    }
  }

  oupFile->Write();
  oupFile->Close();
}
