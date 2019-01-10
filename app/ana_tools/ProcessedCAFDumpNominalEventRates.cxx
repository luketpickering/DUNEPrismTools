#include "CAFReader.hxx"
#include "PlotProcessor.hxx"
#include "ROOTUtility.hxx"

#include "DUNETDRNDHelper.hxx"
#include "OscillationHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

struct FDOscHelper {
  OscillationHelper oh;
  FDOscHelper(fhicl::ParameterSet const &ps) { oh.Setup(ps); }
  std::array<double, 2> operator()(CAFReader const &ev) {
    oh.SetOscillationChannel(oh.GetNuType(ev.nuPDGunosc), oh.GetNuType(14));
    return {ev.Ev_reco_numu, ev.POTWeight * oh.GetWeight(ev.Ev)};
  }
};

// Expect a fhicl like
// Input_ND: ND_CAF.root
// Input_FD: FD_CAF.root
// Oscillation:
int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string inp_ND = ps.get<std::string>("Input_ND");
  CAFReader ND_rdr(inp_ND);

  std::string inp_FD = ps.get<std::string>("Input_FD");
  CAFReader FD_rdr(inp_FD);

  FD_rdr.GetEntry(0);

  bool isFHC = FD_rdr.isFHC;

  size_t FDEvents = FD_rdr.GetEntries();
  Plot1DCAF FDEventRate(FDOscHelper(ps.get<fhicl::ParameterSet>("Oscillation")),
                        "OffAxisEvRate", "", 4500, -5, 40);

  for (size_t fd_it = 0; fd_it < FDEvents; ++fd_it) {
    FD_rdr.GetEntry(fd_it);
    FDEventRate.Process(FD_rdr);
  }

  // Plot2DCAF NDEventRate(
  //     [](CAFReader const &ev) -> std::array<double, 3> {
  //       return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
  //     },
  //     "OffAxisEvRateEnu", "", 20, 0, 6, 90, -5, 40);

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  FDEventRate.Write(oupFile);

  oupFile->Write();
  oupFile->Close();
}
