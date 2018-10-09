#include "FluxCombiner.hxx"

#include "OscillationHelper.hxx"

#include "ROOTUtility.hxx"

#include <string>

std::unique_ptr<TH1D>
GetCombinedFlux(fhicl::ParameterSet const &flux_descriptor,
                std::vector<std::unique_ptr<TH1D>> &components) {

  components.clear();

  std::unique_ptr<TH1D> SummedFlux(nullptr);

  std::string default_input_file =
      flux_descriptor.get<std::string>("InputFluxFile");
  std::string name = flux_descriptor.get<std::string>("Name", "");

  OscillationHelper oh;

  bool doosc = false;
  if (flux_descriptor.has_key("Oscillation")) {
    oh.Setup(flux_descriptor.get<fhicl::ParameterSet>("Oscillation"));
    doosc = true;
  }

  for (fhicl::ParameterSet const &fs :
       flux_descriptor.get<std::vector<fhicl::ParameterSet>>("Fluxes")) {

    std::string input_file =
        fs.get<std::string>("InputFluxFile", default_input_file);
    std::string input_hist = fs.get<std::string>("InputHistogram");

    std::vector<std::unique_ptr<TH1D>> Fluxes;

    std::pair<int, int> OscChannel =
        fs.get<std::pair<int, int>>("Oscillate", {0, 0});

    // Assume 2D
    if (fs.has_key("CombineFluxSlicesDescriptor")) {
      std::vector<std::pair<double, double>> XRanges =
          BuildRangesList(fs.get<std::string>("CombineFluxSlicesDescriptor"));
      std::unique_ptr<TH2D> flux2D =
          GetHistogram_uptr<TH2D>(input_file, input_hist);
      Fluxes = MergeSplitTH2D(flux2D, true, XRanges);
    } else {
      Fluxes.emplace_back(GetHistogram_uptr<TH1D>(input_file, input_hist));
    }

    for (std::unique_ptr<TH1D> &h : Fluxes) {
      if (OscChannel.first) {
        if (!doosc) {
          std::cout
              << "[ERROR]: When building flux, an oscillation was "
                 "requested, but no \"Oscillation\" parameter set was found "
                 "in the flux descriptor element: "
              << flux_descriptor.to_indented_string() << std::endl;
          throw;
        }
        oh.SetOscillationChannel(OscChannel.first, OscChannel.second);
        std::stringstream ss("");
        ss << h->GetName() << "_unosc";
        components.emplace_back(
            static_cast<TH1D *>(h->Clone(ss.str().c_str())));
        components.back()->SetDirectory(nullptr);
        oh.OscillateHistogram(h);
      }
      if (!SummedFlux) {
        if (name.size()) {
          SummedFlux = std::unique_ptr<TH1D>(
              static_cast<TH1D *>(h->Clone(name.c_str())));
        } else {
          SummedFlux = std::unique_ptr<TH1D>(static_cast<TH1D *>(h->Clone()));
        }
        SummedFlux->SetDirectory(nullptr);
      } else {
        SummedFlux->Add(h.get());
      }
      components.push_back(std::move(h));
    }
  }
  return SummedFlux;
}
