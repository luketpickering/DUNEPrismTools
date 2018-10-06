#include "FluxFitResultsTreeReader.hxx"
#include "OscillationParametersTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"

#include "GetUsage.hxx"

#include "FluxCombiner.hxx"
#include "FluxFitter.hxx"

#ifdef USE_FHICL
#include "fhiclcpp/make_ParameterSet.h"

std::string fhicl_file = "";
#endif

bool IsGauss = false;

std::string OutputFile, OutputDirectory = "";
bool UPDATEOutputFile = false;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << GetUsageText(argv[0], "flux_tools")
            << std::endl;
#ifdef USE_FHICL
  std::cout
      << "  --fhicl <config.fcl>          : Configure fits via a FHiCL \n"
         "                                  configuration file. N.B. All \n"
         "                                  other command line options will \n"
         "                                  be ignored."
      << std::endl;
#endif
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-g") {
      opt++; // Skip passed argument
      IsGauss = true;
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
      UPDATEOutputFile = false;
    } else if (std::string(argv[opt]) == "-a") {
      OutputFile = argv[++opt];
      UPDATEOutputFile = true;
    } else if (std::string(argv[opt]) == "-d") {
      OutputDirectory = argv[++opt];
    }
#ifdef USE_FHICL
    else if (std::string(argv[opt]) == "--fhicl") {
      fhicl_file = argv[++opt];
    }
#endif
    else if ((std::string(argv[opt]) == "-?") ||
             std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    }
    opt++;
  }
}
int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);
#ifdef USE_FHICL
  if (fhicl_file.size()) {

    TTree *fit_coeff_tree = nullptr;
    int NCoeffs = 0;
    std::vector<double> Coeffs(1000);

    fhicl::ParameterSet config_ps = fhicl::make_ParameterSet(fhicl_file);

    fhicl::ParameterSet default_flux_inputs =
        config_ps.get<fhicl::ParameterSet>("FluxInputs", fhicl::ParameterSet());
    std::string default_OutputFile =
        config_ps.get<std::string>("OutputFile", "");
    bool RecreateFile = config_ps.get<bool>("RecreateOutputFile", false);

    std::vector<fhicl::ParameterSet> const &fits =
        config_ps.get<std::vector<fhicl::ParameterSet>>("Fits");
    for (size_t fit_it = 0; fit_it < fits.size(); ++fit_it) {
      fhicl::ParameterSet const &fit_ps = fits[fit_it];

      if (!fit_coeff_tree) {
        fit_coeff_tree = new TTree("fit_coefficients", "");
        fit_coeff_tree->Branch("NCoeffs", &NCoeffs, "NCoeffs/I");
        fit_coeff_tree->Branch("Coeffs", Coeffs.data(), "NCoeffs[NCoeffs]/D");
        fit_coeff_tree->SetDirectory(nullptr);
      }

      OutputFile = fit_ps.get<std::string>("OutputFile", default_OutputFile);

      UPDATEOutputFile = fit_ps.get<bool>("UPDATEOutputFile", !RecreateFile);

      if (OutputFile ==
          default_OutputFile) { // If using the default output directory make
                                // sure it is only recreated on the first fit.
        RecreateFile = false;
      }

      OutputDirectory = fit_ps.get<std::string>("OutputDirectory", "");

      fhicl::ParameterSet flux_inputs =
          fit_ps.get<fhicl::ParameterSet>("FluxInputs", default_flux_inputs);

      FluxFitter ff;

      std::cout << "[INFO]: Initializing fitter with: "
                << fit_ps.to_indented_string() << std::endl;

      std::vector<std::unique_ptr<TH1D>> target_flux_components;

      if (fit_ps.get<bool>("IsGauss", false)) {
        ff.InitializeGauss(MakeFluxFitterOptions(fit_ps),
                           MakeGausTargetOptions(fit_ps, flux_inputs));
      } else {
        bool FitterExpectsFluxLater =
            fit_ps.get<bool>("BuildTargetFlux", false);

        ff.InitializeFlux(MakeFluxFitterOptions(fit_ps),
                          MakeFluxTargetOptions(fit_ps, flux_inputs),
                          FitterExpectsFluxLater);

        if (FitterExpectsFluxLater) {
          ff.SetTargetFlux(GetCombinedFlux(fit_ps.get<fhicl::ParameterSet>(
                                               "FluxDescription"),
                                           target_flux_components)
                               .get(),
                           true);
        }
      }

      ff.Fit();

      NCoeffs = ff.GetCoefficients().size();
      std::copy_n(ff.GetCoefficients().begin(), NCoeffs, Coeffs.begin());
      fit_coeff_tree->Fill();

      if (OutputFile.size()) {
        std::unique_ptr<TFile> oupF = std::unique_ptr<TFile>(CheckOpenFile(
            OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE"));
        TDirectory *oupD = oupF.get();

        if (OutputDirectory.length()) {
          oupD = oupF->mkdir(OutputDirectory.c_str());
        }

        ff.Write(oupD);

        if (target_flux_components.size() > 1) {
          TDirectory *compD = oupD->mkdir("InputTargetFluxComponents");
          for (size_t comp_i = 0; comp_i < target_flux_components.size();
               ++comp_i) {
            std::stringstream ss("");
            ss << target_flux_components[comp_i]->GetName() << "_comp_"
               << comp_i;
            TH1D *sf = static_cast<TH1D *>(
                target_flux_components[comp_i]->Clone(ss.str().c_str()));
            sf->SetDirectory(compD);
          }
        }

        if ((fit_it + 1) == fits.size()) {
          fit_coeff_tree->SetDirectory(oupF.get());
        }

        oupF->Write();
      }
    }

#else
  if (false) { // If we don't have FHiCL
#endif
  } else {

    if (!OutputFile.length()) {
      std::cout << "[ERROR]: No output file specified." << std::endl;
      throw;
    }

    FluxFitter ff;

    if (IsGauss) {
      ff.InitializeGauss(MakeFluxFitterOptions(argc, argv),
                         MakeGausTargetOptions(argc, argv));
    } else {
      ff.InitializeFlux(MakeFluxFitterOptions(argc, argv),
                        MakeFluxTargetOptions(argc, argv));
    }

    ff.Fit();

    TFile *oupF =
        CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE");
    TDirectory *oupD = oupF;

    if (OutputDirectory.length()) {
      oupD = oupF->mkdir(OutputDirectory.c_str());
    }

    ff.Write(oupD);

    oupF->Write();
    oupF->Close();
  }
}
