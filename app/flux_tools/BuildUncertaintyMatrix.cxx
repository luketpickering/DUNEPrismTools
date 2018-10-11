#include "GetUsage.hxx"

#include "CovarianceHelper.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TH1D.h"
#include "TMatrixD.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "Eigen/Core"

#include <memory>
#include <string>

std::string fhicl_config = "";

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "flux_tools") << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "--fhicl") {
      fhicl_config = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!fhicl_config.size()) {
    std::cout << "[ERROR]: Expected to be passed a FHiCL configuration file "
                 "with the --fhicl option."
              << std::endl;
    SayUsage(argv);
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(fhicl_config);

  fhicl::ParameterSet flux_uncert_config =
      ps.get<fhicl::ParameterSet>("FluxUncertainty");

  std::vector<std::vector<double>> AllRelativeTweaks;

  for (fhicl::ParameterSet const &twk_ps :
       flux_uncert_config.get<std::vector<std::string>>("Tweaks")) {

    std::cout << "[INFO]: Building tweak: " << twk_ps.get<std::string>("Name")
              << std::endl;

    std::vector<double> flux_pred_nom;
    std::string NominalInputFile = twk_ps.get<std::string>("Nominal.InputFile");
    std::string NominalHistName =
        twk_ps.get<std::string>("Nominal.InputHistName");

    for (std::string const &pred_config :
         flux_uncert_config.get<std::vector<std::string>>("Configurations")) {
      std::cout << "[INFO]:\tFor configuration: " << pred_config << std::endl;
      NominalInputFile = str_replace(NominalInputFile, "%C", pred_config);
      NominalHistName = str_replace(NominalHistName, "%C", pred_config);

      for (std::string const &species :
           flux_uncert_config.get<std::vector<std::string>>("Species")) {
        std::cout << "[INFO]:\t\tFor species: " << species << std::endl;

        NominalInputFile = str_replace(NominalInputFile, "%S", species);
        NominalHistName = str_replace(NominalHistName, "%S", species);

        std::unique_ptr<TH1> nom_hist =
            GetHistogram_uptr<TH1>(NominalInputFile, NominalHistName);
        Mergestdvector(flux_pred_nom, Getstdvector(nom_hist.get()));
      } // End species -- Nominal
    }   // End configurations -- Nominal

    std::string InputFile = twk_ps.get<std::string>("InputFile");
    std::string VariedHistName = twk_ps.get<std::string>("VariedHistName");

    size_t NThrows = twk_ps.get<size_t>("NThrows", 1);
    for (size_t t_it = 0; t_it < NThrows; ++t_it) {
      std::vector<double> flux_pred_i;

      InputFile = str_replace(InputFile, "%i", std::to_string(t_it));
      VariedHistName = str_replace(VariedHistName, "%i", std::to_string(t_it));

      for (std::string const &pred_config :
           flux_uncert_config.get<std::vector<std::string>>("Configurations")) {
        std::cout << "[INFO]:\tFor configuration: " << pred_config << std::endl;

        InputFile = str_replace(InputFile, "%C", pred_config);
        VariedHistName = str_replace(VariedHistName, "%C", pred_config);

        for (std::string const &species :
             flux_uncert_config.get<std::vector<std::string>>("Species")) {
          std::cout << "[INFO]:\t\tFor species: " << species << std::endl;

          InputFile = str_replace(InputFile, "%S", species);
          VariedHistName = str_replace(VariedHistName, "%S", species);
          std::cout << "[INFO]:\t\tReading hist:\"" << VariedHistName
                    << "\" from file \"" << InputFile << "\"" << std::endl;

          std::unique_ptr<TH1> var_hist =
              GetHistogram_uptr<TH1>(InputFile, VariedHistName);
          Mergestdvector(flux_pred_i, Getstdvector(var_hist.get()));
        }
      }

      if (AllRelativeTweaks.size() &&
          (AllRelativeTweaks.back().size() != flux_pred_i.size())) {
        std::cout << "[ERROR]: Mismatched flux vector length: "
                  << AllRelativeTweaks.back().size()
                  << " != " << flux_pred_i.size() << std::endl;
        throw;
      }

      if (flux_pred_i.size() != flux_pred_nom.size()) {
        std::cout
            << "[ERROR]: Mismatched nominal flux vector and throw length: "
            << flux_pred_i.size() << " != " << flux_pred_nom.size()
            << std::endl;
        throw;
      }

      AllRelativeTweaks.push_back(std::move(flux_pred_i));
      for (size_t fbin_i = 0; fbin_i < flux_pred_nom.size(); ++fbin_i) {
        AllRelativeTweaks.back()[fbin_i] /= flux_pred_nom[fbin_i];
      }
    }
  }

  if (!AllRelativeTweaks.size()) {
    std::cout << "[WARN]: Found no configured flux tweaks" << std::endl;
    return 1;
  }

  CovarianceBuilder cb(AllRelativeTweaks.back().size());

  for (std::vector<double> const &flux_tweak : AllRelativeTweaks) {
    cb.AddThrow_MeanCalc(flux_tweak.data());
  }

  cb.FinalizeMeanCalc();

  for (std::vector<double> const &flux_tweak : AllRelativeTweaks) {
    cb.AddThrow_CovMatCalc(flux_tweak.data());
  }

  cb.FinalizeCovMatCalc();

  std::string OutputFile =
      flux_uncert_config.get<std::string>("OutputFile", "");
  bool UPDATEOutputFile =
      !flux_uncert_config.get<bool>("RecreateOutputFile", false);

  std::unique_ptr<TFile> oupF = std::unique_ptr<TFile>(
      CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE"));

  TMatrixD covmat = GetTMatrixD(cb.GetCovMatrix());
  TMatrixD corrmat = GetTMatrixD(cb.GetCorrMatrix());

  covmat.Write("covmat");
  corrmat.Write("corrmat");

  oupF->Write();
}
