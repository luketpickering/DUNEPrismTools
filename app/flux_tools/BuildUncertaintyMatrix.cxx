#include "GetUsage.hxx"

#include "CovarianceHelper.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TH1D.h"
#include "TMatrixD.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#include <SymEigsSolver.h>

#include <cfloat>
#include <cmath>
#include <memory>
#include <string>

// #define DEBUG_BUILDUCERTMATRIX

const char *show_classification(double x) {
  switch (std::fpclassify(x)) {
  case FP_INFINITE:
    return "Inf";
  case FP_NAN:
    return "NaN";
  case FP_NORMAL:
    return "normal";
  case FP_SUBNORMAL:
    return "subnormal";
  case FP_ZERO:
    return "zero";
  default:
    return "unknown";
  }
}

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

  std::string flux_slice_descriptor =
      flux_uncert_config.get<std::string>("FluxSlicesDescriptor", "");

  std::vector<std::vector<double>> AllRelativeTweaks;
  std::vector<std::unique_ptr<TH1>> NominalHistograms;

  for (fhicl::ParameterSet const &twk_ps :
       flux_uncert_config.get<std::vector<std::string>>("Tweaks")) {

    std::string twk_name = twk_ps.get<std::string>("Name");
    std::cout << "[INFO]: Building tweak: " << twk_name << std::endl;

    std::vector<double> flux_pred_nom;
    std::string NominalInputFile = twk_ps.get<std::string>("Nominal.InputFile");
    std::string NominalHistName =
        twk_ps.get<std::string>("Nominal.InputHistName");

    for (std::string const &pred_config :
         flux_uncert_config.get<std::vector<std::string>>("Configurations")) {
#ifdef DEBUG_BUILDUCERTMATRIX
      std::cout << "[INFO]:\tFor configuration: " << pred_config << std::endl;
#endif
      NominalInputFile = str_replace(NominalInputFile, "%C", pred_config);
      NominalHistName = str_replace(NominalHistName, "%C", pred_config);

      for (std::string const &species :
           flux_uncert_config.get<std::vector<std::string>>("Species")) {
#ifdef DEBUG_BUILDUCERTMATRIX
        std::cout << "[INFO]:\t\tFor species: " << species << std::endl;
#endif
        NominalInputFile = str_replace(NominalInputFile, "%S", species);
        NominalHistName = str_replace(NominalHistName, "%S", species);

        if (flux_slice_descriptor.size()) {
          std::vector<std::pair<double, double>> XRanges =
              BuildRangesList(flux_slice_descriptor);
          std::unique_ptr<TH2D> flux2D =
              GetHistogram_uptr<TH2D>(NominalInputFile, NominalHistName);
          std::vector<std::unique_ptr<TH1D>> nom_hists =
              MergeSplitTH2D(flux2D, true, XRanges);
          Mergestdvector(flux_pred_nom, Getstdvector(nom_hists));
          size_t h_ctr = 0;
          for (auto &&h : nom_hists) {
            h->SetName(
                (twk_name + "_Nom_slice_" + std::to_string(h_ctr++)).c_str());
            NominalHistograms.push_back(std::move(h));
          }
        } else {
          std::unique_ptr<TH1> nom_hist =
              GetHistogram_uptr<TH1>(NominalInputFile, NominalHistName);
          Mergestdvector(flux_pred_nom, Getstdvector(nom_hist.get()));
          nom_hist->SetName((twk_name + "_Nom").c_str());
          NominalHistograms.push_back(std::move(nom_hist));
        }
      } // End species -- Nominal
    }   // End configurations -- Nominal

    std::string InputFile_template = twk_ps.get<std::string>("InputFile");
    std::string VariedHistName_template =
        twk_ps.get<std::string>("VariedHistName");

    size_t NThrows = twk_ps.get<size_t>("NThrows", 1);
    size_t NThrowSkip = twk_ps.get<size_t>("NThrowSkip", 0);
    for (size_t t_it = NThrowSkip; t_it < NThrows; ++t_it) {
      std::vector<double> flux_pred_i;

      std::string InputFile_i =
          str_replace(InputFile_template, "%i", std::to_string(t_it));
      std::string VariedHistName_i =
          str_replace(VariedHistName_template, "%i", std::to_string(t_it));

      for (std::string const &pred_config :
           flux_uncert_config.get<std::vector<std::string>>("Configurations")) {
#ifdef DEBUG_BUILDUCERTMATRIX
        std::cout << "[INFO]:\tFor configuration: " << pred_config << std::endl;
#endif
        std::string InputFile_i_c = str_replace(InputFile_i, "%C", pred_config);
        std::string VariedHistName_i_c =
            str_replace(VariedHistName_i, "%C", pred_config);

        for (std::string const &species :
             flux_uncert_config.get<std::vector<std::string>>("Species")) {
#ifdef DEBUG_BUILDUCERTMATRIX
          std::cout << "[INFO]:\t\tFor species: " << species << std::endl;
#endif
          std::string InputFile = str_replace(InputFile_i_c, "%S", species);
          std::string VariedHistName =
              str_replace(VariedHistName_i_c, "%S", species);

#ifdef DEBUG_BUILDUCERTMATRIX
          std::cout << "[INFO]:\t\tReading hist:\"" << VariedHistName
                    << "\" from file \"" << InputFile << "\"" << std::endl;
#endif

          if (flux_slice_descriptor.size()) {
            std::vector<std::pair<double, double>> XRanges =
                BuildRangesList(flux_slice_descriptor);
            std::unique_ptr<TH2D> flux2D =
                GetHistogram_uptr<TH2D>(InputFile, VariedHistName);
            Mergestdvector(flux_pred_i,
                           Getstdvector(MergeSplitTH2D(flux2D, true, XRanges)));
          } else {
            std::unique_ptr<TH1> var_hist =
                GetHistogram_uptr<TH1>(InputFile, VariedHistName);
            Mergestdvector(flux_pred_i, Getstdvector(var_hist.get()));
          };
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
            << flux_pred_nom.size() << " != " << flux_pred_i.size()
            << std::endl;
        throw;
      }

      AllRelativeTweaks.push_back(std::move(flux_pred_i));
      for (size_t fbin_i = 0; fbin_i < flux_pred_nom.size(); ++fbin_i) {
        if (!std::isnormal(AllRelativeTweaks.back()[fbin_i]) ||
            !std::isnormal(flux_pred_nom[fbin_i])) {
          AllRelativeTweaks.back()[fbin_i] = 1;
        } else {
          if (!std::isnormal(AllRelativeTweaks.back()[fbin_i] /
                             flux_pred_nom[fbin_i])) {
            std::cout << "[ERROR]: For flux prediction " << t_it << " for "
                      << twk_name << ", bin " << fbin_i << " was non-normal("
                      << show_classification(AllRelativeTweaks.back()[fbin_i] /
                                             flux_pred_nom[fbin_i])
                      << "): " << AllRelativeTweaks.back()[fbin_i] << "/"
                      << flux_pred_nom[fbin_i] << std::endl;
            throw;
          }
          AllRelativeTweaks.back()[fbin_i] /= flux_pred_nom[fbin_i];
        }
      }
    }
  }

  if (!AllRelativeTweaks.size()) {
    std::cout << "[WARN]: Found no configured flux tweaks" << std::endl;
    return 1;
  }

  CovarianceBuilder cb(AllRelativeTweaks.back().size());

#ifdef DEBUG_BUILDUCERTMATRIX
  size_t t_it = 0;
#endif
  std::cout << "[INFO]: Building " << AllRelativeTweaks.back().size() << "x"
            << AllRelativeTweaks.back().size() << " covariance matrix."
            << std::endl;
  std::cout << "[INFO]: Running mean calculation..." << std::endl;
  for (std::vector<double> const &flux_tweak : AllRelativeTweaks) {
#ifdef DEBUG_BUILDUCERTMATRIX
    std::cout << "[INFO]: Adding mean tweak " << t_it++ << " with "
              << flux_tweak.size() << " entries " << std::endl;
#endif
    cb.AddThrow_MeanCalc(flux_tweak.data());
  }

  cb.FinalizeMeanCalc();
  std::cout << "[INFO]: Done" << std::endl;

#ifdef DEBUG_BUILDUCERTMATRIX
  t_it = 0;
#endif
  std::cout << "[INFO]: Running covmat calculation..." << std::endl;
  for (std::vector<double> const &flux_tweak : AllRelativeTweaks) {
#ifdef DEBUG_BUILDUCERTMATRIX
    std::cout << "[INFO]: Adding covmat tweak " << t_it++ << " with "
              << flux_tweak.size() << " entries " << std::endl;
#endif
    cb.AddThrow_CovMatCalc(flux_tweak.data());
  }
  cb.FinalizeCovMatCalc();
  std::cout << "[INFO]: Done" << std::endl;

  if (flux_uncert_config.get<bool>("use_Eigen")) {
    std::cout << "[INFO]: Attempting decomposition with Eigen" << std::endl;
    Eigen::EigenSolver<Eigen::MatrixXd> es(cb.GetCovMatrix());
    std::cout << "[INFO]: Done" << std::endl;

    Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType evs = es.eigenvalues();

    double ev_tot = 0;
    for (Int_t i = 0; i < evs.size(); ++i) {
      ev_tot += evs(i).real();
    }

    double ev_run = 0;
    for (Int_t i = 0; i < evs.size(); ++i) {
      ev_run += evs(i).real();
      std::cout << "EV[" << i << "] = " << evs(i) << ", All Evs < " << i
                << " contain " << (ev_run / ev_tot) * 100.0
                << "% of the variance." << std::endl;
    }
  } else {
    size_t NEVs = flux_uncert_config.get<size_t>("num_eigenvalues");
    std::cout
        << "[INFO]: Attempting decomposition with Spectra, looking for top "
        << NEVs << " eigenvalues." << std::endl;

    Spectra::DenseSymMatProd<double> op(cb.GetCovMatrix());
    Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                           Spectra::DenseSymMatProd<double>>
        eigs(&op, NEVs, std::min(2 * NEVs, AllRelativeTweaks.back().size()));
    // Initialize and compute
    eigs.init();
    eigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    if (eigs.info() == Spectra::SUCCESSFUL) {
      evalues = eigs.eigenvalues();
    }

    std::cout << "Eigenvalues found:\n" << evalues.size() << std::endl;

    for (Int_t i = 0; i < evalues.size(); ++i) {
      std::cout << "EV[" << i << "] = " << evalues(i) << std::endl;
    }
  }

  std::string OutputFile = flux_uncert_config.get<std::string>("OutputFile");
  bool UPDATEOutputFile =
      !flux_uncert_config.get<bool>("RecreateOutputFile", false);

  std::unique_ptr<TFile> oupF = std::unique_ptr<TFile>(
      CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE"));

  std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
  std::unique_ptr<TMatrixD> corrmat = GetTMatrixD(cb.GetCorrMatrix());

  covmat->Write("covmat");
  corrmat->Write("corrmat");

  for (auto &h : NominalHistograms) {
    h->Write(h->GetName());
  }

  oupF->Write();
}
