#include "GetUsage.hxx"

#include "CovarianceHelper.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TH1D.h"
#include "TMatrixD.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "Eigen/Core"

#include <cfloat>
#include <chrono>
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

  std::string OutputFile =
      flux_uncert_config.get<std::string>("OutputFile", "");
  std::unique_ptr<TFile> oupF(nullptr);
  if (OutputFile.size()) {
    bool UPDATEOutputFile =
        !flux_uncert_config.get<bool>("RecreateOutputFile", false);

    oupF = std::unique_ptr<TFile>(
        CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE"));
  }

  std::string flux_slice_descriptor =
      flux_uncert_config.get<std::string>("FluxSlicesDescriptor", "");

  std::vector<std::unique_ptr<TH1>> NominalHistogramSet;

  bool CovMatInitialized = false;
  Eigen::MatrixXd FullCovarianceMatrix;

  for (fhicl::ParameterSet const &twk_ps :
       flux_uncert_config.get<std::vector<std::string>>("Tweaks")) {

    bool BuildNomHistoSet = !NominalHistogramSet.size();
    bool dump_diagnostics = twk_ps.get<bool>("dump_diagnostics", false);

    std::vector<std::vector<double>> RelativeTweaks;
    std::vector<std::vector<double>>
        diags_Predictions; // only used for diagnostics

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
          if (BuildNomHistoSet) {
            if (nom_hists.size() > 1) {
              NominalHistogramSet.push_back(ReMergeSplitTH2D(
                  nom_hists, XRanges, flux2D->GetName(), flux2D->GetTitle()));
            } else if (nom_hists.size()) {
              NominalHistogramSet.push_back(std::move(nom_hists.front()));
            }
          }
        } else {
          std::unique_ptr<TH1> nom_hist =
              GetHistogram_uptr<TH1>(NominalInputFile, NominalHistName);
          Mergestdvector(flux_pred_nom, Getstdvector(nom_hist.get()));
          nom_hist->SetName((twk_name + "_Nom").c_str());
          if (BuildNomHistoSet) {
            NominalHistogramSet.push_back(std::move(nom_hist));
          }
        }

        if (BuildNomHistoSet) {
          NominalHistogramSet.back()->SetName(
              (pred_config + "_" + species).c_str());
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

      if (RelativeTweaks.size() &&
          (RelativeTweaks.back().size() != flux_pred_i.size())) {
        std::cout << "[ERROR]: Mismatched flux vector length: "
                  << RelativeTweaks.back().size()
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

      if (dump_diagnostics) {
        diags_Predictions.push_back(flux_pred_i);
      }

      RelativeTweaks.push_back(std::move(flux_pred_i));
      for (size_t fbin_i = 0; fbin_i < flux_pred_nom.size(); ++fbin_i) {
        if (!std::isnormal(RelativeTweaks.back()[fbin_i]) ||
            !std::isnormal(flux_pred_nom[fbin_i])) {
          RelativeTweaks.back()[fbin_i] = 1;
        } else {
          if (!std::isnormal(RelativeTweaks.back()[fbin_i] /
                             flux_pred_nom[fbin_i])) {
            std::cout << "[ERROR]: For flux prediction " << t_it << " for "
                      << twk_name << ", bin " << fbin_i << " was non-normal("
                      << show_classification(RelativeTweaks.back()[fbin_i] /
                                             flux_pred_nom[fbin_i])
                      << "): " << RelativeTweaks.back()[fbin_i] << "/"
                      << flux_pred_nom[fbin_i] << std::endl;
            throw;
          }
          RelativeTweaks.back()[fbin_i] =
              1 - (RelativeTweaks.back()[fbin_i] / flux_pred_nom[fbin_i]);
        }
      }
    }
    if (!RelativeTweaks.size()) {
      std::cout << "[WARN]: Found no configured flux tweaks" << std::endl;
      return 1;
    }

    CovarianceBuilder cb(RelativeTweaks);

    if (!CovMatInitialized) {
      FullCovarianceMatrix = cb.GetCovMatrix();
      CovMatInitialized = true;
      std::cout << "[INFO]: Built first covariance component, size: "
                << cb.GetCovMatrix().rows() << " x " << cb.GetCovMatrix().rows()
                << std::endl;
    } else if (FullCovarianceMatrix.rows() == cb.GetCovMatrix().rows()) {
      FullCovarianceMatrix += cb.GetCovMatrix();
    } else {
      std::cout << "[ERROR]: Attempting to add covariance component of size "
                << cb.GetCovMatrix().rows() << " x " << cb.GetCovMatrix().rows()
                << " to full covariance of size: "
                << FullCovarianceMatrix.rows() << " x "
                << FullCovarianceMatrix.rows() << std::endl;
      return 1;
    }

    if (oupF && dump_diagnostics) {
      TDirectory *td =
          oupF->mkdir((std::string("Diagnostics_") + twk_name).c_str());

      std::vector<std::unique_ptr<TH1>> nom_histograms =
          CloneHistVector(NominalHistogramSet, "_nominal_pred");

      FillHistFromstdvector(nom_histograms, flux_pred_nom);

      for (std::unique_ptr<TH1> &h : nom_histograms) {
        h->SetDirectory(td);
        h.release(); // Let root look after the histo again.
      }

      std::vector<std::unique_ptr<TH1>> mean_pred_histograms =
          CloneHistVector(NominalHistogramSet, "_mean_pred");

      std::vector<double> mean_pred_vector;
      std::vector<double> stddev_pred_vector;
      for (int i = 0; i < cb.NRows; ++i) {
        mean_pred_vector.push_back((1 + cb.GetMeanVector()[i]) *
                                   flux_pred_nom[i]);
        stddev_pred_vector.push_back(cb.GetStdDevVector()[i] *
                                     flux_pred_nom[i]);
      }

      FillHistFromstdvector(mean_pred_histograms, mean_pred_vector, 0,
                            stddev_pred_vector);

      for (std::unique_ptr<TH1> &h : mean_pred_histograms) {
        h->SetDirectory(td);
        h.release(); // Let root look after the histo again.
      }

      std::vector<std::unique_ptr<TH1>> stddev_histograms =
          CloneHistVector(NominalHistogramSet, "_stddev");

      TH1D *mean_std_dev = new TH1D("mean_std_dev", "", mean_pred_vector.size(),
                                    0, mean_pred_vector.size());
      mean_std_dev->SetDirectory(td);

      std::vector<double> mean_vector;
      std::vector<double> stddev_vector;
      for (int i = 0; i < cb.NRows; ++i) {
        stddev_vector.push_back(cb.GetStdDevVector()[i]);
        mean_vector.push_back(cb.GetMeanVector()[i]);
      }

      FillHistFromstdvector(stddev_histograms, stddev_vector);
      FillHistFromstdvector(mean_std_dev, mean_vector, 0, stddev_vector);

      for (std::unique_ptr<TH1> &h : stddev_histograms) {
        h->SetDirectory(td);
        h.release(); // Let root look after the histo again.
      }

      std::vector<std::unique_ptr<TH1>> mean_fractional_histograms =
          CloneHistVector(NominalHistogramSet, "_mean_fractional");

      std::vector<double> mean_fractional_vector;
      for (int i = 0; i < cb.NRows; ++i) {
        mean_fractional_vector.push_back(cb.GetMeanVector()[i]);
      }

      FillHistFromstdvector(mean_fractional_histograms, mean_fractional_vector);

      for (std::unique_ptr<TH1> &h : mean_fractional_histograms) {
        h->SetDirectory(td);
        h.release(); // Let root look after the histo again.
      }

      TH2D *variation_distributions =
          new TH2D("ThrowVariations", ";Bin content;Bin number;NThrows",
                   stddev_vector.size(), 0, stddev_vector.size(), 100, -1, 1);
      variation_distributions->SetDirectory(td);

      for (size_t tw_it = 0; tw_it < RelativeTweaks.size(); ++tw_it) {
        TDirectory *t_id =
            td->mkdir((std::string("_tweak_") + std::to_string(tw_it)).c_str());

        std::vector<std::unique_ptr<TH1>> pred_histograms = CloneHistVector(
            NominalHistogramSet, std::string("_pred_") + std::to_string(tw_it));

        FillHistFromstdvector(pred_histograms, diags_Predictions[tw_it]);

        for (std::unique_ptr<TH1> &h : pred_histograms) {
          h->SetDirectory(t_id);
          h.release(); // Let root look after the histo again.
        }

        std::vector<std::unique_ptr<TH1>> tweak_histograms =
            CloneHistVector(NominalHistogramSet, std::string("_fractional_") +
                                                     std::to_string(tw_it));

        FillHistFromstdvector(tweak_histograms, RelativeTweaks[tw_it]);

        for (std::unique_ptr<TH1> &h : tweak_histograms) {
          h->SetDirectory(t_id);
          h.release(); // Let root look after the histo again.
        }

        for (size_t bin_it = 0; bin_it < RelativeTweaks[tw_it].size();
             ++bin_it) {
          variation_distributions->Fill(bin_it, RelativeTweaks[tw_it][bin_it]);
        }
      }

      TDirectory *md = td->mkdir("Matrices");
      md->cd();
      std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
      std::unique_ptr<TMatrixD> corrmat =
          GetTMatrixD(CovToCorr(cb.GetCovMatrix()));

      covmat->Write("covmat");
      corrmat->Write("corrmat");
    }
  }

  size_t NEigvals = flux_uncert_config.get<size_t>("num_eigenvalues", 10);

  EigenvalueHelper eh;
  eh.ComputeFromMatrix(FullCovarianceMatrix,
                       flux_uncert_config.get<bool>("use_Spectra", true),
                       NEigvals);

  NEigvals = eh.EigVals.size();

  Eigen::MatrixXd Tweaks = eh.GetEffectiveParameterVectors();

  if (oupF) {

    TH1D *Evs = new TH1D("eigenvalues", ";Eigen value index;Magnitude",
                         NEigvals, 0, NEigvals);
    Evs->SetDirectory(oupF.get());

    FillHistFromEigenVector(Evs, eh.EigVals);

    TDirectory *td = oupF->mkdir("EffectiveFluxParameters");
    for (int tw_it = 0; tw_it < Tweaks.rows(); ++tw_it) {
      TDirectory *t_id =
          td->mkdir((std::string("param_") + std::to_string(tw_it)).c_str());

      std::vector<std::unique_ptr<TH1>> tweak_histograms =
          CloneHistVector(NominalHistogramSet);

      FillHistFromEigenVector(tweak_histograms, Tweaks.row(tw_it));

      for (std::unique_ptr<TH1> &h : tweak_histograms) {
        h->SetDirectory(t_id);
        h.release(); // Let root look after the histo again.
      }
    }

    TDirectory *n_d = td->mkdir("nominal");
    for (auto &h : NominalHistogramSet) {
      h->SetDirectory(n_d);
      h.release();
    }

    if (flux_uncert_config.get<bool>("WriteMatrices", false)) {
      TDirectory *md = oupF->mkdir("Matrices");
      md->cd();
      std::unique_ptr<TMatrixD> covmat = GetTMatrixD(FullCovarianceMatrix);
      std::unique_ptr<TMatrixD> corrmat =
          GetTMatrixD(CovToCorr(FullCovarianceMatrix));

      covmat->Write("covmat");
      corrmat->Write("corrmat");
    }

    oupF->Write();
  }
}
