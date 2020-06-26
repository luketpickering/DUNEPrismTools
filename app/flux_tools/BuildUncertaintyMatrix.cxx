#include "GetUsage.hxx"

#include "CovarianceHelper.hxx"
#include "VariationBuilders.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TH1D.h"
#include "TList.h"
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

  bool SaveCAFAnaFormat =
      flux_uncert_config.get<bool>("SaveCAFAnaFormat", false);
  bool SaveTH1F = flux_uncert_config.get<bool>("SaveTH1F", false);

  std::string OutputFile =
      flux_uncert_config.get<std::string>("OutputFile", "");
  std::unique_ptr<TFile> oupF(nullptr);
  if (OutputFile.size()) {
    bool UPDATEOutputFile =
        !flux_uncert_config.get<bool>("RecreateOutputFile", false);

    oupF = std::unique_ptr<TFile>(
        CheckOpenFile(OutputFile, UPDATEOutputFile ? "UPDATE" : "RECREATE"));
  }

  std::vector<std::string> Species =
      flux_uncert_config.get<std::vector<std::string>>("Species");
  std::vector<std::string> Configurations =
      flux_uncert_config.get<std::vector<std::string>>("Configurations");

  std::vector<std::unique_ptr<TH1>> NominalHistogramSet;

  bool PCACovMatInitialized = false, FullCovMatInitialized = false;
  Eigen::MatrixXd FullCovarianceMatrix;
  Eigen::MatrixXd PCACovarianceMatrix;

  size_t MaxNDOF = 0;
  std::map<std::string, std::vector<std::unique_ptr<TH1>>> NonPCAComponents;
  for (fhicl::ParameterSet twk_ps :
       flux_uncert_config.get<std::vector<std::string>>("Tweaks")) {
    for (std::string const &Conf : Configurations) {
      bool isjagged = flux_uncert_config.get<bool>(Conf + "_IsJagged", false);
      twk_ps.put(Conf + "_IsJagged", isjagged);
      size_t offaxisbin = flux_uncert_config.get<size_t>(
          Conf + "_OffAxisBin", std::numeric_limits<size_t>::max());
      twk_ps.put(Conf + "_OffAxisBin", offaxisbin);

      // for jagged it is a sequence
      if (isjagged) {
        std::vector<int> const &fluxdescriptor =
            flux_uncert_config.get<std::vector<int>>(
                Conf + "_FluxSlicesDescriptor", {});
        if (fluxdescriptor.size()) {
          twk_ps.put(Conf + "_FluxSlicesDescriptor", fluxdescriptor);
        }

        std::string fluxmerge =
            flux_uncert_config.get<std::string>(Conf + "_FluxSlicesMerge", "");
        if (fluxmerge.size()) {
          twk_ps.put(Conf + "_FluxSlicesMerge", fluxmerge);
        }
      }
      // for legacy TH2 it will be a string
      else {
        std::string fluxdescriptor = flux_uncert_config.get<std::string>(
            Conf + "_FluxSlicesDescriptor", "");
        if (fluxdescriptor.size()) {
          twk_ps.put(Conf + "_FluxSlicesDescriptor", fluxdescriptor);
        }
      }
    }
    twk_ps.put("Species", Species);
    twk_ps.put("Configurations", Configurations);
    std::unique_ptr<VariationBuilder> varb = GetVariationBuilder(twk_ps);

    if (!NominalHistogramSet.size()) {
      NominalHistogramSet = varb->GetNominalHistograms();
    }

    // We might want to not bother PCA-grinding some of the inputs
    if (!twk_ps.get<bool>("include_in_PCA", false)) {
      NonPCAComponents[twk_ps.get<std::string>("Name")] =
          varb->GetOneSigmaVariation();

      // Turn this 1sig variation into a covariance matrix and add it to the
      // running total covariance matrix
      std::vector<double> onesigtwk;
      for (auto &h : NonPCAComponents[twk_ps.get<std::string>("Name")]) {
        Mergestdvector(onesigtwk, Getstdvector(h.get()));
      }

      CovarianceBuilder cb(
          {
              onesigtwk,
          },
          true);

      auto CovarianceComponent = cb.GetCovMatrix();

      if (!FullCovMatInitialized) {
        FullCovarianceMatrix = CovarianceComponent;
        FullCovMatInitialized = true;
      } else {
        FullCovarianceMatrix += CovarianceComponent;
      }

      continue;
    }

    if (oupF) {
      varb->SetDiagnosticDirectory(oupF.get());
    }

    varb->Process();
    MaxNDOF += varb->GetNDOF();

    if (!FullCovMatInitialized) {
      FullCovarianceMatrix = varb->GetCovarianceComponent();
      FullCovMatInitialized = true;
    } else {
      FullCovarianceMatrix += varb->GetCovarianceComponent();
    }

    if (!PCACovMatInitialized) {
      PCACovarianceMatrix = varb->GetCovarianceComponent();
      PCACovMatInitialized = true;
    } else {
      PCACovarianceMatrix += varb->GetCovarianceComponent();
    }
  }

  size_t NEigvals = 0;
  Eigen::MatrixXd Tweaks = Eigen::MatrixXd::Zero(1, 1);
  EigenvalueHelper eh;
  if (PCACovMatInitialized) {
    NEigvals = flux_uncert_config.get<size_t>("num_eigenvalues", 10);

    NEigvals = std::min(NEigvals, MaxNDOF);

    std::cout << "[INFO]: Attempting to invert with " << MaxNDOF << " NDOF ..."
              << std::endl;
    try {
      eh.ComputeFromMatrix(PCACovarianceMatrix,
                           flux_uncert_config.get<bool>("use_Spectra", true),
                           NEigvals);
    } catch (std::runtime_error) {
      std::cout << "[ERROR]: Inversion failed." << std::endl;
      if (oupF) {
        oupF->Write();
        return 1;
      }
    }

    NEigvals = eh.EigVals.size();
    std::cout << "[INFO]: Found " << NEigvals
              << " eigenvalues and eigenvectors." << std::endl;
    double KeptVar = 0;
    for (size_t e_it = 0; e_it < NEigvals; ++e_it) {
      KeptVar += eh.EigVals[e_it];
    }
    std::cout << "[INFO]: Minimum kept variance: "
              << (100.0 * KeptVar /
                  (KeptVar +
                   eh.EigVals[NEigvals - 1] * double(MaxNDOF - NEigvals)))
              << "%" << std::endl;

    Tweaks = eh.GetEffectiveParameterVectors();

    bool ForceEffectiveParametersPositiveOnAxis = flux_uncert_config.get<bool>(
        "ForceEffectiveParametersPositiveOnAxis", false);
    if (ForceEffectiveParametersPositiveOnAxis) {
      // Only want to flip based on on-axis to preclude random noise at large OA
      // angles.
      size_t NBins = NominalHistogramSet.front()->GetNbinsX();

      for (int ev_it = 0; ev_it < Tweaks.cols(); ++ev_it) {
        double min = std::numeric_limits<double>::max();
        double max = -std::numeric_limits<double>::max();

        for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
          double r = Tweaks(bi_it, ev_it);
          min = std::min(r, min);
          max = std::max(r, max);
        }

        // std::cout << "For ev " << ev_it << " searched " << NBins
        //           << " bins, min = " << min << ", max = " << max
        //           << (((max < 0) || (std::abs(min) > std::abs(max)))
        //                   ? " flipping"
        //                   : " not flipping.")
        //           << std::endl;

        if ((max < 0) || (std::abs(min) > std::abs(max))) {
          Tweaks.col(ev_it) *= -1;
        }
      }
    }
  }

  if (oupF) {

    TDirectory *td = nullptr;
    if (!SaveCAFAnaFormat) {
      td = oupF->mkdir("FluxParameters");
    }

    std::vector<std::string> param_names(NonPCAComponents.size() + NEigvals);
    // Do non PCA first
    size_t nNonPCA = 0;
    for (auto &comp : NonPCAComponents) {
      if (!SaveCAFAnaFormat) {

        param_names[nNonPCA] = std::string("param_") + comp.first;
        TDirectory *t_id = td->mkdir(param_names[nNonPCA].c_str());
        for (std::unique_ptr<TH1> &h : comp.second) {
          if (SaveTH1F) {
            std::unique_ptr<TH1> h_f = THToF(h);
            h_f->SetDirectory(t_id);
            h_f.release();
          } else {
            h->SetDirectory(t_id);
            h.release();
          }
        }
      } else {
        std::stringstream dir_caf_ss("");
        dir_caf_ss << "syst" << nNonPCA;
        param_names[nNonPCA] = dir_caf_ss.str();

        TDirectory *td_caf = oupF->mkdir(param_names[nNonPCA].c_str());
        for (std::unique_ptr<TH1> &h : comp.second) {
          // Name translation
          std::string hname = h->GetName();

          std::string species;
          if (hname.find("numubar") != std::string::npos) {
            species = "numubar";
          } else if (hname.find("numu") != std::string::npos) {
            species = "numu";
          } else if (hname.find("nuebar") != std::string::npos) {
            species = "nuebar";
          } else if (hname.find("nue") != std::string::npos) {
            species = "nue";
          }
          std::string det;
          std::string beam_mode;
          if (hname.find("ND_nubar") != std::string::npos) {
            det = "ND";
            beam_mode = "RHC";
          } else if (hname.find("ND_nu") != std::string::npos) {
            det = "ND";
            beam_mode = "FHC";
          } else if (hname.find("FD_nubar") != std::string::npos) {
            det = "FD";
            beam_mode = "RHC";
          } else if (hname.find("FD_nu") != std::string::npos) {
            det = "FD";
            beam_mode = "FHC";
          }

          std::stringstream caf_name("");

          caf_name << det << "_" << species << "_" << beam_mode;
          h->SetName(caf_name.str().c_str());
          if (SaveTH1F) {
            std::unique_ptr<TH1> h_f = THToF(h);
            h_f->SetDirectory(td_caf);
            h_f.release();
          } else {
            h->SetDirectory(td_caf);
            h.release();
          }
        }
      }
      nNonPCA++;
    }

    if (PCACovMatInitialized) {

      TH1D *Evs = new TH1D("pca_eigenvalues", ";Eigen value index;Magnitude",
                           NEigvals, 0, NEigvals);
      Evs->SetDirectory(oupF.get());

      FillHistFromEigenVector(Evs, eh.EigVals);

      for (int tw_it = 0; tw_it < Tweaks.cols(); ++tw_it) {
        if (!SaveCAFAnaFormat) {

          param_names[nNonPCA + tw_it] =
              std::string("param_pca_") + std::to_string(tw_it);

          TDirectory *t_id = td->mkdir(param_names[nNonPCA + tw_it].c_str());

          std::vector<std::unique_ptr<TH1>> tweak_histograms =
              CloneHistVector(NominalHistogramSet);

          FillHistFromEigenVector(tweak_histograms, Tweaks.col(tw_it));

          for (std::unique_ptr<TH1> &h : tweak_histograms) {
            if (SaveTH1F) {
              std::unique_ptr<TH1> h_f = THToF(h);
              h_f->SetDirectory(t_id);
              h_f.release();
            } else {
              h->SetDirectory(t_id);
              h.release();
            }
          }
        } else {

          // Output with CAFAna naming scheme
          std::stringstream dir_caf_ss("");
          dir_caf_ss << "syst" << tw_it + nNonPCA;
          param_names[nNonPCA + tw_it] = dir_caf_ss.str();

          TDirectory *td_caf =
              oupF->mkdir(param_names[nNonPCA + tw_it].c_str());

          std::vector<std::unique_ptr<TH1>> tweak_histograms_caf =
              CloneHistVector(NominalHistogramSet);

          FillHistFromEigenVector(tweak_histograms_caf, Tweaks.col(tw_it));

          for (std::unique_ptr<TH1> &h : tweak_histograms_caf) {
            // Name translation
            std::string hname = h->GetName();

            std::string species;
            if (hname.find("numubar") != std::string::npos) {
              species = "numubar";
            } else if (hname.find("numu") != std::string::npos) {
              species = "numu";
            } else if (hname.find("nuebar") != std::string::npos) {
              species = "nuebar";
            } else if (hname.find("nue") != std::string::npos) {
              species = "nue";
            }
            std::string det;
            std::string beam_mode;
            if (hname.find("ND_nubar") != std::string::npos) {
              det = "ND";
              beam_mode = "RHC";
            } else if (hname.find("ND_nu") != std::string::npos) {
              det = "ND";
              beam_mode = "FHC";
            } else if (hname.find("FD_nubar") != std::string::npos) {
              det = "FD";
              beam_mode = "RHC";
            } else if (hname.find("FD_nu") != std::string::npos) {
              det = "FD";
              beam_mode = "FHC";
            }

            std::stringstream caf_name("");

            caf_name << det << "_" << species << "_" << beam_mode;
            h->SetName(caf_name.str().c_str());
            if (SaveTH1F) {
              std::unique_ptr<TH1> h_f = THToF(h);
              h_f->SetDirectory(td_caf);
              h_f.release();
            } else {
              h->SetDirectory(td_caf);
              h.release();
            }
          }
        }
      }
    }

    Eigen::VectorXd rtvar = FullCovarianceMatrix.diagonal().cwiseSqrt();
    std::vector<std::unique_ptr<TH1>> std_dev_histograms =
        CloneHistVector(NominalHistogramSet);
    FillHistFromEigenVector(std_dev_histograms, rtvar);

    TList param_names_t;

    for (auto &n : param_names) {
      param_names_t.AddLast(new TObjString(n.c_str()));
    }

    if (td) {
      td->WriteTObject(&param_names_t, "param_names");
    } else {
      oupF->WriteTObject(&param_names_t, "param_names");
    }

    TDirectory *e_d = oupF->mkdir("TotalError1D_total");
    for (auto &h : std_dev_histograms) {
      if (SaveTH1F) {
        std::unique_ptr<TH1> h_f = THToF(h);
        h_f->SetDirectory(e_d);
        h_f.release();
      } else {
        h->SetDirectory(e_d);
        h.release();
      }
    }

    if (PCACovMatInitialized) {
      Eigen::VectorXd rtvar_pca = PCACovarianceMatrix.diagonal().cwiseSqrt();

      std::vector<std::unique_ptr<TH1>> std_dev_pca_histograms =
          CloneHistVector(NominalHistogramSet);
      FillHistFromEigenVector(std_dev_pca_histograms, rtvar_pca);

      e_d = oupF->mkdir("TotalError1D_pcaonly");
      for (auto &h : std_dev_pca_histograms) {
        if (SaveTH1F) {
          std::unique_ptr<TH1> h_f = THToF(h);
          h_f->SetDirectory(e_d);
          h_f.release();
        } else {
          h->SetDirectory(e_d);
          h.release();
        }
      }
    }

    TDirectory *n_d = oupF->mkdir("nominal");
    for (auto &h : NominalHistogramSet) {
      if (SaveTH1F) {
        std::unique_ptr<TH1> h_f = THToF(h);
        h_f->SetDirectory(n_d);
        h_f.release();
      } else {
        h->SetDirectory(n_d);
        h.release();
      }
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
