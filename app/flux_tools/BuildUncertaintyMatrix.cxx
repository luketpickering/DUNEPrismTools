#include "GetUsage.hxx"

#include "CovarianceHelper.hxx"
#include "VariationBuilders.hxx"

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

  bool CovMatInitialized = false;
  Eigen::MatrixXd FullCovarianceMatrix;

  size_t MaxNDOF = 0;
  for (fhicl::ParameterSet twk_ps :
       flux_uncert_config.get<std::vector<std::string>>("Tweaks")) {
    for (std::string const &Conf : Configurations) {
      twk_ps.put(std::string("FluxSlicesDescriptor_") + Conf,
                 flux_uncert_config.get<std::string>(
                     std::string("FluxSlicesDescriptor_") + Conf, ""));
    }
    twk_ps.put("Species", Species);
    twk_ps.put("Configurations", Configurations);
    std::unique_ptr<VariationBuilder> varb = GetVariationBuilder(twk_ps);

    if (oupF) {
      varb->SetDiagnosticDirectory(oupF.get());
    }

    varb->Process();
    MaxNDOF += varb->GetNDOF();

    if (!CovMatInitialized) {
      FullCovarianceMatrix = varb->GetCovarianceComponent();
      CovMatInitialized = true;
    } else {
      FullCovarianceMatrix += varb->GetCovarianceComponent();
    }

    if (!NominalHistogramSet.size()) {
      NominalHistogramSet = varb->GetNominalHistograms();
    }
  }

  size_t NEigvals = flux_uncert_config.get<size_t>("num_eigenvalues", 10);

  NEigvals = std::min(NEigvals, MaxNDOF);

  std::cout << "[INFO]: Attempting to invert with " << MaxNDOF << " NDOF ..."
            << std::endl;
  EigenvalueHelper eh;
  try {
    eh.ComputeFromMatrix(FullCovarianceMatrix,
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
  std::cout << "[INFO]: Found " << NEigvals << " eigenvalues and eigenvectors."
            << std::endl;
  double KeptVar = 0;
  for (size_t e_it = 0; e_it < NEigvals; ++e_it) {
    KeptVar += eh.EigVals[e_it];
  }
  std::cout << "[INFO]: Minimum kept variance: "
            << (100.0 * KeptVar /
                (KeptVar +
                 eh.EigVals[NEigvals - 1] * double(MaxNDOF - NEigvals)))
            << "%" << std::endl;

  Eigen::MatrixXd Tweaks = eh.GetEffectiveParameterVectors();

  if (oupF) {

    TH1D *Evs = new TH1D("eigenvalues", ";Eigen value index;Magnitude",
                         NEigvals, 0, NEigvals);
    Evs->SetDirectory(oupF.get());

    FillHistFromEigenVector(Evs, eh.EigVals);

    TDirectory *td;
    if (!SaveCAFAnaFormat) {
      td = oupF->mkdir("EffectiveFluxParameters");
    }
    for (int tw_it = 0; tw_it < Tweaks.cols(); ++tw_it) {
      if (!SaveCAFAnaFormat) {

        TDirectory *t_id =
            td->mkdir((std::string("param_") + std::to_string(tw_it)).c_str());

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
        dir_caf_ss << "syst" << tw_it;
        TDirectory *td_caf = oupF->mkdir(dir_caf_ss.str().c_str());

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
