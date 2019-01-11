#include "VariationBuilders.hxx"

#include "CovarianceHelper.hxx"

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

void VariationBuilder::BaseConfigure(fhicl::ParameterSet const &ps) {
  NDOF = 0;

  paramset = ps;
  DumpDiagnostics = paramset.get<bool>("dump_diagnostics", false);
  DumpMatrices = paramset.get<bool>("dump_matrices", false);
  Name = paramset.get<std::string>("Name");
  std::cout << "[INFO]: Building tweak: " << Name << std::endl;

  Configurations = paramset.get<std::vector<std::string>>("Configurations");
  Species = paramset.get<std::vector<std::string>>("Species");

  std::string NominalInputFile_template =
      paramset.get<std::string>("Nominal.InputFile");
  std::string NominalHistName_template =
      paramset.get<std::string>("Nominal.InputHistName");

  for (std::string const &pred_config : Configurations) {

    std::string FluxSliceDescriptor = paramset.get<std::string>(
        std::string("FluxSlicesDescriptor_") + pred_config, "");

    std::string NominalInputFile_conf =
        str_replace(NominalInputFile_template, "%C", pred_config);
    std::string NominalHistName_conf =
        str_replace(NominalHistName_template, "%C", pred_config);

    for (std::string const &species : Species) {

      std::string NominalInputFile =
          str_replace(NominalInputFile_conf, "%S", species);
      std::string NominalHistName =
          str_replace(NominalHistName_conf, "%S", species);

      if (FluxSliceDescriptor.size()) {
        std::vector<std::pair<double, double>> XRanges =
            BuildRangesList(FluxSliceDescriptor);
        std::unique_ptr<TH2D> flux2D =
            GetHistogram_uptr<TH2D>(NominalInputFile, NominalHistName);

        std::vector<std::unique_ptr<TH1D>> nom_hists =
            MergeSplitTH2D(flux2D, true, XRanges);
        Mergestdvector(NominalPrediction, Getstdvector(nom_hists));
        Mergestdvector(NominalPredictionError, Getstdvectorerror(nom_hists));

        if (nom_hists.size() > 1) {
          NominalHistogramSet.push_back(ReMergeSplitTH2D(
              nom_hists, XRanges, flux2D->GetName(), flux2D->GetTitle()));
        } else if (nom_hists.size()) {
          NominalHistogramSet.push_back(std::move(nom_hists.front()));
        }

      } else {
        std::unique_ptr<TH1> nom_hist =
            GetHistogram_uptr<TH1>(NominalInputFile, NominalHistName);
        Mergestdvector(NominalPrediction, Getstdvector(nom_hist.get()));
        Mergestdvector(NominalPredictionError,
                       Getstdvectorerror(nom_hist.get()));

        nom_hist->SetName((Name + "_Nom").c_str());
        NominalHistogramSet.push_back(std::move(nom_hist));
      }

      NominalHistogramSet.back()->SetName(
          (pred_config + "_" + species).c_str());

    } // End species -- Nominal
  }   // End configurations -- Nominal
}

void ThrownVariations::Configure(fhicl::ParameterSet const &ps) {
  BaseConfigure(ps);
}
void ThrownVariations::Process() {

  fhicl::ParameterSet VariedHist_ps =
      paramset.get<fhicl::ParameterSet>("Varied");
  std::string InputFile_template = VariedHist_ps.get<std::string>("InputFile");
  std::string VariedHistName_template =
      VariedHist_ps.get<std::string>("InputHistName");

  size_t NThrows = paramset.get<size_t>("NThrows");
  size_t NThrowSkip = paramset.get<size_t>("NThrowSkip", 0);

  NDOF = NThrows - NThrowSkip;

  std::cout << "[INFO]: Processing thrown source: " << Name << ", with "
            << NThrows << " throws." << std::endl;

  for (size_t t_it = NThrowSkip; t_it < NThrows; ++t_it) {
    std::vector<double> flux_pred_i;

    std::string InputFile_i =
        str_replace(InputFile_template, "%i", std::to_string(t_it));

    std::string VariedHistName_i =
        str_replace(VariedHistName_template, "%i", std::to_string(t_it));

    for (std::string const &pred_config : Configurations) {

      std::string FluxSliceDescriptor = paramset.get<std::string>(
          std::string("FluxSlicesDescriptor_") + pred_config, "");

      std::string InputFile_i_c = str_replace(InputFile_i, "%C", pred_config);
      std::string VariedHistName_i_c =
          str_replace(VariedHistName_i, "%C", pred_config);

      for (std::string const &species : Species) {

        std::string InputFile = str_replace(InputFile_i_c, "%S", species);
        std::string VariedHistName =
            str_replace(VariedHistName_i_c, "%S", species);

        if (FluxSliceDescriptor.size()) {
          std::vector<std::pair<double, double>> XRanges =
              BuildRangesList(FluxSliceDescriptor);

          std::unique_ptr<TH2D> flux2D =
              GetHistogram_uptr<TH2D>(InputFile, VariedHistName);

          Mergestdvector(flux_pred_i,
                         Getstdvector(MergeSplitTH2D(flux2D, true, XRanges)));
        } else {
          std::unique_ptr<TH1> var_hist =
              GetHistogram_uptr<TH1>(InputFile, VariedHistName);

          Mergestdvector(flux_pred_i, Getstdvector(var_hist.get()));
        }
      }
    }

    if (RelativeTweaks.size() &&
        (RelativeTweaks.back().size() != flux_pred_i.size())) {
      std::cout << "[ERROR]: Mismatched flux vector length: "
                << RelativeTweaks.back().size() << " != " << flux_pred_i.size()
                << std::endl;
      throw;
    }

    if (flux_pred_i.size() != NominalPrediction.size()) {
      std::cout << "[ERROR]: Mismatched nominal flux vector and throw length: "
                << NominalPrediction.size() << " != " << flux_pred_i.size()
                << std::endl;
      throw;
    }

    if (DumpDiagnostics) {
      diags_Predictions.push_back(flux_pred_i);
    }

    RelativeTweaks.push_back(std::move(flux_pred_i));
    for (size_t fbin_i = 0; fbin_i < NominalPrediction.size(); ++fbin_i) {
      if (!std::isnormal(RelativeTweaks.back()[fbin_i]) ||
          !std::isnormal(NominalPrediction[fbin_i])) {
        RelativeTweaks.back()[fbin_i] = 0;
      } else {
        if (!std::isnormal(RelativeTweaks.back()[fbin_i] /
                           NominalPrediction[fbin_i])) {
          std::cout << "[ERROR]: For flux prediction " << t_it << " for "
                    << Name << ", bin " << fbin_i << " was non-normal("
                    << show_classification(RelativeTweaks.back()[fbin_i] /
                                           NominalPrediction[fbin_i])
                    << "): " << RelativeTweaks.back()[fbin_i] << "/"
                    << NominalPrediction[fbin_i] << std::endl;
          throw;
        }
        RelativeTweaks.back()[fbin_i] =
            (RelativeTweaks.back()[fbin_i] / NominalPrediction[fbin_i]) - 1;
      }
    }
  }

  if (!RelativeTweaks.size()) {
    std::cout << "[WARN]: Found no configured flux tweaks" << std::endl;
    return;
  }

  std::cout << "[INFO]: Building " << RelativeTweaks.front().size() << "x"
            << RelativeTweaks.front().size() << " covariance matrix."
            << std::endl;

  CovarianceBuilder cb(RelativeTweaks);

  CovarianceComponent = cb.GetCovMatrix();

  if (!(diagdir && DumpDiagnostics)) {
    return;
  }

  TDirectory *td = diagdir->mkdir((std::string("Diagnostics_") + Name).c_str());

  std::vector<std::unique_ptr<TH1>> nom_histograms =
      CloneHistVector(NominalHistogramSet, "_nominal_pred");

  FillHistFromstdvector(nom_histograms, NominalPrediction);

  for (std::unique_ptr<TH1> &h : nom_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  TH1D *nominal_1D = new TH1D(
      "nominal_pred_1D", ";Bin Number;#Phi_{#nu} cm^{-2} GeV^{-1} per POT",
      NominalPrediction.size(), 0, NominalPrediction.size());
  FillHistFromstdvector(nominal_1D, NominalPrediction);
  nominal_1D->SetDirectory(td);

  std::vector<std::unique_ptr<TH1>> mean_pred_histograms =
      CloneHistVector(NominalHistogramSet, "_mean_pred");

  std::vector<double> mean_pred_vector;
  std::vector<double> stddev_pred_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    mean_pred_vector.push_back((1 + cb.GetMeanVector()[i]) *
                               NominalPrediction[i]);
    stddev_pred_vector.push_back(cb.GetStdDevVector()[i] *
                                 NominalPrediction[i]);
  }

  FillHistFromstdvector(mean_pred_histograms, mean_pred_vector, 0,
                        stddev_pred_vector);

  for (std::unique_ptr<TH1> &h : mean_pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> stddev_histograms =
      CloneHistVector(NominalHistogramSet, "_stddev");

  TH1D *mean_std_dev = new TH1D("mean_std_dev", "", mean_pred_vector.size(), 0,
                                mean_pred_vector.size());
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
        CloneHistVector(NominalHistogramSet,
                        std::string("_fractional_") + std::to_string(tw_it));

    FillHistFromstdvector(tweak_histograms, RelativeTweaks[tw_it]);

    for (std::unique_ptr<TH1> &h : tweak_histograms) {
      h->SetDirectory(t_id);
      h.release(); // Let root look after the histo again.
    }

    for (size_t bin_it = 0; bin_it < RelativeTweaks[tw_it].size(); ++bin_it) {
      variation_distributions->Fill(bin_it, RelativeTweaks[tw_it][bin_it]);
    }
  }

  if (DumpMatrices) {
    TDirectory *md = td->mkdir("Matrices");
    md->cd();
    std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
    std::unique_ptr<TMatrixD> corrmat =
        GetTMatrixD(CovToCorr(cb.GetCovMatrix()));

    covmat->Write("covmat");
    corrmat->Write("corrmat");
  }
}

void DiscreteVariations::Configure(fhicl::ParameterSet const &ps) {
  BaseConfigure(ps);

  std::vector<std::pair<double, std::vector<double>>> DiscreteTweaks;
  std::vector<std::pair<double, std::vector<double>>> DiscreteErrors;

  for (fhicl::ParameterSet const &tweak_ps :
       ps.get<std::vector<fhicl::ParameterSet>>("DiscreteTweaks")) {

    double SigmaValue = tweak_ps.get<double>("Value");

    std::vector<double> pred_twk, pred_error;

    std::string InputFile_template = tweak_ps.get<std::string>("InputFile");
    std::string VariedHistName_template =
        tweak_ps.get<std::string>("InputHistName");

    for (std::string const &pred_config : Configurations) {

      std::string FluxSliceDescriptor = paramset.get<std::string>(
          std::string("FluxSlicesDescriptor_") + pred_config, "");

      std::string InputFile_tw_c =
          str_replace(InputFile_template, "%C", pred_config);
      std::string VariedHistName_tw_c =
          str_replace(VariedHistName_template, "%C", pred_config);

      for (std::string const &species : Species) {

        std::string InputFile = str_replace(InputFile_tw_c, "%S", species);
        std::string VariedHistName =
            str_replace(VariedHistName_tw_c, "%S", species);

        if (FluxSliceDescriptor.size()) {
          std::vector<std::pair<double, double>> XRanges =
              BuildRangesList(FluxSliceDescriptor);

          std::unique_ptr<TH2D> flux2D =
              GetHistogram_uptr<TH2D>(InputFile, VariedHistName);

          auto const &split = MergeSplitTH2D(flux2D, true, XRanges);
          Mergestdvector(pred_twk, Getstdvector(split));
          Mergestdvector(pred_error, Getstdvectorerror(split));
        } else {
          std::unique_ptr<TH1> var_hist =
              GetHistogram_uptr<TH1>(InputFile, VariedHistName);

          Mergestdvector(pred_twk, Getstdvector(var_hist.get()));
          Mergestdvector(pred_error, Getstdvectorerror(var_hist.get()));
        }
      } // end species
    }   // end configurations

    if (DiscreteTweaks.size() &&
        (DiscreteTweaks.back().second.size() != pred_twk.size())) {
      std::cout << "[ERROR]: Mismatched flux vector length: "
                << RelativeTweaks.back().size() << " != " << pred_twk.size()
                << std::endl;
      throw;
    }

    if (pred_twk.size() != NominalPrediction.size()) {
      std::cout << "[ERROR]: Mismatched nominal flux vector and throw length: "
                << NominalPrediction.size() << " != " << pred_twk.size()
                << std::endl;
      throw;
    }

    if (DumpDiagnostics) {
      diags_Predictions.push_back(pred_twk);
      diags_Errors.push_back(pred_error);
    }

    for (size_t fbin_i = 0; fbin_i < NominalPrediction.size(); ++fbin_i) {
      if (!std::isnormal(pred_twk[fbin_i]) ||
          !std::isnormal(NominalPrediction[fbin_i])) {
        pred_twk[fbin_i] = 0;
        pred_error[fbin_i] = 0;
      } else {
        if (!std::isnormal(pred_twk[fbin_i] / NominalPrediction[fbin_i])) {
          std::cout << "[ERROR]: For flux prediction (" << Name
                    << ") at value: " << SigmaValue << " for " << Name
                    << ", bin " << fbin_i << " was non-normal("
                    << show_classification(pred_twk[fbin_i] /
                                           NominalPrediction[fbin_i])
                    << "): " << pred_twk[fbin_i] << "/"
                    << NominalPrediction[fbin_i] << std::endl;
          throw;
        }
        pred_error[fbin_i] = sqrt(
            pow(pred_error[fbin_i] / NominalPrediction[fbin_i], 2) +
            pow((pred_twk[fbin_i] * NominalPredictionError[fbin_i]) /
                    (NominalPrediction[fbin_i] * NominalPrediction[fbin_i]),
                2));

        pred_twk[fbin_i] = (pred_twk[fbin_i] / NominalPrediction[fbin_i]) - 1;
      }
    }
    DiscreteTweaks.emplace_back(SigmaValue, std::move(pred_twk));
    DiscreteErrors.emplace_back(SigmaValue, std::move(pred_error));
    sig_vals.emplace_back(SigmaValue);
  }

  size_t NSigs = DiscreteTweaks.size();

  std::vector<std::tuple<double, double, double>> xyevals;

  for (size_t fbin_i = 0; fbin_i < NominalPrediction.size(); ++fbin_i) {
    xyevals.clear();

    for (size_t twk_i = 0; twk_i < NSigs; ++twk_i) {
      xyevals.emplace_back(DiscreteTweaks[twk_i].first,
                           DiscreteTweaks[twk_i].second[fbin_i],
                           DiscreteErrors[twk_i].second[fbin_i]);
    }

    std::sort(xyevals.begin(), xyevals.end());

    // Add the nominal == 0 point
    if (std::get<0>(xyevals.back()) < 0) {
      xyevals.emplace_back(
          0, 0, NominalPredictionError[fbin_i] / NominalPrediction[fbin_i]);
    } else {
      for (size_t twk_i = 0; twk_i < NSigs; ++twk_i) {
        if (std::get<0>(xyevals[twk_i]) > 0) {
          xyevals.insert(
              xyevals.begin() + twk_i,
              {0, 0,
               NominalPredictionError[fbin_i] / NominalPrediction[fbin_i]});
          break;
        }
      }
    }

    InterpolatedResponses.push_back(PolyResponse<5>(xyevals));
  }
}

void DiscreteVariations::Process() {

  std::random_device r;

  NDOF = 1;

  std::mt19937 RNEngine(r());
  std::normal_distribution<double> RNJesus(0, 1);

  size_t NBins = InterpolatedResponses.size();

  size_t NThrows = paramset.get<size_t>("NThrows");
  for (size_t t_it = 0; t_it < NThrows; ++t_it) {
    std::vector<double> thrw_pred(NBins);
    double val = RNJesus(RNEngine);

#pragma omp parallel for
    for (size_t fbin_i = 0; fbin_i < NBins; ++fbin_i) {
      thrw_pred[fbin_i] = InterpolatedResponses[fbin_i].eval(val);
    }

    RelativeTweaks.emplace_back(std::move(thrw_pred));
  }

  CovarianceBuilder cb(RelativeTweaks);

  CovarianceComponent = cb.GetCovMatrix();

  if (!(diagdir && DumpDiagnostics)) {
    return;
  }

  TDirectory *td = diagdir->mkdir((std::string("Diagnostics_") + Name).c_str());

  TDirectory *td_int = td->mkdir("Interpolations");

  td_int->cd();
  for (size_t bin_it = 0; bin_it < InterpolatedResponses.size(); ++bin_it) {
    TGraph g_resp(50);
    double min = *std::min_element(sig_vals.begin(), sig_vals.end());
    double max = *std::max_element(sig_vals.begin(), sig_vals.end());
    double step = (max - min) / double(50);
    for (size_t p_it = 0; p_it < 50; ++p_it) {
      double val = min + step * p_it;
      g_resp.SetPoint(p_it, val, InterpolatedResponses[bin_it].eval(val));
    }

    std::stringstream ss("");
    ss << "resp_" << bin_it << std::endl;
    g_resp.Write(ss.str().c_str());

    TGraphErrors g_pred(sig_vals.size() + 1);
    std::vector<std::tuple<double, double, double>> xyevals;
    for (size_t v_it = 0; v_it < sig_vals.size(); ++v_it) {

      std::cout << "[" << bin_it << "] , var = " << v_it
                << ",  bc = " << diags_Predictions[v_it][bin_it]
                << ", be = " << diags_Errors[v_it][bin_it]
                << ", NP = " << NominalPrediction[bin_it]
                << ", NE = " << NominalPredictionError[bin_it] << std::endl;

      double err =
          sqrt(pow(diags_Errors[v_it][bin_it] / NominalPrediction[bin_it], 2) +
               pow((diags_Predictions[v_it][bin_it] *
                    NominalPredictionError[bin_it]) /
                       (NominalPrediction[bin_it] * NominalPrediction[bin_it]),
                   2));

      xyevals.emplace_back(
          sig_vals[v_it],
          (diags_Predictions[v_it][bin_it] / NominalPrediction[bin_it]) - 1,
          err);
    }

    double Nom_err =
        (NominalPredictionError[bin_it] / NominalPrediction[bin_it]);

    xyevals.emplace_back(0, 0, Nom_err);
    std::sort(xyevals.begin(), xyevals.end());
    for (size_t v_it = 0; v_it < xyevals.size(); ++v_it) {
      g_pred.SetPoint(v_it, std::get<0>(xyevals[v_it]),
                      std::get<1>(xyevals[v_it]));
      g_pred.SetPointError(v_it, 0, std::get<2>(xyevals[v_it]));

      std::cout << "[" << bin_it << "] , var = " << v_it << ", "
                << std::get<0>(xyevals[v_it]) << ", "
                << std::get<1>(xyevals[v_it]) << ", "
                << std::get<2>(xyevals[v_it]) << std::endl;
    }
    ss.str("");
    ss << "binc_" << bin_it << std::endl;
    g_pred.Write(ss.str().c_str());
  }

  std::vector<std::unique_ptr<TH1>> nom_histograms =
      CloneHistVector(NominalHistogramSet, "_nominal_pred");

  FillHistFromstdvector(nom_histograms, NominalPrediction);

  for (std::unique_ptr<TH1> &h : nom_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  TH1D *nominal_1D = new TH1D(
      "nominal_pred_1D", ";Bin Number;#Phi_{#nu} cm^{-2} GeV^{-1} per POT",
      NominalPrediction.size(), 0, NominalPrediction.size());
  FillHistFromstdvector(nominal_1D, NominalPrediction);
  nominal_1D->SetDirectory(td);

  std::vector<std::unique_ptr<TH1>> mean_pred_histograms =
      CloneHistVector(NominalHistogramSet, "_mean_pred");

  std::vector<double> mean_pred_vector;
  std::vector<double> stddev_pred_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    mean_pred_vector.push_back((1 + cb.GetMeanVector()[i]) *
                               NominalPrediction[i]);
    stddev_pred_vector.push_back(cb.GetStdDevVector()[i] *
                                 NominalPrediction[i]);
  }

  FillHistFromstdvector(mean_pred_histograms, mean_pred_vector, 0,
                        stddev_pred_vector);

  for (std::unique_ptr<TH1> &h : mean_pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> stddev_histograms =
      CloneHistVector(NominalHistogramSet, "_stddev");

  TH1D *mean_std_dev = new TH1D("mean_std_dev", "", mean_pred_vector.size(), 0,
                                mean_pred_vector.size());
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

  for (size_t tw_it = 0; tw_it < sig_vals.size(); ++tw_it) {

    std::string sig_val_str = std::string(sig_vals[tw_it] < 0 ? "m" : "") +
                              std::to_string(fabs(sig_vals[tw_it]));

    TDirectory *t_id =
        td->mkdir((std::string("_tweak_") + sig_val_str).c_str());

    std::vector<std::unique_ptr<TH1>> pred_histograms = CloneHistVector(
        NominalHistogramSet, std::string("_pred_") + sig_val_str);

    FillHistFromstdvector(pred_histograms, diags_Predictions[tw_it]);

    for (std::unique_ptr<TH1> &h : pred_histograms) {
      h->SetDirectory(t_id);
      h.release(); // Let root look after the histo again.
    }

    std::vector<std::unique_ptr<TH1>> tweak_histograms = CloneHistVector(
        NominalHistogramSet, std::string("_fractional_") + sig_val_str);

    std::vector<double> frac(diags_Predictions[tw_it].size());
    for (size_t bin_it = 0; bin_it < diags_Predictions[tw_it].size();
         ++bin_it) {
      frac[bin_it] =
          ((diags_Predictions[tw_it][bin_it] / NominalPrediction[bin_it]) - 1);
    }

    FillHistFromstdvector(tweak_histograms, frac);

    for (std::unique_ptr<TH1> &h : tweak_histograms) {
      h->SetDirectory(t_id);
      h.release(); // Let root look after the histo again.
    }

    for (size_t bin_it = 0; bin_it < diags_Predictions[tw_it].size();
         ++bin_it) {
      variation_distributions->Fill(bin_it, diags_Predictions[tw_it][bin_it]);
    }
  }

  if (DumpMatrices) {
    TDirectory *md = td->mkdir("Matrices");
    md->cd();
    std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
    std::unique_ptr<TMatrixD> corrmat =
        GetTMatrixD(CovToCorr(cb.GetCovMatrix()));

    covmat->Write("covmat");
    corrmat->Write("corrmat");
  }
}

void DirectVariations::Configure(fhicl::ParameterSet const &ps) {
  BaseConfigure(ps);

  NDOF = 1;

  std::vector<double> pred_twk;

  std::string InputFile_template =
      paramset.get<std::string>("Varied.InputFile");
  std::string VariedHistName_template =
      paramset.get<std::string>("Varied.InputHistName");

  for (std::string const &pred_config : Configurations) {

    std::string FluxSliceDescriptor = paramset.get<std::string>(
        std::string("FluxSlicesDescriptor_") + pred_config, "");

    std::string InputFile_tw_c =
        str_replace(InputFile_template, "%C", pred_config);

    std::string VariedHistName_tw_c =
        str_replace(VariedHistName_template, "%C", pred_config);

    for (std::string const &species : Species) {

      std::string InputFile = str_replace(InputFile_tw_c, "%S", species);
      std::string VariedHistName =
          str_replace(VariedHistName_tw_c, "%S", species);

      if (FluxSliceDescriptor.size()) {
        std::vector<std::pair<double, double>> XRanges =
            BuildRangesList(FluxSliceDescriptor);

        std::unique_ptr<TH2D> flux2D =
            GetHistogram_uptr<TH2D>(InputFile, VariedHistName);

        Mergestdvector(pred_twk,
                       Getstdvector(MergeSplitTH2D(flux2D, true, XRanges)));
      } else {
        std::unique_ptr<TH1> var_hist =
            GetHistogram_uptr<TH1>(InputFile, VariedHistName);

        Mergestdvector(pred_twk, Getstdvector(var_hist.get()));
      }
    } // end species
  }   // end configurations

  if (pred_twk.size() != NominalPrediction.size()) {
    std::cout << "[ERROR]: Mismatched nominal flux vector and throw length: "
              << NominalPrediction.size() << " != " << pred_twk.size()
              << std::endl;
    throw;
  }

  if (DumpDiagnostics) {
    diags_Predictions.push_back(pred_twk);
  }

  for (size_t fbin_i = 0; fbin_i < NominalPrediction.size(); ++fbin_i) {
    if (!std::isnormal(pred_twk[fbin_i]) ||
        !std::isnormal(NominalPrediction[fbin_i])) {
      pred_twk[fbin_i] = 0;
    } else {
      if (!std::isnormal(pred_twk[fbin_i] / NominalPrediction[fbin_i])) {
        std::cout << "[ERROR]: For flux prediction: " << Name << ", bin "
                  << fbin_i << " was non-normal("
                  << show_classification(pred_twk[fbin_i] /
                                         NominalPrediction[fbin_i])
                  << "): " << pred_twk[fbin_i] << "/"
                  << NominalPrediction[fbin_i] << std::endl;
        throw;
      }
      pred_twk[fbin_i] = (pred_twk[fbin_i] / NominalPrediction[fbin_i]) - 1;
    }
  }

  RelativeTweaks.push_back(std::move(pred_twk));
}

void DirectVariations::Process() {
  CovarianceBuilder cb(RelativeTweaks);

  CovarianceComponent = cb.GetCovMatrix();

  if (!(diagdir && DumpDiagnostics)) {
    return;
  }

  TDirectory *td = diagdir->mkdir((std::string("Diagnostics_") + Name).c_str());

  std::vector<std::unique_ptr<TH1>> nom_histograms =
      CloneHistVector(NominalHistogramSet, "_nominal_pred");

  FillHistFromstdvector(nom_histograms, NominalPrediction);

  for (std::unique_ptr<TH1> &h : nom_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  TH1D *nominal_1D = new TH1D(
      "nominal_pred_1D", ";Bin Number;#Phi_{#nu} cm^{-2} GeV^{-1} per POT",
      NominalPrediction.size(), 0, NominalPrediction.size());
  FillHistFromstdvector(nominal_1D, NominalPrediction);
  nominal_1D->SetDirectory(td);

  std::vector<std::unique_ptr<TH1>> mean_pred_histograms =
      CloneHistVector(NominalHistogramSet, "_mean_pred");

  std::vector<double> mean_pred_vector;
  std::vector<double> stddev_pred_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    mean_pred_vector.push_back((1 + cb.GetMeanVector()[i]) *
                               NominalPrediction[i]);
    stddev_pred_vector.push_back(cb.GetStdDevVector()[i] *
                                 NominalPrediction[i]);
  }

  FillHistFromstdvector(mean_pred_histograms, mean_pred_vector, 0,
                        stddev_pred_vector);

  for (std::unique_ptr<TH1> &h : mean_pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> stddev_histograms =
      CloneHistVector(NominalHistogramSet, "_stddev");

  std::vector<double> stddev_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    stddev_vector.push_back(cb.GetStdDevVector()[i]);
  }

  FillHistFromstdvector(stddev_histograms, stddev_vector);

  for (std::unique_ptr<TH1> &h : stddev_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> pred_histograms =
      CloneHistVector(NominalHistogramSet, "_pred");

  FillHistFromstdvector(pred_histograms, diags_Predictions.front());

  for (std::unique_ptr<TH1> &h : pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> pred_fractional_histograms =
      CloneHistVector(NominalHistogramSet, "_pred_fractional");

  std::vector<double> frac(diags_Predictions.front().size());
  for (size_t bin_it = 0; bin_it < diags_Predictions.front().size(); ++bin_it) {
    frac[bin_it] =
        ((diags_Predictions.front()[bin_it] / NominalPrediction[bin_it]) - 1);
  }

  FillHistFromstdvector(pred_fractional_histograms, frac);

  for (std::unique_ptr<TH1> &h : pred_fractional_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  if (DumpMatrices) {
    TDirectory *md = td->mkdir("Matrices");
    md->cd();
    std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
    std::unique_ptr<TMatrixD> corrmat =
        GetTMatrixD(CovToCorr(cb.GetCovMatrix()));

    covmat->Write("covmat");
    corrmat->Write("corrmat");
  }
}

void UniformVariations::Configure(fhicl::ParameterSet const &ps) {
  BaseConfigure(ps);

  NDOF = 1;

  double Uncertainty_pc = ps.get<double>("Uncertainty_pc");

  std::vector<double> pred_twk;

  for (double &np : NominalPrediction) {
    pred_twk.push_back(np * (1 + Uncertainty_pc / 100.0));
  }

  if (DumpDiagnostics) {
    diags_Predictions.push_back(pred_twk);
  }

  for (size_t fbin_i = 0; fbin_i < NominalPrediction.size(); ++fbin_i) {
    if (!std::isnormal(pred_twk[fbin_i]) ||
        !std::isnormal(NominalPrediction[fbin_i])) {
      pred_twk[fbin_i] = 0;
    } else {
      if (!std::isnormal(pred_twk[fbin_i] / NominalPrediction[fbin_i])) {
        std::cout << "[ERROR]: For flux prediction: " << Name << ", bin "
                  << fbin_i << " was non-normal("
                  << show_classification(pred_twk[fbin_i] /
                                         NominalPrediction[fbin_i])
                  << "): " << pred_twk[fbin_i] << "/"
                  << NominalPrediction[fbin_i] << std::endl;
        throw;
      }
      pred_twk[fbin_i] = (pred_twk[fbin_i] / NominalPrediction[fbin_i]) - 1;
    }
  }

  RelativeTweaks.push_back(std::move(pred_twk));
}

void UniformVariations::Process() {
  CovarianceBuilder cb(RelativeTweaks);

  CovarianceComponent = cb.GetCovMatrix();

  if (!(diagdir && DumpDiagnostics)) {
    return;
  }

  TDirectory *td = diagdir->mkdir((std::string("Diagnostics_") + Name).c_str());

  std::vector<std::unique_ptr<TH1>> nom_histograms =
      CloneHistVector(NominalHistogramSet, "_nominal_pred");

  FillHistFromstdvector(nom_histograms, NominalPrediction);

  for (std::unique_ptr<TH1> &h : nom_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  TH1D *nominal_1D = new TH1D(
      "nominal_pred_1D", ";Bin Number;#Phi_{#nu} cm^{-2} GeV^{-1} per POT",
      NominalPrediction.size(), 0, NominalPrediction.size());
  FillHistFromstdvector(nominal_1D, NominalPrediction);
  nominal_1D->SetDirectory(td);

  std::vector<std::unique_ptr<TH1>> mean_pred_histograms =
      CloneHistVector(NominalHistogramSet, "_mean_pred");

  std::vector<double> mean_pred_vector;
  std::vector<double> stddev_pred_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    mean_pred_vector.push_back((1 + cb.GetMeanVector()[i]) *
                               NominalPrediction[i]);
    stddev_pred_vector.push_back(cb.GetStdDevVector()[i] *
                                 NominalPrediction[i]);
  }

  FillHistFromstdvector(mean_pred_histograms, mean_pred_vector, 0,
                        stddev_pred_vector);

  for (std::unique_ptr<TH1> &h : mean_pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> stddev_histograms =
      CloneHistVector(NominalHistogramSet, "_stddev");

  std::vector<double> stddev_vector;
  for (int i = 0; i < cb.NRows; ++i) {
    stddev_vector.push_back(cb.GetStdDevVector()[i]);
  }

  FillHistFromstdvector(stddev_histograms, stddev_vector);

  for (std::unique_ptr<TH1> &h : stddev_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> pred_histograms =
      CloneHistVector(NominalHistogramSet, "_pred");

  FillHistFromstdvector(pred_histograms, diags_Predictions.front());

  for (std::unique_ptr<TH1> &h : pred_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  std::vector<std::unique_ptr<TH1>> pred_fractional_histograms =
      CloneHistVector(NominalHistogramSet, "_pred_fractional");

  std::vector<double> frac(diags_Predictions.front().size());
  for (size_t bin_it = 0; bin_it < diags_Predictions.front().size(); ++bin_it) {
    frac[bin_it] =
        ((diags_Predictions.front()[bin_it] / NominalPrediction[bin_it]) - 1);
  }

  FillHistFromstdvector(pred_fractional_histograms, frac);

  for (std::unique_ptr<TH1> &h : pred_fractional_histograms) {
    h->SetDirectory(td);
    h.release(); // Let root look after the histo again.
  }

  if (DumpMatrices) {
    TDirectory *md = td->mkdir("Matrices");
    md->cd();
    std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
    std::unique_ptr<TMatrixD> corrmat =
        GetTMatrixD(CovToCorr(cb.GetCovMatrix()));

    covmat->Write("covmat");
    corrmat->Write("corrmat");
  }
}
