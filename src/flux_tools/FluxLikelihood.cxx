#include "FluxLikelihood.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

// #define DEBUG_LHOOD

void FluxLinearCombiner::Initialize(FluxSliceOptions const &fso) {
  fFSO = fso;
  LoadFluxes();

  for (size_t i = 0; i < fFSO.XRanges.size(); ++i) {
    if (fFSO.RegFactor != 0xdeadbeef) {
      if (i != 0) {
        // Don't regularise across gaps
        if (fabs(fFSO.XRanges[i].first - fFSO.XRanges[i - 1].second) > 1E-5) {
          ApplyRegularization.back() = false;
        }
      }
      ApplyRegularization.push_back(true);
    }
  }
}

void FluxLinearCombiner::LoadFluxes() {
  std::unique_ptr<TH2D> Flux2D =
      GetHistogram_uptr<TH2D>(fFSO.InputFluxFile, fFSO.InputFluxName);

  if (!Flux2D) {
    std::cout << "[ERROR]: Found no input flux with name: \""
              << fFSO.InputFluxName << "\" in file: \"" << fFSO.InputFluxFile
              << "\"." << std::endl;
    throw;
  }

  if (fFSO.MergeENuBins && fFSO.MergeOAPBins) {
    Flux2D->Rebin2D(fFSO.MergeENuBins, fFSO.MergeOAPBins);
    Flux2D->Scale(1.0 / double(fFSO.MergeENuBins));
  } else if (fFSO.MergeENuBins) {
    Flux2D->RebinX(fFSO.MergeENuBins);
    Flux2D->Scale(1.0 / double(fFSO.MergeENuBins));
  } else if (fFSO.MergeOAPBins) {
    Flux2D->RebinY(fFSO.MergeOAPBins);
    Flux2D->Scale(1.0 / double(fFSO.MergeOAPBins));
  }

  if (fFSO.XRanges.size()) {
    FluxSlices = MergeSplitTH2D(Flux2D, true, fFSO.XRanges);
  } else {
    std::vector<std::pair<std::pair<double, double>, TH1D *>> SplitFlux2D =
        SplitTH2D(Flux2D.get(), true);
    for (size_t i = 0; i < SplitFlux2D.size(); ++i) {
      fFSO.XRanges.push_back(SplitFlux2D[i].first);
      FluxSlices.emplace_back(SplitFlux2D[i].second);
    }
  }

  if (!FluxSlices.size()) {
    std::cout << "[ERROR]: Found no input fluxes." << std::endl;
    throw;
  }

  SummedFlux = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(FluxSlices.front()->Clone("SummedFlux")));
  SummedFlux->SetDirectory(nullptr);

  SummedFlux->Reset();
}

size_t FluxLinearCombiner::GetNearestPeakingFluxIndex(double enu) const {
  size_t flux_with_closest_peak_it = -1;
  double flux_peak_dist = std::numeric_limits<double>::max();
  for (size_t flux_it = 0; flux_it < FluxSlices.size(); flux_it++) {
    double maxenu = FluxSlices[flux_it]->GetXaxis()->GetBinCenter(
        FluxSlices[flux_it]->GetMaximumBin());
    double peak_difference = fabs(maxenu - enu);
    if (peak_difference < flux_peak_dist) {
      flux_peak_dist = peak_difference;
      flux_with_closest_peak_it = flux_it;
    }
  }
  return flux_with_closest_peak_it;
}

void FluxLinearCombiner::SumFluxes(double const *FluxCoefficients) const {
  SumHistograms(SummedFlux.get(), FluxCoefficients, FluxSlices);
}

double FluxLinearCombiner::GetRegularizationPenalty(
    double const *FluxCoefficients) const {
  if (fFSO.RegFactor == 0xdeadbeef) {
    return 0;
  }

  double reg = 0;

  for (size_t i = 0; i < FluxSlices.size(); i++) {
    if (!ApplyRegularization[i]) {
      continue;
    }
    reg += pow((FluxCoefficients[i] - FluxCoefficients[i + 1]) / fFSO.RegFactor,
               2);
  }

  return reg;
}

void FluxTargetLikelihood::Initialize(
    FluxTargetLikelihood::FluxTargetOptions const &fto) {
  fLHoodOpts = fto;

  FluxLinearCombiner::Initialize(fLHoodOpts.input_flux_opts);

  TargetFlux = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(FluxSlices.front()->Clone("TargetFlux")));
  TargetFlux->SetDirectory(nullptr);
  TargetFlux->Reset();
}

size_t FluxTargetLikelihood::GetNearestPeakingFluxIndex() const {
  return GetNearestPeakingFluxIndex(
      FluxSlices[0]->GetXaxis()->GetBinCenter(TargetFlux->GetMaximumBin()));
}

void FluxTargetLikelihood::SetTargetFlux(TH1D const *tgt) {
  InputFlux =
      std::unique_ptr<TH1D>(static_cast<TH1D *>(tgt->Clone("InputFlux")));
  InputFlux->SetDirectory(nullptr);

  TargetFlux->Reset();

  for (Int_t bi_it = 1; bi_it < TargetFlux->GetXaxis()->GetNbins() + 1;
       ++bi_it) {
    if (!std::isnormal(tgt->GetBinContent(bi_it))) {
      std::cout << "[ERROR]: Bad bin content found in new TargetFlux, bin "
                << bi_it << std::endl;
      throw;
    }
    TargetFlux->SetBinContent(bi_it, tgt->GetBinContent(bi_it));
  }

  if (fLHoodOpts.FitBetweenFoundPeaks) {
    FindTH1Peaks(TargetFlux.get(), fFitBinLow, fFitBinHigh, 3);
    if (fFitBinLow == 0) {
      std::cout << "[ERROR]: Failed to find the expected number of "
                   "peaks, "
                << std::endl;
      throw;
    }
    fLHoodOpts.FitBetween.first =
        TargetFlux->GetXaxis()->GetBinLowEdge(fFitBinLow);
    fLHoodOpts.FitBetween.second =
        TargetFlux->GetXaxis()->GetBinUpEdge(fFitBinHigh);
  } else {
    if (fLHoodOpts.FitBetween.first == 0xdeadbeef) {
      fFitBinLow = 1;
      fFitBinHigh = TargetFlux->GetXaxis()->GetNbins();
    } else {
      fFitBinLow =
          TargetFlux->GetXaxis()->FindFixBin(fLHoodOpts.FitBetween.first);
      fFitBinHigh =
          TargetFlux->GetXaxis()->FindFixBin(fLHoodOpts.FitBetween.second);
    }
  }

#ifdef DEBUG_LHOOD

  std::cout << "fFitBinLow: " << fFitBinLow << std::endl;
  std::cout << "fFitBinHigh: " << fFitBinHigh << std::endl;
  std::cout << "FitBetween.first: " << fLHoodOpts.FitBetween.first << std::endl;
  std::cout << "FitBetween.second: " << fLHoodOpts.FitBetween.second
            << std::endl;
  std::cout << "OutOfRangeMode: " << fLHoodOpts.OutOfRangeMode << std::endl;
#endif

  // Build target to the left of the fit region
  for (Int_t bi_it = 1; bi_it < fFitBinLow; ++bi_it) {
    if ((fLHoodOpts.OutOfRangeMode != FluxTargetOptions::kIgnore) &&
        (fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kBoth ||
         fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kLeft)) {
      double target = 0;
      double enu_first_counted_bin =
          TargetFlux->GetXaxis()->GetBinCenter(fFitBinLow);
      double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
      double content_first_counted_bin = TargetFlux->GetBinContent(fFitBinLow);
      double enu_bottom_bin = TargetFlux->GetXaxis()->GetBinCenter(1);
      double sigma5_range = enu_first_counted_bin - enu_bottom_bin;

      if (fLHoodOpts.OutOfRangeMode == FluxTargetOptions::kExponentialDecay) {
        target = content_first_counted_bin *
                 exp(-fLHoodOpts.ExpDecayRate * (enu_first_counted_bin - enu) /
                     sigma5_range);
      } else if (fLHoodOpts.OutOfRangeMode ==
                 FluxTargetOptions::kGaussianDecay) {
        target =
            content_first_counted_bin *
            exp(-fLHoodOpts.ExpDecayRate * (enu_first_counted_bin - enu) *
                (enu_first_counted_bin - enu) / (sigma5_range * sigma5_range));
      }

      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, fLHoodOpts.TargetFractionalFreedom *
                                       TargetFlux->GetBinContent(bi_it));
  }

  // Set the target errors within the fit region
  for (Int_t bi_it = fFitBinLow; bi_it < (fFitBinHigh + 1); ++bi_it) {
    TargetFlux->SetBinError(bi_it, fLHoodOpts.TargetFractionalFreedom *
                                       TargetFlux->GetBinContent(bi_it));
  }

  // Build the target above the fit region
  for (Int_t bi_it = (fFitBinHigh + 1);
       bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
    if ((fLHoodOpts.OutOfRangeMode != FluxTargetOptions::kIgnore) &&
        (fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kBoth ||
         fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kRight)) {
      double target = 0;
      double enu_last_counted_bin =
          TargetFlux->GetXaxis()->GetBinCenter(fFitBinHigh);
      double enu = TargetFlux->GetXaxis()->GetBinCenter(bi_it);
      double content_last_counted_bin = TargetFlux->GetBinContent(fFitBinHigh);
      double enu_top_bin = TargetFlux->GetXaxis()->GetBinCenter(
          TargetFlux->GetXaxis()->GetNbins());
      double sigma5_range = enu_top_bin - enu_last_counted_bin;

      if (fLHoodOpts.OutOfRangeMode == FluxTargetOptions::kExponentialDecay) {
        target = content_last_counted_bin *
                 exp(-fLHoodOpts.ExpDecayRate * (enu - enu_last_counted_bin) /
                     sigma5_range);
      } else if (fLHoodOpts.OutOfRangeMode ==
                 FluxTargetOptions::kGaussianDecay) {
        target =
            content_last_counted_bin *
            exp(-fLHoodOpts.ExpDecayRate * (enu - enu_last_counted_bin) *
                (enu - enu_last_counted_bin) / (sigma5_range * sigma5_range));
      }
      TargetFlux->SetBinContent(bi_it, target);
    } else {
      TargetFlux->SetBinContent(bi_it, 0);
    }

    TargetFlux->SetBinError(bi_it, fLHoodOpts.TargetFractionalFreedom *
                                       TargetFlux->GetBinContent(bi_it));
  }

  double OnAxisPeak = FluxSlices[0]->GetMaximum();
  double TargetMax = TargetFlux->GetMaximum();
  fNDOverFDFitScaleFactor = OnAxisPeak / TargetMax;

  if (!std::isnormal(fNDOverFDFitScaleFactor)) {
    std::cout << "[ERROR]: Peak scaling failed: OnAxisPeak = "
              << FluxSlices[0]->GetMaximum()
              << ", TargetMax = " << TargetFlux->GetMaximum() << std::endl;
    throw;
  }

  TargetFlux->Scale(fNDOverFDFitScaleFactor);

  fTargetPeakNorm = FluxSlices[GetNearestPeakingFluxIndex()]->GetMaximum();

#ifdef DEBUG_LHOOD
  for (Int_t bi_it = 1; bi_it < TargetFlux->GetXaxis()->GetNbins() + 1;
       ++bi_it) {
    std::cout << "[TARGET] Bin[" << bi_it
              << "] Enu = " << TargetFlux->GetXaxis()->GetBinCenter(bi_it)
              << ", content = " << TargetFlux->GetBinContent(bi_it)
              << std::endl;
  }
#endif
}

double FluxTargetLikelihood::BinDiff(Int_t bin_it) const {

  double diff, err;
  if (fLHoodOpts.UseNuPrismChi2) {
    err = SummedFlux->GetBinContent(bin_it)
              ? (0.0001 * SummedFlux->GetBinContent(bin_it) *
                 SummedFlux->GetBinContent(bin_it))
              : (0.0001 * fTargetPeakNorm * fTargetPeakNorm);
    diff = ((TargetFlux->GetBinContent(bin_it) -
             SummedFlux->GetBinContent(bin_it)) *
            (TargetFlux->GetBinContent(bin_it) -
             SummedFlux->GetBinContent(bin_it)) /
            err);
  } else {
    err = (TargetFlux->GetBinError(bin_it) * TargetFlux->GetBinError(bin_it) +
           SummedFlux->GetBinError(bin_it) * SummedFlux->GetBinError(bin_it));

    diff = ((TargetFlux->GetBinContent(bin_it) -
             SummedFlux->GetBinContent(bin_it)) *
            (TargetFlux->GetBinContent(bin_it) -
             SummedFlux->GetBinContent(bin_it)) /
            err);
  }

  if (!err) {
    return 0;
  }

#ifdef DEBUG_LHOOD
  if (diff && !std::isnormal(diff)) {
    std::cout << "[ERROR]: Found invalid diff, bin " << bin_it
              << " [ORR LowE], targ: " << TargetFlux->GetBinContent(bin_it)
              << ":" << TargetFlux->GetBinError(bin_it)
              << ", sum: " << TargetFlux->GetBinContent(bin_it) << ":"
              << TargetFlux->GetBinError(bin_it) << ", err = " << err
              << ", diff = " << diff << std::endl;
    throw;
  }
#endif

  return diff;
}

FluxTargetLikelihood::LHood FluxTargetLikelihood::GetLikelihoodComponents(
    double const *coefficients) const {
  SumFluxes(coefficients);

  LHood lh;

  lh.InFitRegionChi2 = 0;
  lh.OutOfFitRegionLH = 0;

  if ((fLHoodOpts.OutOfRangeMode != FluxTargetOptions::kIgnore)) {

    if ((fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kBoth ||
         fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kLeft)) {
      for (Int_t bi_it = 1; bi_it < fFitBinLow; ++bi_it) {
        double diff = fLHoodOpts.OORFactor * BinDiff(bi_it);
        lh.OutOfFitRegionLH += diff;
      }
    }

    if ((fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kBoth ||
         fLHoodOpts.OutOfRangeSide == FluxTargetOptions::kRight)) {
      for (Int_t bi_it = (fFitBinHigh + 1);
           bi_it < TargetFlux->GetXaxis()->GetNbins() + 1; ++bi_it) {
        double diff = fLHoodOpts.OORFactor * BinDiff(bi_it);
        lh.OutOfFitRegionLH += diff;
      }
    }
  }

  for (Int_t bi_it = fFitBinLow; bi_it < (fFitBinHigh + 1); ++bi_it) {

    double sumdiff = BinDiff(bi_it);
    lh.InFitRegionChi2 += sumdiff;

#ifdef DEBUG_LHOOD
    if (sumdiff && !std::isnormal(sumdiff)) {
      std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                << ", targ: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it)
                << ", sum: " << TargetFlux->GetBinContent(bi_it) << ":"
                << TargetFlux->GetBinError(bi_it) << std::endl;
      throw;
    }
#endif
  }

  lh.RegLH = GetRegularizationPenalty(coefficients);

  return lh;
}

void FluxTargetLikelihood::Write(TDirectory *oupD) {
  TH1D *sf = static_cast<TH1D *>(SummedFlux->Clone("BestFit"));
  sf->SetDirectory(oupD);

  TH1D *tf = static_cast<TH1D *>(TargetFlux->Clone("Target"));
  tf->SetDirectory(oupD);

  TH1D *inpf = static_cast<TH1D *>(InputFlux->Clone("InputFlux"));
  inpf->SetDirectory(oupD);

  TH1D *sf_on = static_cast<TH1D *>(SummedFlux->Clone("BestFit_OriginalNorm"));
  sf_on->Scale(1.0 / fNDOverFDFitScaleFactor);
  sf_on->SetDirectory(oupD);

  TH1D *tf_on = static_cast<TH1D *>(TargetFlux->Clone("Target_OriginalNorm"));
  tf_on->Scale(1.0 / fNDOverFDFitScaleFactor);
  tf_on->SetDirectory(oupD);

  TH1D *if_tn = static_cast<TH1D *>(InputFlux->Clone("InputFlux_TargetNorm"));
  if_tn->Scale(fNDOverFDFitScaleFactor);
  if_tn->SetDirectory(oupD);
}

void GausTargetLikelihood::Initialize(
    GausTargetLikelihood::GausTargetOptions const &gto) {
  fLHoodOpts = gto;

  FluxLinearCombiner::Initialize(fLHoodOpts.input_flux_opts);

  TargetGauss = std::make_unique<TF1>("tFunc", "gaus", 0, 20);

  SetGaussParameters(fLHoodOpts.GaussC, fLHoodOpts.GaussW);
}

size_t GausTargetLikelihood::GetNearestPeakingFluxIndex() const {
  return GetNearestPeakingFluxIndex(
      FluxSlices[0]->GetXaxis()->GetBinCenter(fLHoodOpts.GaussC));
}

void GausTargetLikelihood::SetGaussParameters(double GaussC, double GaussW) {
  fLHoodOpts.GaussC = GaussC;
  fLHoodOpts.GaussW = GaussW;

  if ((fLHoodOpts.GaussC <= 0) || (fLHoodOpts.GaussW <= 0)) {
    std::cout << "[ERROR]: Pass invalid argument for -g, expected "
                 "\"<GaussMean>,<GaussWidth>\", where both values are > 0."
              << std::endl;
    throw;
  }

  fTargetPeakNorm = FluxSlices[GetNearestPeakingFluxIndex()]->GetMaximum();

  TargetGauss->SetParameter(0, fTargetPeakNorm);
  TargetGauss->SetParameter(1, fLHoodOpts.GaussC);
  TargetGauss->SetParameter(2, fLHoodOpts.GaussW);

  if (fLHoodOpts.FitBetween.first == 0xdeadbeef) {
    fFitBinLow = 1;
    fFitBinHigh = SummedFlux->GetXaxis()->GetNbins();
  } else {
    fFitBinLow =
        SummedFlux->GetXaxis()->FindFixBin(fLHoodOpts.FitBetween.first);
    fFitBinHigh =
        SummedFlux->GetXaxis()->FindFixBin(fLHoodOpts.FitBetween.second);
  }
}

double GausTargetLikelihood::GetLikelihood(double const *coefficients) const {

  SumFluxes(coefficients);

  double sumdiff = 0;

  for (Int_t bi_it = fFitBinLow; bi_it < fFitBinHigh + 1; ++bi_it) {
    double bi_c_E = SummedFlux->GetXaxis()->GetBinCenter(bi_it);

    double GaussEval = TargetGauss->Eval(bi_c_E);
    double SummedBinContent = SummedFlux->GetBinContent(bi_it);

    double uncert;
    if (fLHoodOpts.UseNuPrismChi2) {
      uncert = pow(0.0001 * GaussEval, 2) + pow(0.00005 * fTargetPeakNorm, 2);
    } else {
      uncert = pow(GaussEval / 20.0, 2) + pow(fTargetPeakNorm / 30.0, 2);
    }

    sumdiff += (pow(GaussEval - SummedBinContent, 2) / uncert);

#ifdef DEBUG_LHOOD
    if (sumdiff && !std::isnormal(sumdiff)) {
      std::cout << "[ERROR]: Found invalid diff, bin " << bi_it
                << ", gauss: " << GaussEval << ":" << fTargetPeakNorm
                << ", sum: " << SummedFlux->GetBinContent(bi_it) << ":"
                << SummedFlux->GetBinError(bi_it) << std::endl;
      throw;
    }
#endif
  }
  return sumdiff + GetRegularizationPenalty(coefficients);
}

#ifdef USE_FHICL_CPP
FluxLinearCombiner::FluxSliceOptions
MakeFluxSliceOptions(fhicl::ParameterSet const &);
FluxTargetLikelihood::FluxTargetOptions
MakeFluxTargetOptions(fhicl::ParameterSet const &);
GausTargetLikelihood::GausTargetOptions
MakeGausTargetOptions(fhicl::ParameterSet const &);
#endif

FluxLinearCombiner::FluxSliceOptions MakeFluxSliceOptions(int argc,
                                                          char const *argv[]) {
  FluxLinearCombiner::FluxSliceOptions fso;
  fso.InputFluxFile = "";
  fso.InputFluxName = "";
  fso.MergeOAPBins = 0;
  fso.MergeENuBins = 0;
  fso.RegFactor = 0xdeadbeef;

  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-f") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -f, expected 2." << std::endl;
        throw;
      }
      fso.InputFluxFile = params[0];
      fso.InputFluxName = params[1];
    } else if (std::string(argv[opt]) == "-rg") {
      fso.RegFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-MX") {
      fso.MergeENuBins = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-MY") {
      fso.MergeOAPBins = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-M") {
      fso.XRanges = BuildRangesList(argv[++opt]);
    }
    opt++;
  }
  return fso;
}

FluxTargetLikelihood::FluxTargetOptions
MakeFluxTargetOptions(int argc, char const *argv[]) {
  FluxTargetLikelihood::FluxTargetOptions fto;

  fto.input_flux_opts = MakeFluxSliceOptions(argc, argv);

  fto.FitBetweenFoundPeaks = false;
  fto.FitBetween = {0xdeadbeef, 0xdeadbeef};
  fto.OutOfRangeMode = FluxTargetLikelihood::FluxTargetOptions::kZero;
  fto.OutOfRangeSide = FluxTargetLikelihood::FluxTargetOptions::kBoth;
  fto.OORFactor = 1;
  fto.ExpDecayRate = 3;
  fto.TargetFractionalFreedom = 0.01;
  fto.UseNuPrismChi2 = false;

  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for l, expected 2." << std::endl;
        throw;
      }
      fto.FitBetween = std::make_pair(params[0], params[1]);
    } else if (std::string(argv[opt]) == "-p") {
      fto.FitBetweenFoundPeaks = true;
    } else if (std::string(argv[opt]) == "-m") {
      fto.OutOfRangeMode = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ed") {
      fto.ExpDecayRate = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-of") {
      fto.OORFactor = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-ms") {
      fto.OutOfRangeSide = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-C") {
      fto.UseNuPrismChi2 = true;
    }
    opt++;
  }

  return fto;
}

GausTargetLikelihood::GausTargetOptions
MakeGausTargetOptions(int argc, char const *argv[]) {
  GausTargetLikelihood::GausTargetOptions gto;

  gto.GaussC = -1;
  gto.GaussW = -1;
  gto.FitBetween = {0xdeadbeef, 0xdeadbeef};
  gto.UseNuPrismChi2 = false;

  gto.input_flux_opts = MakeFluxSliceOptions(argc, argv);

  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-g") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for g, expected 2." << std::endl;
        throw;
      }
      gto.GaussC = params[0];
      gto.GaussW = params[1];
    } else if (std::string(argv[opt]) == "-l") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -l, expected 2." << std::endl;
        throw;
      }
      gto.FitBetween = std::make_pair(params[0], params[1]);
    } else if (std::string(argv[opt]) == "-C") {
      gto.UseNuPrismChi2 = true;
    }
    opt++;
  }

  return gto;
}
