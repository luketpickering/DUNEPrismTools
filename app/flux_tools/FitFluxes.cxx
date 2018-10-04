#include "FluxFitResultsTreeReader.hxx"
#include "OscillationParametersTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"

#include "GetUsage.hxx"

#include "FluxFitter.hxx"

bool IsGauss = false;

std::string OutputFile, OutputDirectory = "";
bool UPDATEOutputFile = false;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << GetUsageText(argv[0], "flux_tools")
            << std::endl;
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
    } else if ((std::string(argv[opt]) == "-?") ||
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

  //
  // TDirectory *wD = oupD->mkdir("weighted_fluxes");
  // wD->cd();
  //
  // for (size_t flux_it = 0; flux_it < Fluxes.size(); flux_it++) {
  //   TH1D *cl = static_cast<TH1D *>(Fluxes[flux_it]->Clone());
  //   cl->Scale(coeffs[flux_it]);
  //   cl->SetDirectory(wD);
  // }
  //
  // if (minimizer) {
  //   std::pair<std::vector<double>, std::vector<double>> XRangeBinCoeffs =
  //       SliceConfig::BuildXRangeBinsCoeffs(XRanges, coeffs);
  //   TH1D *coeffsH = new TH1D("Coeffs", "Coeffs;Off-axis position;Weight",
  //                            (XRangeBinCoeffs.first.size() - 1),
  //                            XRangeBinCoeffs.first.data());
  //
  //   for (size_t coeff_it = 0; coeff_it < (XRangeBinCoeffs.first.size() - 1);
  //        coeff_it++) {
  //     coeffsH->SetBinContent(coeff_it + 1, XRangeBinCoeffs.second[coeff_it]);
  //     coeffsH->SetBinError(coeff_it + 1, 0);
  //   }
  // }
  //
  // oupD->cd();
  //
  // SummedFlux->SetDirectory(oupD);
  // SummedFlux->SetName("BestFit");
  // if (IsGauss) {
  //   TGraph *tGauss = new TGraph(1);
  //   tGauss->Set(1E4 - 1);
  //
  //   double min = SummedFlux->GetXaxis()->GetBinLowEdge(1);
  //   double step = (SummedFlux->GetXaxis()->GetBinUpEdge(
  //                      SummedFlux->GetXaxis()->GetNbins()) -
  //                  SummedFlux->GetXaxis()->GetBinLowEdge(1)) /
  //                 double(1E4);
  //   for (size_t i = 1; i < 1E4; ++i) {
  //     double enu = min + i * step;
  //     double val = TargetGauss->Eval(enu);
  //     if (val != val) {
  //       continue;
  //     }
  //     tGauss->SetPoint(i - 1, enu, val);
  //   }
  //
  //   tGauss->Write("Target_gauss");
  //
  // } else {
  //   TargetFlux->SetDirectory(oupD);
  //   TargetFlux->SetName("Target");
  //
  //   TargetFlux_orignorm->SetDirectory(oupD);
  //   TargetFlux_orignorm->SetName("Target_input_normalisation");
  //
  //   OscFlux->SetDirectory(oupD);
  //   OscFlux->SetName("InputFlux");
  //
  //   OscFlux_orignorm->SetDirectory(oupD);
  //   OscFlux_orignorm->SetName("InputFlux_input_normalisation");
  //
  //   BestFit_FDNorm = static_cast<TH1D
  //   *>(SummedFlux->Clone("BestFit_FDNorm")); BestFit_FDNorm->Scale(1.0 /
  //   NDOverFDFitScaleFactor); BestFit_FDNorm->SetDirectory(oupD);
  // }

  // if (fitstatus > 0) {
  //   std::cout << "[WARN]: Failed to find minimum (STATUS: " << fitstatus <<
  //   ")."
  //             << std::endl;
  // }
  //
  // // Write out configuration trees.
  // SliceConfig *sc = SliceConfig::MakeTreeWriter();
  // for (size_t i = 0; i < Fluxes.size(); ++i) {
  //   sc->XRange[0] =
  //       XRanges[i].first * 100.0; // Assumes tree written in off-axis
  //       position
  //   sc->XRange[1] =
  //       XRanges[i].second * 100.0; // Assumes tree written in off-axis
  //       position
  //   sc->Coeff = coeffs[i];
  //
  //   sc->Fill();
  // }
  //
  // FluxFitResultsTreeReader *fr =
  //     FluxFitResultsTreeReader::MakeTreeWriter(IsGauss);
  //
  // fr->NFluxes = Fluxes.size();
  // fr->NIterations = NCalls;
  // fr->Chi2 = Chi2_last;
  // fr->RegularisationPenalty = reg_last;
  // fr->FitRange[0] = FitBetween_low;
  // fr->FitRange[1] = FitBetween_high;
  //
  // if (!IsGauss) {
  //   fr->OutOfRangePenalty = OOR_last;
  //   fr->NDOverFDFitScaleFactor = NDOverFDFitScaleFactor;
  //   // (SummedFlux->GetMaximum() / OscFlux_orignorm->GetMaximum());
  // } else {
  //   fr->GaussCenter_GeV = GaussC;
  //   fr->GaussWidth_GeV = GaussW;
  // }
  // fr->Fill();
  //
  // if (!IsGauss) {
  //   try {
  //     OscillationParameters opTree(InputTargetFluxFile);
  //     TTree *opTree_copy = opTree.tree->CloneTree();
  //     opTree_copy->SetDirectory(oupF);
  //   } catch (...) {
  //     std::cout
  //         << "[INFO]: Not fitting oscillated flux, ignoring summary tree
  //         copy."
  //         << std::endl;
  //   }
  // }
  //
  oupF->Write();
  oupF->Close();
}
