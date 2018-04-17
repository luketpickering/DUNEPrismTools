#include "StopDimensions.hxx"

#include "DetectorStop.hxx"
#include "ROOTUtility.hxx"

#include "GetUsage.hxx"

#include <cmath>
#include <vector>

std::vector<double> StopVetoGap;
std::vector<double> ExtraVertexSelectionPadding;

std::string runPlanCfg = "", runPlanName = "";

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "sim_tools") << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
    } else if (std::string(argv[opt]) == "-V") {
      StopVetoGap = ParseToVect<double>(argv[++opt], ",");
      if (StopVetoGap.size() != 3) {
        std::cout << "[ERROR]: -V option contained " << StopVetoGap.size()
                  << " entries, expected 3." << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-FV") {
      ExtraVertexSelectionPadding = ParseToVect<double>(argv[++opt], ",");
      if (ExtraVertexSelectionPadding.size() != 3) {
        std::cout << "[ERROR]: -FV option contained "
                  << ExtraVertexSelectionPadding.size()
                  << " entries, expected 3." << std::endl;
        SayUsage(argv);
        exit(1);
      }
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

  if (!runPlanCfg.size()) {
    std::cout << "[INFO]: Expected to find at least -r option." << std::endl;
    SayUsage(argv);
    exit(1);
  }

  if (!ExtraVertexSelectionPadding.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {argv[0], "-FV", "0,0,0"};
    handleOpts(argc_dum, argv_dum);
  }

  if (!StopVetoGap.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {argv[0], "-V", "50,50,50"};
    handleOpts(argc_dum, argv_dum);
  }

  Double_t MaxAbsX = -std::numeric_limits<double>::max();
  Double_t MinAbsX = std::numeric_limits<double>::max();

  // Read in det stops
  std::vector<DetectorStop> DetectorStops =
      ReadDetectorStopConfig(runPlanCfg, runPlanName);
  std::vector<std::pair<double, double>> FVXRanges;

  for (size_t d_it = 0; d_it < DetectorStops.size(); ++d_it) {
    double XOffset = DetectorStops[d_it].CenterPosition[0] * 100.0;
    double XWidth_det = DetectorStops[d_it].ActiveExent[0] * 100.0;
    // double YWidth_det = DetectorStops[d_it].ActiveExent[1] * 100.0;
    // double ZWidth_det = DetectorStops[d_it].ActiveExent[2] * 100.0;
    double XWidth_nonveto = XWidth_det - 2 * StopVetoGap[0];
    // double YWidth_nonveto = YWidth_det - 2 * StopVetoGap[1];
    // double ZWidth_nonveto = ZWidth_det - 2 * StopVetoGap[2];

    double Detlow = XOffset - XWidth_det / 2.0;
    double DetHigh = XOffset + XWidth_det / 2.0;

    double X_Range_nonveto[2];
    X_Range_nonveto[0] = XOffset - XWidth_nonveto / 2.0;
    X_Range_nonveto[1] = XOffset + XWidth_nonveto / 2.0;

    // double Y_Range_nonveto[0] = -YWidth_nonveto / 2.0;
    // double Y_Range_nonveto[1] = YWidth_nonveto / 2.0;
    //
    // double Z_Range_nonveto[0] = -ZWidth_nonveto / 2.0;
    // double Z_Range_nonveto[1] = ZWidth_nonveto / 2.0;

    double X_Range_FV[2];
    X_Range_FV[0] =
        XOffset - (XWidth_nonveto / 2.0) + ExtraVertexSelectionPadding[0];
    X_Range_FV[1] =
        XOffset + (XWidth_nonveto / 2.0) - ExtraVertexSelectionPadding[0];

    // double Y_Range_FV[0] =
    //     -YWidth_nonveto + ExtraVertexSelectionPadding[1] / 2.0;
    // double Y_Range_FV[1] =
    //     YWidth_nonveto - ExtraVertexSelectionPadding[1] / 2.0;
    //
    // double Z_Range_FV[0] =
    //     -ZWidth_nonveto + ExtraVertexSelectionPadding[2] / 2.0;
    // double Z_Range_FV[1] =
    //     ZWidth_nonveto - ExtraVertexSelectionPadding[2] / 2.0;

    FVXRanges.push_back(std::make_pair(X_Range_FV[0], X_Range_FV[1]));

    MaxAbsX = std::max(MaxAbsX, DetHigh);
    MinAbsX = std::min(MinAbsX, Detlow);

    std::cout << "[INFO]: Det stop[" << d_it << "]: { Active: [" << Detlow
              << ", " << DetHigh << "], non-veto: [" << X_Range_nonveto[0]
              << ", " << X_Range_nonveto[1] << "], FV: [" << X_Range_FV[0]
              << ", " << X_Range_FV[1] << "] }" << std::endl;
  }

  size_t NSteps = lrint(MaxAbsX - MinAbsX) / 2;

  std::cout << "[INFO]: Checking coverage [ " << MinAbsX << " -- " << MaxAbsX
            << " : " << ((MaxAbsX - MinAbsX) / double(NSteps))
            << " ], NSteps = " << NSteps << std::endl;

  TH1D *XRangeCoverageHelper =
      new TH1D("XRangeCoverageHelper", "", NSteps, MinAbsX, MaxAbsX);
  XRangeCoverageHelper->SetDirectory(nullptr);

  for (std::pair<double, double> const &xr : FVXRanges) {

    for (Int_t bi_it = 0; bi_it < XRangeCoverageHelper->GetXaxis()->GetNbins();
         ++bi_it) {

      bool BinCenterIn =
          ((XRangeCoverageHelper->GetXaxis()->GetBinCenter(bi_it + 1) >
            xr.first) &&
           (XRangeCoverageHelper->GetXaxis()->GetBinCenter(bi_it + 1) <
            xr.second));

      if (BinCenterIn) {
        XRangeCoverageHelper->AddBinContent(bi_it + 1, 1);
      }
    }
  }

  bool FoundFirstWithContent = false;
  for (Int_t bi_it = 0; bi_it < XRangeCoverageHelper->GetXaxis()->GetNbins();
       ++bi_it) {

    if ((!FoundFirstWithContent) &&
        XRangeCoverageHelper->GetBinContent(bi_it + 1)) {
      FoundFirstWithContent = true;
      std::cout << "[INFO]: Lowest position with selected events @ "
                << XRangeCoverageHelper->GetXaxis()->GetBinLowEdge(bi_it + 1)
                << std::endl;
    }

    if (FoundFirstWithContent &&
        (!XRangeCoverageHelper->GetBinContent(bi_it + 1))) {
      std::cout << "[INFO]: Lowest position without coverage @ "
                << XRangeCoverageHelper->GetXaxis()->GetBinLowEdge(bi_it + 1)
                << std::endl;
      FoundFirstWithContent = false;
      break;
    }
  }
  if (FoundFirstWithContent) {
    std::cout << "[INFO]: Lowest position without coverage @ "
              << XRangeCoverageHelper->GetXaxis()->GetBinUpEdge(
                     XRangeCoverageHelper->GetXaxis()->GetNbins())
              << std::endl;
  }
}
