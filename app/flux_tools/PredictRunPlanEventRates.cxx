#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "BargerPropagator.h"

#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

std::vector<DetectorStop> detStops;
std::string runPlanCfg, runPlanName = "", FluxesFile;
std::vector<std::pair<std::string, std::string> > XSecComponentInputs;
std::map<std::string, TH1D *> XSecComponents;

double DipAngle = 5.8;
double OscParams[6] = {0.825, 0.10, 1.0, 7.9e-5, 2.5e-3, 0.0};

int NuSpecies = 14;
int OscNuSpecies = 0;

std::string oupTableFile, oupFile;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-f <ROOT file>                     : The input flux "
               "histograms.\n"
               "\n"
               "\t-r <RunPlan.XML[,Plan Name]>       : An XML file specifying "
               "a run plan  \n"
               "\t                                     to make event rate "
               "predictions and \n"
               "\t                                     measured spectra for. "
               "See          \n"
               "\t                                     documentation for XML "
               "structure.   \n"
               "\n"
               "\t-o <ROOT file>                     : The output root file. \n"
               "\n"
               "\t-t <table file>                    : Output event rate table."
               "\n\n"
               "\n"
               "\t-x <ROOT file, hist name>          : Add xsec component for "
               "making event\n"
               "\t                                     rate predictions.       "
               "           \n"

               "\t-p "
               "<sin2(theta12)=0.825>,<sin2(theta13)=0.10>,<sin2(theta23)=1.0>,"
               "<dm12=7.9e-5>,<dm23=2.5e-3>,<dcp=0.0>\n"
               "\t                                    : Oscillation parameters."
               "By default \n"
               "\t                                     no oscillations occur.\n"
               "\n"
               "\t-d <dipangle=5.8>                   : Beam dip angle. Used to"
               " calculate path\n"
               "\t                                     length for oscillation."
               "\n\n"
               "\t-n <neutrino PDG=14>                : Species to run.\n"
               "\t-u <neutrino PDG=14>                : Species to oscillate to"
               ".\n"
               "\n"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-t") {
      oupTableFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-f") {
      FluxesFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      oupFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
    } else if (std::string(argv[opt]) == "-x") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() < 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -x, expected at least 2." << std::endl;
        exit(1);
      }
      for (size_t xs_it = 1; xs_it < params.size(); ++xs_it) {
        XSecComponentInputs.push_back(std::make_pair(params[0], params[xs_it]));
      }
    } else if (std::string(argv[opt]) == "-d") {
      DipAngle = str2T<double>(argv[++opt]);

      const static double REarth_cm = 6371.393 * 1.0E5;
      const static double ProductionHeight_cm = 0;
      static const double deg2rad = asin(1) / 90.0;
      double cz = cos((90.0 + DipAngle) * deg2rad);
      double PathLength = sqrt((REarth_cm + ProductionHeight_cm) *
                                   (REarth_cm + ProductionHeight_cm) -
                               (REarth_cm * REarth_cm) * (1 - cz * cz)) -
                          REarth_cm * cz;

      std::cout << "Calculated path length: " << (PathLength / 1.0E5) << " km."
                << std::endl;

    } else if (std::string(argv[opt]) == "-p") {
      std::vector<double> params = ParseToVect<double>(argv[++opt], ",");
      if (params.size() != 6) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -p, expected 6." << std::endl;
        exit(1);
      }

      OscParams[0] = params[0];
      OscParams[1] = params[1];
      OscParams[2] = params[2];
      OscParams[3] = params[3];
      OscParams[4] = params[4];
      OscParams[5] = params[5];

      std::cout << "Sin^2(Theta_12) = " << OscParams[0] << std::endl;
      std::cout << "Sin^2(Theta_13) = " << OscParams[1] << std::endl;
      std::cout << "Sin^2(Theta_23) = " << OscParams[2] << std::endl;

      std::cout << "Dm^2_21 = " << OscParams[3] << std::endl;
      std::cout << "|Dm^2_Atm| = " << OscParams[4] << std::endl;

      std::cout << "dcp = " << OscParams[5] << std::endl;

    } else if (std::string(argv[opt]) == "-n") {
      NuSpecies = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-u") {
      OscNuSpecies = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

enum nuTypes {
  kNuebarType = -1,
  kNumubarType = -2,
  kNutaubarType = -3,
  kNueType = 1,
  kNumuType = 2,
  kNutauType = 3,
};

nuTypes GetNuType(int pdg) {
  switch (pdg) {
    case 16:
      return kNutauType;
    case 14:
      return kNumuType;
    case 12:
      return kNueType;
    case -16:
      return kNutaubarType;
    case -14:
      return kNumubarType;
    case -12:
      return kNuebarType;
    default: {
      std::cout << "[ERROR]: Attempting to convert \"neutrino pdg\": " << pdg
                << std::endl;
      exit(1);
    }
  }
}

double OscWeight(double enu) {
  BargerPropagator bp;
  int NuType = GetNuType(NuSpecies);
  bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
            OscParams[4], OscParams[5], enu, true, NuType);

  static const double deg2rad = asin(1) / 90.0;
  double lengthParam = cos((90.0 + DipAngle) * deg2rad);
  bp.DefinePath(lengthParam, 0);
  bp.propagate(NuType);
  return bp.GetProb(NuType, GetNuType(OscNuSpecies));
}

void Oscillate(TH1D *in) {
  for (Int_t bi_it = 1; bi_it < in->GetXaxis()->GetNbins() + 1; ++bi_it) {
    double ow = OscWeight(in->GetXaxis()->GetBinCenter(bi_it));
    in->SetBinContent(bi_it, in->GetBinContent(bi_it) * ow);
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!runPlanCfg.size()) {
    std::cout << "[ERROR]: Expected to find an input run plan, -r" << std::endl;
    SayUsage(argv);
    return 1;
  }

  detStops = ReadDetectorStopConfig(runPlanCfg, runPlanName);

  if (!detStops.size()) {
    std::cout << "[ERROR]: Could not read any detector stops from the run plan."
              << std::endl;
    SayUsage(argv);
    return 1;
  }

  if (!FluxesFile.size()) {
    std::cout << "[ERROR]: Expected to find an input flux file, -f."
              << std::endl;
    SayUsage(argv);
    return 1;
  }

  if (!oupFile.size()) {
    std::cout << "[ERROR]: Expected to find a output file, -o." << std::endl;
    SayUsage(argv);
    return 1;
  }

  if (!XSecComponentInputs.size()) {
    std::cout
        << "[ERROR]: Expected to find at least one input cross-section, -x"
        << std::endl;
    SayUsage(argv);
    return 1;
  }

  TFile *ifl = CheckOpenFile(FluxesFile.c_str());

  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    detStops[ds_it].Read(ifl);
  }

  for (std::pair<std::string, std::string> hdescript : XSecComponentInputs) {
    TH1D *hist = GetHistogram<TH1D>(hdescript.first, hdescript.second);
    std::cout << "[INFO]: Got XSec component: " << hist->GetName() << std::endl;
    XSecComponents[hist->GetName()] = static_cast<TH1D *>(hist->Clone());
    XSecComponents[hist->GetName()]->SetDirectory(nullptr);
  }

  TFile *ofl = new TFile(oupFile.c_str(), "RECREATE");

  std::vector<double> EvrAll = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double TPOT = 0;
  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    if (OscNuSpecies) {
      for (size_t ms_it = 0; ms_it < detStops[ds_it].GetNMeasurementSlices();
           ++ms_it) {
        Oscillate(detStops[ds_it].GetFluxForSpecies(ms_it, NuSpecies));
      }
    }

    detStops[ds_it].PredictEventRates(XSecComponents, NuSpecies);
    detStops[ds_it].Write();

    if (!ds_it) {
      std::cout << "Detector fiducial mass: " << (detStops[ds_it].Mass() / 1.E3)
                << " T LAr." << std::endl;
      std::cout
          << "Offset & POT/$1\\times10^{21}$ & CCInc & CC0$\\pi$ / CCInc & "
             "CC1$\\pi$ / CCInc & CCN$\\pi$ / CCInc & "
             "NCInc & $\\nu-$e El. "
             "\\\\"
          << std::endl;
    }

    std::vector<double> Evr = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (size_t ms_it = 0; ms_it < detStops[ds_it].GetNMeasurementSlices();
         ++ms_it) {
      Evr[0] += detStops[ds_it]
                    .GetPredictedEventRate("CC_0pi_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[4] += detStops[ds_it]
                    .GetPredictedEventRate("CC_0pi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[1] += detStops[ds_it]
                    .GetPredictedEventRate("CC_1cpi_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[4] += detStops[ds_it]
                    .GetPredictedEventRate("CC_1cpi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[2] += detStops[ds_it]
                    .GetPredictedEventRate("CC_1pi0_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[4] += detStops[ds_it]
                    .GetPredictedEventRate("CC_1pi0_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[3] += detStops[ds_it]
                    .GetPredictedEventRate("CC_Npi_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[4] += detStops[ds_it]
                    .GetPredictedEventRate("CC_Npi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[5] += detStops[ds_it]
                    .GetPredictedEventRate("NC_El_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[9] += detStops[ds_it]
                    .GetPredictedEventRate("NC_El_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[6] += detStops[ds_it]
                    .GetPredictedEventRate("NC_1cpi_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[9] += detStops[ds_it]
                    .GetPredictedEventRate("NC_1cpi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[7] += detStops[ds_it]
                    .GetPredictedEventRate("NC_1pi0_xsec", NuSpecies, ms_it)
                    ->Integral("width");
      Evr[9] += detStops[ds_it]
                    .GetPredictedEventRate("NC_1pi0_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[8] += detStops[ds_it]
                    .GetPredictedEventRate("NC_Npi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      Evr[9] += detStops[ds_it]
                    .GetPredictedEventRate("NC_Npi_xsec", NuSpecies, ms_it)
                    ->Integral("width");

      if (abs(NuSpecies) == 14) {
        Evr[10] += detStops[ds_it]
                       .GetPredictedEventRate("nu_e_el_xsec", NuSpecies, ms_it)
                       ->Integral("width");
      } else {
        Evr[10] += detStops[ds_it]
                       .GetPredictedEventRate("M_xsec", NuSpecies, ms_it)
                       ->Integral("width");
      }
    }

    std::cout << std::setprecision(3) << detStops[ds_it].LateralOffset
              << "~m & " << (detStops[ds_it].POTExposure / 1.E21) << " & "
              << Evr[4] << " & " << Evr[0] / Evr[4] << " & "
              << (Evr[1] + Evr[2]) / Evr[4] << " & " << Evr[3] / Evr[4] << " & "
              << Evr[9] << " & " << Evr[10] << " \\\\" << std::endl;

    TPOT += (detStops[ds_it].POTExposure / 1.E21);

    for (size_t i = 0; i < 11; ++i) {
      EvrAll[i] += Evr[i];
    }
  }

  std::cout << std::setprecision(3) << "All"
            << " & " << TPOT << " & " << EvrAll[4] << " & "
            << EvrAll[0] / EvrAll[4] << " & "
            << (EvrAll[1] + EvrAll[2]) / EvrAll[4] << " & "
            << EvrAll[3] / EvrAll[4] << " & " << EvrAll[9] << " & "
            << EvrAll[10] << " \\\\" << std::endl;
  ofl->Write();
  ofl->Close();
  delete ofl;
}
