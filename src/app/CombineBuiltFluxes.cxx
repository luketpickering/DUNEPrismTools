#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "TFile.h"
#include "TRegexp.h"
#include "TString.h"

// Unix
#include <dirent.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

std::string runPlanCfg = "", runPlanName = "";
std::vector<std::string> input_patterns;
std::string outputfile = "";

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << std::endl;
  std::cout << "\t-f <RunPlan[,RunPlanName]> : Run plan configuration xml, "
               "optionally with a specified plan. If none is specified, the "
               "first one in the file is used."
            << std::endl;
  std::cout << "\t-i <Input search pattern>  : Search pattern to find input "
               "files. Can be specified multiple times."
            << std::endl;
  std::cout
      << "\t-o <Output file name>      : File to write combined output to."
      << std::endl;
  std::cout << "\t-?                         : Display this message."
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
    } else if (std::string(argv[opt]) == "-i") {
      input_patterns.push_back(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      outputfile = argv[++opt];
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

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!runPlanCfg.size()) {
    std::cout << "[ERROR]: Expected to find an input run plan, -r" << std::endl;
    SayUsage(argv);
    return 1;
  }

  if (!input_patterns.size()) {
    std::cout << "[ERROR]: Expected to find an input file pattern, -i"
              << std::endl;
    SayUsage(argv);
    return 1;
  }

  if (!outputfile.size()) {
    std::cout << "[ERROR]: Expected to find an out file, -o" << std::endl;
    SayUsage(argv);
    return 1;
  }

  std::vector<DetectorStop> detStops =
      ReadDetectorStopConfig(runPlanCfg, runPlanName);

  if (!detStops.size()) {
    std::cout << "[ERROR]: Could not read any detector stops from the run plan."
              << std::endl;
    SayUsage(argv);
    return 1;
  }

  size_t NFilesAdded = 0;

  for (size_t ip_it = 0; ip_it < input_patterns.size(); ++ip_it) {
    std::string input_pattern = input_patterns[ip_it];
    size_t AsteriskPos = input_pattern.find_last_of('*');
    if (AsteriskPos == std::string::npos) {
      std::cout << "[ERROR]: Expected to find a wildcard in the argument of "
                   "the -i parameter. This executable just joins the output of "
                   "dp_BuildFluxes, if you only have a single input, something "
                   "is wrong."
                << std::endl;
      return 1;
    }

    DIR *dir;
    struct dirent *ent;
    size_t lastFSlash = input_pattern.find_last_of('/');
    std::string matchPat = input_pattern.substr(lastFSlash + 1);
    std::string dirpath = "";
    if (lastFSlash == std::string::npos) {
      char *cwd = new char[1000];
      getcwd(cwd, sizeof(char) * 1000);
      std::cout << "\t--Looking in current directory (" << cwd
                << ") for matching (\"" << matchPat << "\") files."
                << std::endl;
      dirpath = "./";
      delete cwd;
    } else {
      if (AsteriskPos < lastFSlash) {
        std::cout << "[ERROR]: Currently cannot handle a wildcard in the "
                     "directory structure. Please put input files in the same "
                     "directory. Expected -i \"../some/rel/path/*.root\""
                  << std::endl;
        return 1;
      }
      dirpath = input_pattern.substr(0, lastFSlash + 1);
      std::cout << "\t--Looking in directory (" << dirpath
                << ") for matching files." << std::endl;
    }
    dir = opendir(dirpath.c_str());

    if (dir != NULL) {
      TRegexp matchExp(matchPat.c_str(), true);
      /* print all the files and directories within directory */
      Ssiz_t len = 0;
      while ((ent = readdir(dir)) != NULL) {
        if (matchExp.Index(TString(ent->d_name), &len) != Ssiz_t(-1)) {
          std::cout << "\t\t\tAdding matching file: " << ent->d_name
                    << std::endl;

          TFile *ifl = CheckOpenFile(dirpath + ent->d_name);

          for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
            detStops[ds_it].Read(ifl, true);
          }

          NFilesAdded++;
          ifl->Close();
          delete ifl;
        }
      }
      closedir(dir);
    } else {
      /* could not open directory */
      perror("");
      return false;
    }
  }

  std::cout << "[INFO]: Added " << NFilesAdded << " files." << std::endl;

  TFile *ofl = CheckOpenFile(outputfile, "RECREATE");
  ofl->cd();
  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    detStops[ds_it].Write();
  }

  std::stringstream ss("");

  for (int species : {-14, -12, 12, 14}) {
    std::vector<TH1D *> Fluxes;
    std::vector<double> OAABinEdges;

    for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
      for (size_t ms_it = 0; ms_it < detStops[ds_it].GetNMeasurementSlices();
           ++ms_it) {
        double BinLowEdge = detStops[ds_it].GetAbsoluteOffsetOfSlice(ms_it) -
                            (detStops[ds_it].MeasurementRegionWidth / 2.0);

        if (!OAABinEdges.size() ||
            (fabs(OAABinEdges.back() - BinLowEdge) > 1E-5)) {
          OAABinEdges.push_back(BinLowEdge);
        }

        Fluxes.push_back(detStops[ds_it].GetFluxForSpecies(ms_it, species));
      }
      double BinUpEdge = (detStops[ds_it].GetAbsoluteOffsetOfSlice(
                              detStops[ds_it].GetNMeasurementSlices() - 1) +
                          (detStops[ds_it].MeasurementRegionWidth / 2.0));

      OAABinEdges.push_back(BinUpEdge);
    }

    if (!Fluxes.size()) {
      continue;
    }

    if ((Fluxes.size() + 1) != OAABinEdges.size()) {
      std::cout << "[ERROR]: Found " << Fluxes.size() << " and built "
                << OAABinEdges.size() << " flux bin edges. Incompatible."
                << std::endl;
      throw;
    }

    Int_t NXBins = Fluxes.front()->GetXaxis()->GetNbins();
    Int_t NYBins = Fluxes.size();

    double *XBins = new double[NXBins + 1];

    std::cout << "[INFO]: Building 2D histo with " << NYBins
              << " off axis slices." << std::endl;

    for (Int_t enu_bi_it = 0; enu_bi_it < NXBins; ++enu_bi_it) {
      XBins[enu_bi_it] =
          Fluxes.front()->GetXaxis()->GetBinLowEdge(enu_bi_it + 1);
    }
    XBins[NXBins] = Fluxes.front()->GetXaxis()->GetBinUpEdge(NXBins);

    ss.str("");
    ss << GetSpeciesName(species) << "_flux_2D";

    TH2D *ENuOAA =
        new TH2D(ss.str().c_str(),
                 ";Off-axis Position (m);#it{E}_{#nu} (GeV);#Phi_{#nu} "
                 "(m^{-2} per 1 GeV per 1 POT)",
                 NXBins, XBins, NYBins, OAABinEdges.data());

    Int_t dum = 0;
    for (Int_t fl_it = 0; fl_it < NYBins; ++fl_it) {
      for (Int_t enu_bi_it = 0; enu_bi_it < NXBins; ++enu_bi_it) {
        Int_t GBin = ENuOAA->GetBin(enu_bi_it + 1, fl_it + 1, dum);
        ENuOAA->SetBinContent(GBin,
                              Fluxes[fl_it]->GetBinContent(enu_bi_it + 1));
        ENuOAA->SetBinError(GBin, Fluxes[fl_it]->GetBinError(enu_bi_it + 1));
      }
    }

    ENuOAA->SetName(ss.str().c_str());
    ENuOAA->Write(ENuOAA->GetName(), TObject::kOverwrite);
    ENuOAA->SetDirectory(nullptr);

    delete[] XBins;

    delete ENuOAA;
  }

  ofl->Write();
  ofl->Close();
  delete ofl;
}
