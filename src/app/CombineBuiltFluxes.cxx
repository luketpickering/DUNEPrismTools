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

std::string runPlanCfg = "", runPlanName = "";
std::string input_pattern = "";
std::string outputfile = "";

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << std::endl;
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
      input_pattern = argv[++opt];
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

  if (!input_pattern.size()) {
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
              << ") for matching (\"" << matchPat << "\") files." << std::endl;
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
    size_t NFilesAdded = 0;
    while ((ent = readdir(dir)) != NULL) {
      if (matchExp.Index(TString(ent->d_name), &len) != Ssiz_t(-1)) {
        std::cout << "\t\t\tAdding matching file: " << ent->d_name << std::endl;

        TFile *ifl = CheckOpenFile(ent->d_name);

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

  TFile *ofl = CheckOpenFile(outputfile, "RECREATE");
  ofl->cd();
  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    detStops[ds_it].Write();
  }
  ofl->Write();
  ofl->Close();
  delete ofl;
}
