#include "ROOTUtility.hxx"
#include "PhysicsUtility.hxx"

#include "TFile.h"
#include "TRegexp.h"
#include "TString.h"

// Unix
#include <dirent.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

std::vector<std::string> input_patterns;
std::string outputfile = "";

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << std::endl;
  std::cout
      << "\t-i <Input search pattern>  : Search pattern to find input "
         "files. Can be specified multiple times.\n"
      << std::endl;
  std::cout
      << "\t-o <Output file name>      : File to write combined output to."
      << std::endl;
  std::cout
      << "\t-?                         : Display this message."
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-i") {
      input_patterns.push_back(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-?"||
               std::string(argv[opt]) == "--help") {
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

  std::vector<int> NuPDGTargets = {-12, 12, -14, 14};
  std::map< int, std::vector<TH1 *> > InputHists;

  size_t NFilesAdded = 0, NHistogramsAdded = 0;

  for (size_t ip_it = 0; ip_it < input_patterns.size(); ++ip_it) {
    std::string input_pattern = input_patterns[ip_it];
    size_t AsteriskPos = input_pattern.find_last_of('*');

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
        std::cout
            << "[ERROR]: Currently cannot handle a wildcard in the "
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

          // Read histograms
          for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
            std::stringstream hist_name("");
            hist_name << "LBNF_" << GetSpeciesName(NuPDGTargets[nuPDG_it]) << "_flux";
            TH1 * ih = dynamic_cast<TH1 *>(ifl->Get(hist_name.str().c_str()));
            if(!ih){
              continue;
            }
            if(!InputHists.count(NuPDGTargets[nuPDG_it])){
              InputHists[NuPDGTargets[nuPDG_it]] = std::vector<TH1 *>();
            }
            InputHists[NuPDGTargets[nuPDG_it]].push_back(static_cast<TH1 *>(ih->Clone()));
            InputHists[NuPDGTargets[nuPDG_it]].back()->SetDirectory(nullptr);
            NHistogramsAdded++;
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
  std::cout << "[INFO]: Added " << NFilesAdded << " files (" << NHistogramsAdded
    << " Histograms)." << std::endl;

  TFile *ofl = CheckOpenFile(outputfile, "RECREATE");
  ofl->cd();

  for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
    for (size_t hist_it = 1;
      hist_it < InputHists[NuPDGTargets[nuPDG_it]].size(); ++hist_it) {
      InputHists[NuPDGTargets[nuPDG_it]][0]->Add(
          InputHists[NuPDGTargets[nuPDG_it]][hist_it]);
    }
    if(InputHists[NuPDGTargets[nuPDG_it]].size()){
      InputHists[NuPDGTargets[nuPDG_it]][0]->Scale(
          1.0/double(InputHists[NuPDGTargets[nuPDG_it]].size()));
      InputHists[NuPDGTargets[nuPDG_it]][0]->SetDirectory(ofl);
    }
  }

  ofl->Write();
  ofl->Close();
  delete ofl;
}
