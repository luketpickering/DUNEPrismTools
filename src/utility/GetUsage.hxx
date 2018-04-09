#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>

// #define DEBUG

std::string GetUsageText(std::string const &appname, std::string const &appdir){

  char const *PRISMENV = getenv("DUNEPRISMTOOLSROOT");
  if(!PRISMENV){
    std::cout << "[ERROR]: No DUNEPRISMTOOLSROOT environment variable found. Is"
      " the environment set up correctly?" << std::endl;
    throw;
  }
  std::string penv(PRISMENV);
  if(penv.back() != '/'){ penv += "/"; }
  std::string filename = penv + "dox/" + appdir + ".md";
  std::ifstream mdfile (filename);

  if (!mdfile.is_open()){
    std::cout << "[WARN]: Couldn't find usage text for app " << appname
      << " from directory " << appdir << ". Looking for " << filename
      << std::endl;
    return "<<NO USAGE TEXT>>";
  }

  bool FoundAppTitle = false;
  bool OpeningTics = false;
  std::stringstream UsageText("");
  std::string str;
  while (std::getline(mdfile, str)) {
    if(!str.size()){
      continue;
    }
#ifdef DEBUG
      std::cout << "--" << str << std::endl;
#endif

    if((str[0] == '#') && (str.find(appname) != std::string::npos)){
      FoundAppTitle = true;
#ifdef DEBUG
      std::cout << "!Found AppTitle" << std::endl;
#endif
      continue;
    }

    if(FoundAppTitle){
      if(!OpeningTics && (str.find("```") == 0)){
        OpeningTics = true;
#ifdef DEBUG
        std::cout << "!Found opening tics" << std::endl;
#endif
        continue;
      }

      if(!OpeningTics){
        continue;
      }

      if(str.find("```") == 0){
#ifdef DEBUG
        std::cout << "!Found closing tics" << std::endl;
#endif
        break;
      }
      UsageText << str << std::endl;
    }
  }

  return UsageText.str();
}
