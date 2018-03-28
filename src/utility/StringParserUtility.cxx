#include "StringParserUtility.hxx"

// Converts "5_10:1" into a vector containing: 5,6,7,8,9,10
std::vector<double> BuildDoubleList(std::string const &str) {
  std::vector<std::string> steps = ParseToVect<std::string>(str, ":");
  if (steps.size() != 2) {
    return ParseToVect<double>(str, ",");
  }
  double step = str2T<double>(steps[1]);

  std::vector<double> range = ParseToVect<double>(steps[0], "_");
  if (steps.size() != 2) {
    std::cout
        << "[ERROR]: When attempting to parse bin range descriptor: \" " << str
        << "\", couldn't determine range. Expect form: <bin1low>_<binXUp>:step"
        << std::endl;
    exit(1);
  }

  int nsteps = (range[1] - range[0]) / step;

  std::vector<double> rtn;
  for (int step_it = 0; step_it <= nsteps; ++step_it) {
    rtn.push_back(range[0] + step * step_it);
  }
  return rtn;
}
std::vector<double> ParseDoubleListDescriptor(std::string const &str){
  std::vector<std::string> splitDescriptors =
      ParseToVect<std::string>(str, ",");
  std::vector<double> list;
  for (size_t vbd_it = 0; vbd_it < splitDescriptors.size(); ++vbd_it) {
    AppendVect(list, BuildDoubleList(splitDescriptors[vbd_it]));
  }
  return list;
}

std::vector<double> BuildBinEdges(std::string const &str){
  std::vector<double> varbins = ParseDoubleListDescriptor(str);

  for (size_t bin_it = 1; bin_it < varbins.size(); ++bin_it) {
    if (varbins[bin_it] == varbins[bin_it - 1]) {
      std::cout << "[INFO]: Removing duplciate bin low edge " << bin_it
                << " low edge: " << varbins[bin_it] << std::endl;
      varbins.erase(varbins.begin() + bin_it);
    }
  }

  for (size_t bin_it = 1; bin_it < varbins.size(); ++bin_it) {
    if (varbins[bin_it] < varbins[bin_it - 1]) {
      std::cout << "[ERROR]: Bin " << bin_it
                << " low edge: " << varbins[bin_it]
                << " is smaller than bin " << (bin_it - 1)
                << " low edge: " << varbins[bin_it - 1] << std::endl;
      throw;
    }
  }
  return varbins;
}

void chomp(std::string &str) {
  size_t lnf = str.find_last_not_of("\r\n");
  if (lnf != std::string::npos) {
    str = str.substr(0, lnf + 1);
  }
}
