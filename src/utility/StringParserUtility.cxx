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
std::vector<double> ParseDoubleListDescriptor(std::string const &str) {
  std::vector<std::string> splitDescriptors =
      ParseToVect<std::string>(str, ",");
  std::vector<double> list;
  for (size_t vbd_it = 0; vbd_it < splitDescriptors.size(); ++vbd_it) {
    AppendVect(list, BuildDoubleList(splitDescriptors[vbd_it]));
  }
  return list;
}

std::vector<double> SanitizeBinEdges(std::vector<double> varbins) {
  for (size_t bin_it = 1; bin_it < varbins.size(); ++bin_it) {
    if (varbins[bin_it] == varbins[bin_it - 1]) {
      varbins.erase(varbins.begin() + bin_it);
    }
  }

  for (size_t bin_it = 1; bin_it < varbins.size(); ++bin_it) {
    if (varbins[bin_it] < varbins[bin_it - 1]) {
      std::cout << "[ERROR]: Bin " << bin_it << " low edge: " << varbins[bin_it]
                << " is smaller than bin " << (bin_it - 1)
                << " low edge: " << varbins[bin_it - 1] << std::endl;
      throw;
    }
  }
  return varbins;
}

std::vector<double> BuildBinEdges(std::string const &str) {
  return SanitizeBinEdges(ParseDoubleListDescriptor(str));
}

std::vector<std::pair<double, double>> BuildRangesList(std::string const &str) {
  std::vector<std::string> listDescriptor = ParseToVect<std::string>(str, ",");
  std::vector<std::pair<double, double>> RangesList;

  for (size_t l_it = 0; l_it < listDescriptor.size(); ++l_it) {
    // If this includes a range to build
    if (listDescriptor[l_it].find("_") != std::string::npos) {
      std::vector<std::string> rangeDescriptor =
          ParseToVect<std::string>(listDescriptor[l_it], ":");

      if (rangeDescriptor.size() != 2) {
        std::cout
            << "[ERROR]: Range descriptor: \"" << str
            << "\" contained bad descriptor: \"" << listDescriptor[l_it]
            << "\", expected <RangeCenterLow>_<RangeCenterHigh>:<RangeWidths>."
            << std::endl;
        exit(0);
      }

      std::vector<double> rangeCenters = BuildDoubleList(listDescriptor[l_it]);
      double width = str2T<double>(rangeDescriptor[1]);

      for (size_t sp_it = 0; sp_it < rangeCenters.size(); ++sp_it) {
        RangesList.push_back(
            std::make_pair(rangeCenters[sp_it] - (width / 2.0),
                           rangeCenters[sp_it] + (width / 2.0)));
      }

    } else {
      std::vector<double> rangeDescriptor =
          ParseToVect<double>(listDescriptor[l_it], ":");
      if (rangeDescriptor.size() != 2) {
        std::cout << "[ERROR]: Range descriptor: \"" << str
                  << "\" contained bad descriptor: \"" << listDescriptor[l_it]
                  << "\", expected <RangeCenter>:<RangeWidth>." << std::endl;
        exit(0);
      }
      RangesList.push_back(
          std::make_pair(rangeDescriptor[0] - (rangeDescriptor[1] / 2.0),
                         rangeDescriptor[0] + (rangeDescriptor[1] / 2.0)));
    }
  }

  for (size_t range_it = 1; range_it < RangesList.size(); ++range_it) {
    if ((RangesList[range_it - 1].second - RangesList[range_it].first) > 1E-5) {
      std::cout << "[ERROR]: Range #" << range_it << " = {"
                << RangesList[range_it].first << " -- "
                << RangesList[range_it].second << "}, but #" << (range_it - 1)
                << " = {" << RangesList[range_it - 1].first << " -- "
                << RangesList[range_it - 1].second << "}." << std::endl;
      exit(1);
    }
  }
  return RangesList;
}

void chomp(std::string &str) {
  size_t lnf = str.find_last_not_of("\r\n");
  if (lnf != std::string::npos) {
    str = str.substr(0, lnf + 1);
  }
}
