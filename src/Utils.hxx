#ifndef DP_UTILS_HXX_SEEN
#define DP_UTILS_HXX_SEEN

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

inline double str2d(std::string str) {
  std::istringstream stream(str);
  double d;
  stream >> d;
  return d;
}

inline int str2i(std::string str) {
  std::istringstream stream(str);
  int d;
  stream >> d;
  return d;
}

inline std::vector<double> ParseToDbl(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<double> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    std::istringstream stream(temp_string);
    double entry;
    stream >> entry;

    vals.push_back(entry);
  }
  return vals;
}

inline std::vector<std::string> ParseToStr(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<std::string> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    vals.push_back(temp_string);
  }

  return vals;
}

inline void AppendDblVect(std::vector<double> &target,
                   std::vector<double> const &toApp) {
  for (size_t i = 0; i < toApp.size(); ++i) {
    target.push_back(toApp[i]);
  }
}

// Converts "5_10:1" into a vector containing: 5,6,7,8,9,10
inline std::vector<double> BuildDoubleList(std::string const &str) {
  std::vector<std::string> steps = ParseToStr(str, ":");
  if (steps.size() != 2) {
    return ParseToDbl(str, ",");
  }
  double step = str2d(steps[1]);

  std::vector<double> range = ParseToDbl(steps[0], "_");
  if (!steps.size() == 2) {
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
#endif
