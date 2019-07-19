#ifndef STRINGPARSERUTILITY_HXX_SEEN
#define STRINGPARSERUTILITY_HXX_SEEN

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template <typename T> inline T str2T(std::string const &str) {
  std::istringstream stream(str);
  T d;
  stream >> d;

  if (stream.fail()) {
    std::cerr << "[WARN]: Failed to parse string: " << str
              << " as requested type." << std::endl;
    return T();
  }

  return d;
}

template <typename T> inline std::string to_str(T const &inp) {
  std::stringstream stream("");
  stream << inp;
  return stream.str();
}

template <> inline bool str2T<bool>(std::string const &str) {
  if ((str == "true") || (str == "True") || (str == "TRUE") || (str == "1")) {
    return true;
  }

  if ((str == "false") || (str == "False") || (str == "FALSE") ||
      (str == "0")) {
    return false;
  }

  std::istringstream stream(str);
  bool d;
  stream >> d;

  if (stream.fail()) {
    std::cerr << "[WARN]: Failed to parse string: " << str
              << " as requested type." << std::endl;
    return false;
  }

  return d;
}

template <typename T>
inline std::vector<T> ParseToVect(std::string const &str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<T> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) {
      continue;
    }
    vals.push_back(str2T<T>(temp_string));
  }
  return vals;
}

template <typename T>
inline std::vector<T> strsToTs(std::vector<std::string> const &strs) {
  std::vector<T> vals;

  for (std::string const &s : strs) {
    if (s.empty()) {
      continue;
    }
    vals.push_back(str2T<T>(s));
  }
  return vals;
}

template <typename T>
inline void AppendVect(std::vector<T> &target, std::vector<T> const &toApp) {
  for (size_t i = 0; i < toApp.size(); ++i) {
    target.push_back(toApp[i]);
  }
}

// Converts "5_10:1" into a vector containing: 5,6,7,8,9,10
std::vector<double> BuildDoubleList(std::string const &str);
// Converts "1,5_10:1,15" into a vector containing: 1,5,6,7,8,9,10,15
std::vector<double> ParseDoubleListDescriptor(std::string const &str);
// Converts "1_5:2,5_10:1,15" into a vector containing: 1,3,5,6,7,8,9,10,15
std::vector<double> BuildBinEdges(std::string const &str);
// Takes 1,1,1,2,3,4,5,5,6 to 1,2,3,4,5,6 and throws on 1,2,1,3 (bins must be in
// increasing order)
std::vector<double> SanitizeBinEdges(std::vector<double> varbins);
std::vector<std::pair<double, double>> BuildRangesList(std::string const &str);
void chomp(std::string &str);

inline std::string str_replace(std::string const &inp, std::string const &from,
                               std::string const &to) {
  std::stringstream ss("");

  size_t nextOccurence = 0;
  size_t prevOccurence = 0;
  bool AtEnd = false;
  while (!AtEnd) {
    nextOccurence = inp.find(from, prevOccurence);
    if (nextOccurence == std::string::npos) {
      if (prevOccurence == inp.length()) {
        break;
      }
      AtEnd = true;
    }
    if ((nextOccurence != prevOccurence) || (nextOccurence == 0)) {
      ss << inp.substr(prevOccurence, (nextOccurence - prevOccurence));
      if (!AtEnd) {
        ss << to;
      }
    }
    prevOccurence = nextOccurence + from.size();
  }
  return ss.str();
}

#endif
