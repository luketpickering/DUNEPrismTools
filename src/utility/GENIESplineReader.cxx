#include "GENIESplineReader.hxx"

#include "XMLUtility.hxx"

#include <cctype>
#include <cstring>
#include <iostream>

SSAX::SSAX()
    : EventOnAllTags(true), LeavingAlready(false), InterestingTags(),
      buf_size(1), buf_usage(0), buf(NULL), depth(0) {
  buf = new char[buf_size];
};

SSAX::~SSAX() { delete buf; }

void SSAX::OnStartDocument() {}
void SSAX::OnEndDocument() {}
void SSAX::OnStartElement(
    const char *name,
    std::vector<std::pair<const char *, const char *>> &attributes) {
  std::cout << "[SSAX]: OnStartElement(\"" << name << "\", " << std::endl;
  for (size_t a_it = 0; a_it < attributes.size(); ++a_it) {
    std::cout << "\t\"" << attributes[a_it].first << "\" = \""
              << attributes[a_it].second << "\"" << std::endl;
  }
  std::cout << "\t)" << std::endl;
}
void SSAX::OnCharacters(const char *characters) {
  std::cout << "[SSAX]: OnCharacters(\"" << characters << "\")" << std::endl;
}
void SSAX::OnEndElement(const char *name) {
  std::cout << "[SSAX]: OnEndElement(\"" << name << "\")" << std::endl;
}
void SSAX::OnComment(const char *) {}

void SSAX::NotifyInterestingTags(std::vector<std::string> it) {
  InterestingTags = it;
  EventOnAllTags = !InterestingTags.size();
}
void SSAX::NotifyEarlyExit() { LeavingAlready = true; }

void SSAX::ExpandBuf() {
  size_t new_buf_size = (buf_size <= 1E6) ? buf_size * 10 : (buf_size + 1E6);

  std::cout << "[SSAX]: Expanding internal buffer to: " << new_buf_size
            << " bytes." << std::endl;

  char *new_buf = new char[new_buf_size];
  memcpy(new_buf, buf, buf_size);
  delete[] buf;

  buf = new_buf;
  buf_size = new_buf_size;
}

void SSAX::AddToBuf(char character) {
  if (buf_usage == buf_size) {
    ExpandBuf();
  }
  buf[buf_usage++] = character;
}

void SSAX::NullTermBuf() { AddToBuf('\0'); }

bool SSAX::IsNameStartChar(char c) {
  return (isalpha(c) || (c == ':') || (c == '_'));
}

bool SSAX::IsNameChar(char c) {
  return (IsNameStartChar(c) || isdigit(c) || (c == '-') || (c == '.'));
}

void SSAX::IgnoreAllToRightBracket(std::ifstream &ifs) {
  char prev, next;
  while (EOF != (next = ifs.get())) {
    if (next == '>') {
      if (prev == '/') { // Empty element
        depth--;
      }

      return;
    }
    prev = next;
  }
  std::cerr << "[SSAX - ERROR]: encountered EOF before ending '>' at depth: "
            << depth << std::endl;
  throw;
}

char SSAX::IgnoreSpace(std::ifstream &ifs) {
  char next;
  while (EOF != (next = ifs.get())) {
    if (isspace(next)) {
      continue;
    } else {
      return next;
    }
  }
  std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
            << std::endl;
  throw;
}

char SSAX::BufNameToBoundary(std::ifstream &ifs) {
  char next;
  while (EOF != (next = ifs.get())) {
    if (isspace(next) || (next == '>') || (next == '/') || (next == '=')) {
      return next;
    }
    if (IsNameChar(next)) {
      AddToBuf(next);
    } else {
      std::cerr << "[SSAX - ERROR]: encountered bad character " << next
                << std::endl;
      throw;
    }
  }
  std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
            << std::endl;
  throw;
}

char SSAX::BufAttrValueToMatchingQuote(std::ifstream &ifs, char quot) {
  char next;
  while (EOF != (next = ifs.get())) {
    if (next == quot) {
      return next;
    }
    AddToBuf(next);
  }
  std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
            << std::endl;
  throw;
}

std::pair<char const *, char const *>
SSAX::HandleAttribute(std::ifstream &ifs) {
  return std::pair<char const *, char const *>(nullptr, nullptr);
}

void SSAX::HandleOpenTag(std::ifstream &ifs) {
  char next = ifs.get();
  size_t AttrNameStart = 0;
  size_t AttrValueStart = 0;
  attrs.clear();
  bool IsInteresting = true;

  if (next == EOF) {
    std::cerr << "[SSAX - ERROR]: encountered EOF directly after element "
                 "opening at depth: "
              << depth << std::endl;
    throw;
  }

  if (!IsNameStartChar(next)) {
    std::cerr << "[SSAX - ERROR]: Bad element name first character: " << next
              << std::endl;
    throw;
  }
  AddToBuf(next);

  next = BufNameToBoundary(ifs);
  NullTermBuf();
  AttrNameStart = buf_usage;

  if (!EventOnAllTags) {
    IsInteresting = false;
    for (size_t it_it = 0; it_it < InterestingTags.size(); ++it_it) {
      if (InterestingTags[it_it] == buf) {
        IsInteresting = true;
        break;
      }
    }
  }

  if ((next == '>')) { // Handles elements like <bla>
    if (IsInteresting) {
      OnStartElement(buf, attrs);
    }
    depth++;
    return;
  }
  if ((next == '/')) { // Handles elements like <bla/>
    next = ifs.get();
    if (EOF == next) {
      std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
                << std::endl;
      throw;
    }

    if ((next == '>')) {
      if (IsInteresting) {
        OnStartElement(buf, attrs);
      }
      return;
    } else {
      std::cerr << "[SSAX - ERROR]: encountered unexpected character " << next
                << " after /. Expected >. Exiting early at depth: " << depth
                << std::endl;
      throw;
    }
  }

  while (EOF != (next = ifs.get())) {
    if (isspace(next)) {
      next = IgnoreSpace(ifs);
    }

    if (IsNameStartChar(next)) {
      AddToBuf(next);
      next = BufNameToBoundary(ifs);
      if (next != '=') {
        std::cerr << "[SSAX - ERROR]: encountered unexpected character \""
                  << next
                  << "\". Expected \"=\". Exiting early at depth: " << depth
                  << std::endl;
        throw;
      }
      NullTermBuf();
      AttrValueStart = buf_usage;

      next = ifs.get();
      if ((next != '\"') && (next != '\'')) {
        std::cerr << "[SSAX - ERROR]: encountered unexpected character \""
                  << next << "\" (" << int(next) << "). Expected \"\'\"("
                  << int('\'') << ") or \"\"\"(" << int('\"')
                  << "). Exiting early at depth: " << depth << std::endl;
        throw;
      }
      next = BufAttrValueToMatchingQuote(ifs, next);
      NullTermBuf();
      attrs.push_back(
          std::make_pair(&buf[AttrNameStart], &buf[AttrValueStart]));
      AttrNameStart = buf_usage;
      continue;
    }

    if ((next == '>')) { // Handles elements like <bla>
      if (IsInteresting) {
        OnStartElement(buf, attrs);
      }
      depth++;
      return;
    }
    if ((next == '/')) { // Handles elements like <bla/>
      next = ifs.get();
      if (EOF == next) {
        std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
                  << std::endl;
        throw;
      }

      if ((next == '>')) {
        if (IsInteresting) {
          OnStartElement(buf, attrs);
        }
        // No depth ++;
        return;
      } else {
        std::cerr << "[SSAX - ERROR]: encountered unexpected character " << next
                  << " after /. Expected >. Exiting early at depth: " << depth
                  << std::endl;
        throw;
      }
    }

    std::cerr << "[SSAX - ERROR]: encountered unexpected character " << next
              << " after /. Expected >. Exiting early at depth: " << depth
              << std::endl;
    throw;
  }
  std::cerr << "[SSAX - ERROR]: encountered EOF early at depth: " << depth
            << std::endl;
  throw;
}

void SSAX::HandleCloseTag(std::ifstream &ifs) {
  ifs.ignore();
  char next = ifs.get();
  bool IsInteresting = true;

  if (!IsNameStartChar(next)) {
    std::cerr << "[SSAX - ERROR]: Bad element name first character: " << next
              << std::endl;
    throw;
  }
  AddToBuf(next);

  next = BufNameToBoundary(ifs);
  NullTermBuf();

  if (!EventOnAllTags) {
    IsInteresting = false;
    for (size_t it_it = 0; it_it < InterestingTags.size(); ++it_it) {
      if (InterestingTags[it_it] == buf) {
        IsInteresting = true;
        break;
      }
    }
  }

  if (IsInteresting) {
    OnEndElement(buf);
  }
}

void SSAX::HandleCommentTag(std::ifstream &ifs) {
  IgnoreAllToRightBracket(ifs);
}

void SSAX::FlushBufferToCharacterHandler() {
  NullTermBuf();
  OnCharacters(buf);
}
void SSAX::ClearBuf() { buf_usage = 0; }

void SSAX::ParseFile(const char *filename) {
  std::ifstream ifs(filename);

  if (!ifs.good()) {
    std::cerr << "[SSAX]: Failed to open " << filename << " for reading."
              << std::endl;
    throw;
  }
  ifs.clear();
  ifs.seekg(0);

  OnStartDocument();

  char next;
  while (EOF != (next = ifs.get())) {
    if (LeavingAlready) {
      return;
    }
    if (next == '<') {
      FlushBufferToCharacterHandler();
      ClearBuf();
      switch (ifs.peek()) {
      case '!':
      case '?': {
        HandleCommentTag(ifs);
        ClearBuf();
        break;
      }
      case '/': {
        HandleCloseTag(ifs);
        ClearBuf();
        break;
      }
      default: {
        HandleOpenTag(ifs);
        ClearBuf();
      }
      }
      continue;
    }
    AddToBuf(next);
  }
  OnEndDocument();
}

GENIESplineGetter::GENIESplineGetter(std::vector<std::string> sn) : SSAX() {
  SplineNames = sn;

  DoneReadingSplineElements = false;
  NSplinesToRead = SplineNames.size();
  InSplineElement = NSplinesToRead;
  NSplinesRead = 0;

  ReadingDocument = false;

  InEElement = false;
  InxsecElement = false;

  for (size_t sn_it = 0; sn_it < NSplinesToRead; ++sn_it) {
    SplinesRead.push_back(false);
    EValues.push_back(std::vector<double>());
    XSecValues.push_back(std::vector<double>());
  }
}

void GENIESplineGetter::OnStartDocument() {
  ReadingDocument = true;
  std::cout << "[GENIESplineGetter]: Begin document..." << std::endl;
}
void GENIESplineGetter::OnEndDocument() {
  ReadingDocument = false;
  std::cout << "[GENIESplineGetter]: Finished document." << std::endl;
  if (!(NSplinesRead == NSplinesToRead)) {
    std::cout << "[GENIESplineGetter]: Error only found " << NSplinesRead << "/"
              << NSplinesToRead << " splines:" << std::endl;
    for(size_t s_it = 0; s_it < SplineNames.size(); ++s_it){
      if(!SplinesRead[s_it]){
        std::cout << "\tFailed to read \"" << SplineNames[s_it] << "\"." << std::endl;
      }
    }
    throw;
  }
}

void GENIESplineGetter::OnStartElement(
    const char *name,
    std::vector<std::pair<const char *, const char *>> &attributes) {
  std::string name_str = std::string(name);

  if (InSplineElement != NSplinesToRead) {
    if (name_str == "E") {
      InEElement = true;
    }
    if (name_str == "xsec") {
      InxsecElement = true;
    }
  }

  if (name_str == "spline") {
    for (size_t attr_it = 0; attr_it < attributes.size(); ++attr_it) {
      if (std::string(attributes[attr_it].first) == "name") {
        std::string attr_val = std::string(attributes[attr_it].second);
        for (size_t sn_it = 0; sn_it < NSplinesToRead; ++sn_it) {
          if (SplinesRead[sn_it]) {
            continue;
          }
          if (attr_val == SplineNames[sn_it]) {
            InSplineElement = sn_it;
            return;
          }
        }
      }
    }
  }
}

void GENIESplineGetter::OnCharacters(const char *characters) {
  if (InEElement) {
    std::string ch(characters);
    chomp(ch);
    EValues[InSplineElement].push_back(str2T<double>(ch));
    InEElement = false;
  }
  if (InxsecElement) {
    std::string ch(characters);
    chomp(ch);
    XSecValues[InSplineElement].push_back(str2T<double>(ch));
    InxsecElement = false;
  }
}

void GENIESplineGetter::OnEndElement(const char *name) {
  if ((InSplineElement != NSplinesToRead) && std::string(name) == "spline") {
    std::cout << "[GENIESplineGetter]: Finished reading spline element: "
              << SplineNames[InSplineElement] << std::endl;

    SplinesRead[InSplineElement] = true;
    InSplineElement = NSplinesToRead;
    NSplinesRead++;
    if ((NSplinesRead == NSplinesToRead)) {
      std::cout << "[GENIESplineGetter]: Leaving after reading "
                << NSplinesToRead << " splines." << std::endl;
      NotifyEarlyExit();
    }
  }
}

void GENIESplineGetter::OnComment(const char *) {}

std::vector<TGraph> GENIESplineGetter::GetTGraphs() {
  std::vector<TGraph> graphs;
  for (size_t sp_it = 0; sp_it < NSplinesToRead; ++sp_it) {
    if (!SplinesRead[sp_it]) {
      continue;
    }

    graphs.emplace_back(EValues[sp_it].size(), EValues[sp_it].data(),
                        XSecValues[sp_it].data());
    graphs.back().SetName(SplineNames[sp_it].c_str());
  }
  return graphs;
}

GENIEXSecReader::~GENIEXSecReader() { delete saxParser; }
GENIEXSecReader::GENIEXSecReader(
    std::vector<std::pair<std::string, std::vector<std::string>>> const
        &components) {
  XSecComponents = components;
  std::vector<std::string> splines;

  for (size_t c_it = 0; c_it < XSecComponents.size(); ++c_it) {
    for (size_t s_it = 0; s_it < XSecComponents[c_it].second.size(); ++s_it) {
      splines.push_back(XSecComponents[c_it].second[s_it]);
    }
  }

  saxParser = new GENIESplineGetter(splines);
}
GENIEXSecReader::GENIEXSecReader(std::string const &XSecInputsFile)
    : saxParser(nullptr) {
  TXMLEngine xE;
  xE.SetSkipComments(true);
  XMLDocPointer_t doc = xE.ParseFile(XSecInputsFile.c_str());
  if (!doc) {
    std::cout << "[ERROR]: Attempted to parse XML file: " << XSecInputsFile
              << ", but failed." << std::endl;
    exit(1);
  }
  XMLNodePointer_t rootNode = xE.DocGetRootElement(doc);

  XMLNodePointer_t root_child = xE.GetChild(rootNode);
  while (root_child) { // Look for category node
    if (std::string("category") == xE.GetNodeName(root_child)) {
      bool found = false;
      std::string categ_name =
          GetXMLAttributeValue<std::string>(root_child, "name", found);
      if (!found) {
        std::cout << "[ERROR]: Found unnamed xsec category." << std::endl;
        throw;
      }

      XSecComponents.push_back(
          std::make_pair(categ_name, std::vector<std::string>()));

      std::cout << "[GENIEXSecReader]: Building category: " << categ_name
                << std::endl;

      XMLNodePointer_t categ_child = xE.GetChild(root_child);

      while (categ_child) { // Look for category node
        if (std::string("spline") == xE.GetNodeName(categ_child)) {
          found = false;
          std::string spline_name =
              GetXMLAttributeValue<std::string>(categ_child, "name", found);
          if (!found) {
            std::cout << "[ERROR]: Found spline tag with no name." << std::endl;
            throw;
          }

          XSecComponents.back().second.push_back(spline_name);
          std::cout << "[GENIEXSecReader]: \tReading GENIE spline: "
                    << spline_name << std::endl;
        }

        categ_child = xE.GetNext(categ_child);
      }
    }

    root_child = xE.GetNext(root_child);
  }

  std::vector<std::string> splines;

  for (size_t c_it = 0; c_it < XSecComponents.size(); ++c_it) {
    for (size_t s_it = 0; s_it < XSecComponents[c_it].second.size(); ++s_it) {
      splines.push_back(XSecComponents[c_it].second[s_it]);
    }
  }

  saxParser = new GENIESplineGetter(splines);
}

void GENIEXSecReader::Read(std::string const &gxmlfile) {
  saxParser->ParseFile(gxmlfile.c_str());
}

std::vector<std::pair<std::string, TGraph>>
GENIEXSecReader::GetTGraphs(double min, double max) {
  std::vector<TGraph> splines = saxParser->GetTGraphs();

  std::vector<std::pair<std::string, TGraph>> rtn;

  size_t tc = 0;
  for (size_t c_it = 0; c_it < XSecComponents.size(); ++c_it) {
    size_t Np = 1000;

    TGraph Sum(1000);

    double step = (max - min) / double(Np);

    for (Int_t p_it = 0; p_it < Sum.GetN(); ++p_it) {
      Sum.SetPoint(p_it, min + p_it * step, 0);
    }

    for (size_t s_it = 0; s_it < XSecComponents[c_it].second.size(); ++s_it) {
      for (Int_t p_it = 0; p_it < Sum.GetN(); ++p_it) {
        double x1, y1, ev;
        Sum.GetPoint(p_it, x1, y1);

        //Copied from GENIE/Convents/Units.h
        static const double gigaelectronvolt  = 1.;
        static const double GeV               = gigaelectronvolt;
        static const double meter    = 5.07e+15 / GeV;
        static const double centimeter  = 0.01*meter;
        static const double centimeter2 = centimeter*centimeter;

        ev = splines[tc].Eval(x1) / centimeter2;

        Sum.SetPoint(p_it, x1, y1 + ev);
      }

      tc++;
    }

    rtn.push_back(std::make_pair(XSecComponents[c_it].first, Sum));
  }

  return rtn;
}
