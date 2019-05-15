#ifndef GENIESPLINEREADER_HXX_SEEN
#define GENIESPLINEREADER_HXX_SEEN

#include "TGraph.h"

#include <fstream>
#include <string>
#include <utility>
#include <vector>

class SSAX {
  bool EventOnAllTags;
  bool LeavingAlready;
  std::vector<std::string> InterestingTags;

  size_t buf_size;
  size_t buf_usage;
  char *buf;
  std::vector<std::pair<const char *, const char *>> attrs;

  size_t depth;

public:
  SSAX();

  virtual ~SSAX();

  virtual void OnStartDocument();
  virtual void OnEndDocument();
  virtual void OnStartElement(
      const char *name,
      std::vector<std::pair<const char *, const char *>> &attributes);
  virtual void OnCharacters(const char *characters);
  virtual void OnEndElement(const char *name);
  virtual void OnComment(const char *);

  void NotifyInterestingTags(std::vector<std::string> it);
  void NotifyEarlyExit();

private:
  void ExpandBuf();

  void AddToBuf(char character);
  void NullTermBuf();

  bool IsNameStartChar(char c);
  bool IsNameChar(char c);

  void IgnoreAllToRightBracket(std::ifstream &ifs);

  char IgnoreSpace(std::ifstream &ifs);

  char BufNameToBoundary(std::ifstream &ifs);

  char BufAttrValueToMatchingQuote(std::ifstream &ifs, char quot);

  std::pair<char const *, char const *> HandleAttribute(std::ifstream &ifs);

  void HandleOpenTag(std::ifstream &ifs);
  void HandleCloseTag(std::ifstream &ifs);
  void HandleCommentTag(std::ifstream &ifs);

  void FlushBufferToCharacterHandler();

  void ClearBuf();

public:
  void ParseFile(const char *filename);
};

class GENIESplineGetter : public SSAX {
public:
  std::vector<std::vector<double>> EValues;
  std::vector<std::vector<double>> XSecValues;

  std::vector<std::string> SplineNames;
  std::vector<bool> SplinesRead;
  size_t NSplinesRead;

  bool DoneReadingSplineElements;

  size_t NSplinesToRead;
  size_t InSplineElement;

  bool InEElement;
  bool InxsecElement;

  bool ReadingDocument;

  GENIESplineGetter(std::vector<std::string> sn);

  void OnStartDocument();
  void OnEndDocument();
  void OnStartElement(
      const char *name,
      std::vector<std::pair<const char *, const char *>> &attributes);
  void OnCharacters(const char *characters);
  void OnEndElement(const char *name);
  void OnComment(const char *);

  std::vector<TGraph> GetTGraphs();
};

struct GENIEXSecReader {
  std::vector<std::pair<std::string, std::vector<std::string>>> XSecComponents;

  GENIESplineGetter *saxParser;

  ~GENIEXSecReader();

  GENIEXSecReader(std::string const &XSecInputsFile);
  GENIEXSecReader(std::vector<std::pair<std::string, std::vector<std::string>>> const &components);

  void Read(std::string const &gxmlfile);

  std::vector<std::pair<std::string, TGraph>> GetTGraphs(double min,
                                                         double max);
};

#endif
