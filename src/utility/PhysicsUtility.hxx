#ifndef PHYSICSUTILTIY_HXX_SEEN
#define PHYSICSUTILTIY_HXX_SEEN

#include "InteractionModel.hxx"

#include <string>

std::string GetSpeciesName(int pdg);

struct GENIECodeStringParser {
  int nu_PDG;
  InteractionModel::TrueChannel channel;
  bool IsCC;

  GENIECodeStringParser(std::string const &evc);
};

#endif
