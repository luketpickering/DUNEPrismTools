#include "PhysicsUtility.hxx"

std::string GetSpeciesName(int pdg) {
  switch (pdg) {
    case -14:
      return "numubar";
    case -12:
      return "nuebar";
    case 12:
      return "nue";
    case 14:
      return "numu";
    default:
      return "unknown";
  }
}
