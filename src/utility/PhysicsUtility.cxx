#include "PhysicsUtility.hxx"

#include "StringParserUtility.hxx"

#include <vector>

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

GENIECodeStringParser::GENIECodeStringParser(std::string const &evc) {
  std::vector<std::string> split_evc = ParseToVect<std::string>(evc, ",");

  std::vector<std::string> split_for_nu =
      ParseToVect<std::string>(split_evc.front(), ";");
  std::vector<std::string> split_for_nu_pdg =
      ParseToVect<std::string>(split_for_nu.front(), ":");

  nu_PDG = str2T<int>(split_for_nu_pdg[1]);

  if (evc.find("MEC") != std::string::npos) {
    channel = InteractionModel::TrueChannel::k2p2h;
  } else if (evc.find("QES") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kQES;
  } else if (evc.find("RES") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kRES;
  } else if (evc.find("DIS") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kDIS;
  } else if (evc.find("COH") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kCOH;
  } else if (evc.find("NuEEL") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kNuEEL;
  } else if (evc.find("IMD") != std::string::npos) {
    channel = InteractionModel::TrueChannel::kIMD;
  } else {
    std::cout << "[ERROR]: Unaccounted for channel string in: " << evc
              << std::endl;
    throw;
  }

  if (evc.find("[CC]") != std::string::npos) {
    IsCC = true;
  } else if (evc.find("[NC]") != std::string::npos) {
    IsCC = false;
  } else {
    std::cout << "[WARN]: Couldn't find CC/NC in ev code: " << evc
              << std::endl;
    IsCC = false;
  }
}
