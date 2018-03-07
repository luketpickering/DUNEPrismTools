#include "Utils.hxx"

#include <iostream>
#include <utility>
#include <vector>

namespace VALORModel {
enum class TrueClass;
enum class TrueChannel;
}

#define VARLIST          \
  X(kNu_QE_1, 0.082)     \
  X(kNu_QE_2, 0.23)      \
  X(kNu_QE_3, 0.48)      \
  X(kNuBar_QE_1, 0.087)  \
  X(kNuBar_QE_2, 0.24)   \
  X(kNuBar_QE_3, 0.40)   \
  X(kNu_MEC, 1)          \
  X(kNuBar_MEC, 1)       \
  X(kNu_1Pi0_1, 0.13)    \
  X(kNu_1Pi0_2, 0.23)    \
  X(kNu_1Pi0_3, 0.35)    \
  X(kNu_1PiC_1, 0.13)    \
  X(kNu_1PiC_2, 0.24)    \
  X(kNu_1PiC_3, 0.40)    \
  X(kNuBar_1Pi0_1, 0.16) \
  X(kNuBar_1Pi0_2, 0.27) \
  X(kNuBar_1Pi0_3, 0.35) \
  X(kNuBar_1PiC_1, 0.16) \
  X(kNuBar_1PiC_2, 0.3)  \
  X(kNuBar_1PiC_3, 0.4)  \
  X(kNu_2Pi, 0.22)       \
  X(kNuBar_2Pi, 0.22)    \
  X(kNu_DIS_1, 0.035)    \
  X(kNu_DIS_2, 0.035)    \
  X(kNu_DIS_3, 0.027)    \
  X(kNuBar_DIS_1, 0.01)  \
  X(kNuBar_DIS_2, 0.017) \
  X(kNuBar_DIS_3, 0.017) \
  X(kNu_COH, 1.28)       \
  X(kNuBar_COH, 1.34)    \
  X(kNu_NC, 0.16)        \
  X(kNuBar_NC, 0.16)     \
  X(kNuMu_E_Ratio, 0.03) \
  X(kNVALORDials,0)

#define X(A, B) A,
enum class VALORModel::TrueClass { VARLIST };
#undef X
#define X(A, B)                    \
  case VALORModel::TrueClass::A: { \
    return os << #A;               \
  }

inline std::ostream &operator<<(std::ostream &os, VALORModel::TrueClass te) {
  switch (te) { VARLIST }
  return os;
}
#undef X
#define X(A, B)                    \
  case VALORModel::TrueClass::A: { \
    return B;                      \
  }

inline double GetDefaultDialTweak(VALORModel::TrueClass te) {
  switch (te) { VARLIST }
  return 0;
}
#undef X
#define X(A, B) dv.push_back(std::make_tuple(VALORModel::TrueClass::A, 1, B));
inline std::vector<std::tuple<VALORModel::TrueClass, double, double> >
GetDialValueVector(VALORModel::TrueClass te) {
  std::vector<std::tuple<VALORModel::TrueClass, double, double> > dv;
  VARLIST
  return dv;
}
#undef X
#undef VARLIST

#define VARLIST \
  Y(kQES, 1)    \
  X(kMEC)       \
  X(kRES)       \
  X(kDIS)       \
  X(kCOH)       \
  X(kNuEEL)     \
  X(kIMD)       \
  X(kUnknown)

#define X(A) A,
#define Y(A, B) A = B,
enum class VALORModel::TrueChannel { VARLIST };
#undef X
#undef Y
#define X(A)                         \
  case VALORModel::TrueChannel::A: { \
    return os << #A;                 \
  }
#define Y(A, B)                      \
  case VALORModel::TrueChannel::A: { \
    return os << #A;                 \
  }
inline std::ostream &operator<<(std::ostream &os, VALORModel::TrueChannel te) {
  switch (te) { VARLIST }
  return os;
}
#undef X
#undef Y
#undef VARLIST

struct GENIECodeStringParser {
  int nu_PDG;
  VALORModel::TrueChannel channel;
  bool IsCC;

  GENIECodeStringParser(std::string const &evc) {
    std::vector<std::string> split_evc = ParseToVect<std::string>(evc, ",");

    std::vector<std::string> split_for_nu =
        ParseToVect<std::string>(split_evc.front(), ";");
    std::vector<std::string> split_for_nu_pdg =
        ParseToVect<std::string>(split_for_nu.front(), ":");

    nu_PDG = str2T<int>(split_for_nu_pdg[1]);

    if (evc.find("MEC") != std::string::npos) {
      channel = VALORModel::TrueChannel::kMEC;
    } else if (evc.find("QES") != std::string::npos) {
      channel = VALORModel::TrueChannel::kQES;
    } else if (evc.find("RES") != std::string::npos) {
      channel = VALORModel::TrueChannel::kRES;
    } else if (evc.find("DIS") != std::string::npos) {
      channel = VALORModel::TrueChannel::kDIS;
    } else if (evc.find("COH") != std::string::npos) {
      channel = VALORModel::TrueChannel::kCOH;
    } else if (evc.find("NuEEL") != std::string::npos) {
      channel = VALORModel::TrueChannel::kNuEEL;
    } else if (evc.find("IMD") != std::string::npos) {
      channel = VALORModel::TrueChannel::kIMD;
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
      std::cout << "[ERROR]: Couldn't find CC/NC in ev code: " << evc
                << std::endl;
      throw;
    }
  }
};

std::vector<VALORModel::TrueClass> GetApplicableDials(EDep const &ed) {
  std::vector<VALORModel::TrueClass> AppDials;

  AppDials.push_back(VALORModel::TrueClass::kNuMu_E_Ratio);

#if DEBUG
  GENIECodeStringParser gcp(ed.EventCode->GetString().Data());

  if (gcp.nu_PDG != ed.nu_PDG) {
    std::cout << "[ERROR]: Failed GENIE passthrough: nu-PDG from G4Ar = "
              << ed.nu_PDG << " from GENIE event string = " << gcp.nu_PDG
              << std::endl;
    throw;
  }

  if (gcp.IsCC != ed.IsCC) {
    std::cout << "[ERROR]: Failed GENIE passthrough: IsCC from G4Ar = "
              << ed.IsCC << "(nu: " << ed.nu_PDG
              << ", lep: " << ed.PrimaryLepPDG
              << ") from GENIE event string = " << gcp.IsCC << std::endl;
    throw;
  }

  if (ed.IsAntinu != (gcp.nu_PDG < 0)) {
    std::cout << "[ERROR]: Failed GENIE passthrough: IsAntiNu from G4Ar = "
              << ed.IsAntinu << " from GENIE event string = " << gcp.nu_PDG
              << std::endl;
    throw;
  }
#endif

  // Not good enough to get post FSI NPi. Need to work out from channel.

  // Q2_True
  // IsAntinu
  switch (static_cast<VALORModel::TrueChannel>(ed.GENIEInteractionTopology)) {
    case VALORModel::TrueChannel::kQES: {
      if (ed.Q2_True < 0.2) {
        AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_QE_1
                                       : VALORModel::TrueClass::kNu_QE_1);
      } else if (ed.Q2_True < 0.55) {
        AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_QE_2
                                       : VALORModel::TrueClass::kNu_QE_2);
      } else {
        AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_QE_3
                                       : VALORModel::TrueClass::kNu_QE_3);
      }
      break;
    }
    case VALORModel::TrueChannel::kRES: {
      if (ed.Is1Pi0) {
        if (ed.Q2_True < 0.35) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1Pi0_1
                                         : VALORModel::TrueClass::kNu_1Pi0_1);
        } else if (ed.Q2_True < 0.90) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1Pi0_2
                                         : VALORModel::TrueClass::kNu_1Pi0_2);
        } else {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1Pi0_3
                                         : VALORModel::TrueClass::kNu_1Pi0_3);
        }
      } else {
        if (ed.Q2_True < 0.35) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1PiC_1
                                         : VALORModel::TrueClass::kNu_1PiC_1);
        } else if (ed.Q2_True < 0.80) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1PiC_2
                                         : VALORModel::TrueClass::kNu_1PiC_2);
        } else {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_1PiC_3
                                         : VALORModel::TrueClass::kNu_1PiC_3);
        }
      }
      break;
    }
    case VALORModel::TrueChannel::kDIS: {
      if (((ed.NPi0 + ed.NPiC) == 2) &&
          ((ed.NGamma + ed.NBaryonicRes + ed.NOther) == 0)) {
        AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_2Pi
                                       : VALORModel::TrueClass::kNu_2Pi);
      } else {
        if (ed.nu_4mom[3] < 7.5) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_DIS_1
                                         : VALORModel::TrueClass::kNu_DIS_1);
        } else if (ed.nu_4mom[3] < 15) {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_DIS_2
                                         : VALORModel::TrueClass::kNu_DIS_2);
        } else {
          AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_DIS_3
                                         : VALORModel::TrueClass::kNu_DIS_3);
        }
      }
      break;
    }
    case VALORModel::TrueChannel::kCOH: {
      AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_COH
                                     : VALORModel::TrueClass::kNu_COH);
      break;
    }
    case VALORModel::TrueChannel::kNuEEL: {
      break;
    }
    case VALORModel::TrueChannel::kIMD: {
      break;
    }
    default: { throw; }
  }

  if (!ed.IsCC) {
    AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_NC
                                   : VALORModel::TrueClass::kNu_NC);
  }

  return AppDials;
}

double GetVALORWeight(VALORModel::TrueClass te, double value, EDep const &ed) {
  switch (te) {
    case VALORModel::TrueClass::kNu_QE_1:
    case VALORModel::TrueClass::kNu_QE_2:
    case VALORModel::TrueClass::kNu_QE_3:
    case VALORModel::TrueClass::kNuBar_QE_1:
    case VALORModel::TrueClass::kNuBar_QE_2:
    case VALORModel::TrueClass::kNuBar_QE_3:
    case VALORModel::TrueClass::kNu_MEC:
    case VALORModel::TrueClass::kNuBar_MEC:
    case VALORModel::TrueClass::kNu_1Pi0_1:
    case VALORModel::TrueClass::kNu_1Pi0_2:
    case VALORModel::TrueClass::kNu_1Pi0_3:
    case VALORModel::TrueClass::kNu_1PiC_1:
    case VALORModel::TrueClass::kNu_1PiC_2:
    case VALORModel::TrueClass::kNu_1PiC_3:
    case VALORModel::TrueClass::kNuBar_1Pi0_1:
    case VALORModel::TrueClass::kNuBar_1Pi0_2:
    case VALORModel::TrueClass::kNuBar_1Pi0_3:
    case VALORModel::TrueClass::kNuBar_1PiC_1:
    case VALORModel::TrueClass::kNuBar_1PiC_2:
    case VALORModel::TrueClass::kNuBar_1PiC_3:
    case VALORModel::TrueClass::kNu_2Pi:
    case VALORModel::TrueClass::kNuBar_2Pi:
    case VALORModel::TrueClass::kNu_DIS_1:
    case VALORModel::TrueClass::kNu_DIS_2:
    case VALORModel::TrueClass::kNu_DIS_3:
    case VALORModel::TrueClass::kNuBar_DIS_1:
    case VALORModel::TrueClass::kNuBar_DIS_2:
    case VALORModel::TrueClass::kNuBar_DIS_3:
    case VALORModel::TrueClass::kNu_COH:
    case VALORModel::TrueClass::kNuBar_COH:
    case VALORModel::TrueClass::kNu_NC:
    case VALORModel::TrueClass::kNuBar_NC: {
      return 1.0 + value;
    }
    // Splits the uncertainty between numu and nue normalisations that go
    // opposite ways.
    case VALORModel::TrueClass::kNuMu_E_Ratio: {
      return ed.IsNumu ? 1.0 + 0.5*value
                       : 1.0 - 0.5*value;
    }
    default: { throw; }
  }
}
