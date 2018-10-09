#include "VALORModelClassifier.hxx"

#include "PhysicsUtility.hxx"

std::vector<VALORModel::TrueClass> GetApplicableDials(DepositsSummary const &ed) {
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
  switch (static_cast<InteractionModel::TrueChannel>(ed.GENIEInteractionTopology)) {
    case InteractionModel::TrueChannel::kQES: {
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
    case InteractionModel::TrueChannel::k2p2h:{
      AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_MEC
                                       : VALORModel::TrueClass::kNu_MEC);
      break;
    }
    case InteractionModel::TrueChannel::kRES: {
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
    case InteractionModel::TrueChannel::kDIS: {
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
    case InteractionModel::TrueChannel::kCOH: {
      AppDials.push_back(ed.IsAntinu ? VALORModel::TrueClass::kNuBar_COH
                                     : VALORModel::TrueClass::kNu_COH);
      break;
    }
    case InteractionModel::TrueChannel::kNuEEL: {
      break;
    }
    case InteractionModel::TrueChannel::kIMD: {
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

double GetVALORWeight(VALORModel::TrueClass te, double value, DepositsSummary const &ed) {
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
