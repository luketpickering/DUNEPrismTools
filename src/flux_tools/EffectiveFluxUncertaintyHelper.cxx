#include "EffectiveFluxUncertaintyHelper.hxx"

#include "ROOTUtility.hxx"

void EffectiveFluxUncertaintyHelper::Initialize(fhicl::ParameterSet const &ps) {

  NDIs2D = -1;
  NDTweaks.clear();
  FDTweaks.clear();

  bool IsCAFAnaFormat = ps.get<bool>("IsCAFAnaFormat", false);

  size_t NParams = ps.get<size_t>("NEffectiveParametersToRead");
  std::string ND_detector_tag = ps.get<std::string>("ND_detector_tag", "");
  std::string FD_detector_tag = ps.get<std::string>("FD_detector_tag", "");
  std::string nu_mode_beam_tag = ps.get<std::string>("nu_mode_beam_tag", "");
  std::string nubar_mode_beam_tag =
      ps.get<std::string>("nubar_mode_beam_tag", "");
  std::string numu_species_tag = ps.get<std::string>("numu_species_tag", "");
  std::string nue_species_tag = ps.get<std::string>("nue_species_tag", "");
  std::string numubar_species_tag =
      ps.get<std::string>("numubar_species_tag", "");
  std::string nuebar_species_tag =
      ps.get<std::string>("nuebar_species_tag", "");

  std::string input_file = ps.get<std::string>("InputFile");
  std::string input_dir = ps.get<std::string>(
      "InputDir", IsCAFAnaFormat ? "" : "EffectiveFluxParameters");

  for (size_t p_it = 0; p_it < NParams; ++p_it) {
    NDTweaks.emplace_back();
    FDTweaks.emplace_back();
    std::string input_dir_i = input_dir + (input_dir.size() ? "/" : "") +
                              (IsCAFAnaFormat ? "syst" : "param_") +
                              std::to_string(p_it) + "/";

    size_t NHistsLoaded = 0;
    int nucf = kND_numu_numode;
    for (std::string const &location_tag : {ND_detector_tag, FD_detector_tag}) {
      for (std::string const &beam_mode_tag :
           {nu_mode_beam_tag, nubar_mode_beam_tag}) {
        for (std::string const &species_tag :
             {numu_species_tag, nue_species_tag, numubar_species_tag,
              nuebar_species_tag}) {

          NDTweaks.back().emplace_back(nullptr);
          FDTweaks.back().emplace_back(nullptr);

          std::string hname =
              (IsCAFAnaFormat ? input_dir_i + location_tag + "_" + species_tag +
                                    "_" + beam_mode_tag
                              : input_dir_i + location_tag + "_" +
                                    beam_mode_tag + "_" + species_tag);

          if (nucf < kFD_numu_numode) { // Is ND
            NDTweaks.back().back() =
                GetHistogram_uptr<TH1>(input_file, hname, true);

            if ((NDIs2D == -1) && NDTweaks.back().back()) {
              NDIs2D = NDTweaks.back().back()->GetDimension() - 1;
            }

            NHistsLoaded += bool(NDTweaks.back().back());
            // if (!NDTweaks.back().back()) {
            //   std::cout << "failed to read " << hname << std::endl;
            // }
          } else {
            FDTweaks.back().back() =
                GetHistogram_uptr<TH1>(input_file, hname, true);

            NHistsLoaded += bool(FDTweaks.back().back());
            // if (!FDTweaks.back().back()) {
            //   std::cout << "failed to read " << hname << std::endl;
            // }
          }
          nucf += 1;
        }
      }
    }
    // std::cout << "[EffectiveFluxParameters]: Loaded " << NHistsLoaded
    //           << " tweak histograms for parameter " << p_it << std::endl;

    if (!NHistsLoaded) {
      NDTweaks.pop_back();
      FDTweaks.pop_back();
    }
  }
}
int EffectiveFluxUncertaintyHelper::GetNuConfig(int nu_pdg, bool IsND,
                                                bool IsNuMode) {

  int nucf;

  switch (nu_pdg) {
  case 14: {
    nucf = IsNuMode ? kND_numu_numode : kND_numu_nubarmode;
    break;
  }
  case -14: {
    nucf = IsNuMode ? kND_numubar_numode : kND_numubar_nubarmode;
    break;
  }
  case 12: {
    nucf = IsNuMode ? kND_nue_numode : kND_nue_nubarmode;
    break;
  }
  case -12: {
    nucf = IsNuMode ? kND_nuebar_numode : kND_nuebar_nubarmode;
    break;
  }
  }

  return nucf + (IsND ? 0 : 8);
}

int EffectiveFluxUncertaintyHelper::GetNuConfig(int nu_pdg, double enu_GeV,
                                                double off_axix_pos_m,
                                                size_t param_id, bool IsND,
                                                bool IsNuMode) {
  if (GetBin(nu_pdg, enu_GeV, off_axix_pos_m, param_id, IsND, IsNuMode) ==
      kInvalidBin) {
    return kUnhandled;
  }
  return GetNuConfig(nu_pdg, enu_GeV, off_axix_pos_m);
}

int EffectiveFluxUncertaintyHelper::GetBin(int nu_pdg, double enu_GeV,
                                           double off_axix_pos_m,
                                           size_t param_id, bool IsND,
                                           bool IsNuMode) {
  int nucf = GetNuConfig(nu_pdg, IsND, IsNuMode);

  if (nucf < kFD_numu_numode) { // Is ND
    if (!NDTweaks[param_id][nucf]) {
      return kInvalidBin;
    } else {
      TH2JaggedF *jag;
      if ((jag = dynamic_cast<TH2JaggedF *>(NDTweaks[param_id][nucf].get()))) {
        Int_t gbin = jag->FindFixBin(enu_GeV, off_axix_pos_m);
        return jag->IsFlowBin(gbin) ? kInvalidBin : gbin;
      } else {
        Int_t xbin = NDTweaks[param_id][nucf]->GetXaxis()->FindFixBin(enu_GeV);
        if ((xbin == 0) ||
            (xbin == NDTweaks[param_id][nucf]->GetXaxis()->GetNbins() + 1)) {
          return kInvalidBin;
        }
        if (NDIs2D == 1) {
          Int_t ybin =
              NDTweaks[param_id][nucf]->GetYaxis()->FindFixBin(off_axix_pos_m);
          if ((ybin == 0) ||
              (ybin == NDTweaks[param_id][nucf]->GetYaxis()->GetNbins() + 1)) {
            return kInvalidBin;
          }
          return NDTweaks[param_id][nucf]->GetBin(xbin, ybin);
        } else {
          return NDTweaks[param_id][nucf]->GetBin(xbin);
        }
      }
    }
  } else {
    if (!FDTweaks[param_id][nucf]) {
      return kInvalidBin;
    } else {
      Int_t bin = FDTweaks[param_id][nucf]->FindFixBin(enu_GeV);
      if ((bin == 0) ||
          (bin == FDTweaks[param_id][nucf]->GetXaxis()->GetNbins() + 1)) {
        return kInvalidBin;
      }
      return bin;
    }
  }
}

double EffectiveFluxUncertaintyHelper::GetFluxWeight(size_t param_id,
                                                     double param_val, int bin,
                                                     int nucf) const {
  if (nucf == kUnhandled) {
    return 1;
  }

  if (bin == kInvalidBin) {
    return 1;
  }

  if (nucf < kFD_numu_numode) {
    return 1 + param_val * (NDTweaks[param_id][nucf]->GetBinContent(bin));
  } else {
    return 1 + param_val * (FDTweaks[param_id][nucf]->GetBinContent(bin));
  }
}
