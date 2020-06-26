#include "EffectiveFluxUncertaintyHelper.hxx"

#include "ROOTUtility.hxx"

void EffectiveFluxUncertaintyHelper::Initialize(fhicl::ParameterSet const &ps) {

  NDIs2D = -1;
  NDSpecHCRunIs2D = -1;
  NDTweaks.clear();
  NDSpecHCRunTweaks.clear();
  FDTweaks.clear();

  bool IsCAFAnaFormat = ps.get<bool>("IsCAFAnaFormat", false);

  std::string ND_detector_tag = ps.get<std::string>("ND_detector_tag", "");
  std::string ND_SpecHCRun_detector_tag =
      ps.get<std::string>("ND_SpecHCRun_detector_tag", "");
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
  std::string input_dir =
      ps.get<std::string>("InputDir", IsCAFAnaFormat ? "" : "FluxParameters");

  TFile *f = TFile::Open(input_file.c_str(), "READ");

  if (!f || !f->IsOpen()) {
    std::cout << "[ERROR]: Failed to open TFile: " << input_file << std::endl;
    abort();
  }

  TDirectory *d = f->GetDirectory(input_dir.c_str());

  if (!d) {
    std::cout << "[ERROR]: Failed to open TDirectory: " << input_dir
              << " in file: " << input_file << std::endl;
    abort();
  }

  std::vector<std::string> param_names =
      ps.get<std::vector<std::string>>("param_names", {});

  // Try loading them from
  if (!param_names.size()) {
    size_t NMaxParams =
        ps.get<size_t>("nmax_params", std::numeric_limits<size_t>::max());

    TList *param_names_list = nullptr;
    d->GetObject("param_names", param_names_list);
    if (!param_names_list) {
      std::cout << "[ERROR]: Failed to find param_name list in TDirectory: "
                << input_dir << " in file: " << input_file
                << " and was passed no param_names FHiCL list to try to load."
                << std::endl;
      abort();
    }
    for (size_t p_it = 0;
         p_it < std::min(NMaxParams, size_t(param_names_list->GetSize()));
         ++p_it) {
      param_names.push_back(
          static_cast<TObjString *>(param_names_list->At(p_it))
              ->String()
              .Data());
    }
  }

  size_t NParams = param_names.size();

  for (size_t p_it = 0; p_it < NParams; ++p_it) {
    NDTweaks.emplace_back();
    NDSpecHCRunTweaks.emplace_back();
    FDTweaks.emplace_back();
    std::string input_dir_i =
        input_dir + (input_dir.size() ? "/" : "") + param_names[p_it] + "/";

    size_t NHistsLoaded = 0;
    int nucf = kND_numu_numode;
    for (std::string const &location_tag :
         {ND_detector_tag, ND_SpecHCRun_detector_tag, FD_detector_tag}) {
      for (std::string const &beam_mode_tag :
           {nu_mode_beam_tag, nubar_mode_beam_tag}) {
        for (std::string const &species_tag :
             {numu_species_tag, nue_species_tag, numubar_species_tag,
              nuebar_species_tag}) {

          NDTweaks.back().emplace_back(nullptr);
          NDSpecHCRunTweaks.back().emplace_back(nullptr);
          FDTweaks.back().emplace_back(nullptr);

          std::string hname =
              (IsCAFAnaFormat ? input_dir_i + location_tag + "_" + species_tag +
                                    "_" + beam_mode_tag
                              : input_dir_i + location_tag + "_" +
                                    beam_mode_tag + "_" + species_tag);

          if (nucf < kND_SpecHCRun_numu_numode) { // Is ND
            NDTweaks.back().back() =
                GetHistogram_uptr<TH1>(input_file, hname, true);
            if (!NDTweaks.back().back()) {
              std::cout << "[WARN]: Failed to load tweak: " << hname << " from "
                        << input_file << std::endl;
            }
            if ((NDIs2D == -1) && NDTweaks.back().back()) {
              NDIs2D = NDTweaks.back().back()->GetDimension() - 1;
            }
            NHistsLoaded += bool(NDTweaks.back().back());
          } else if (nucf < kFD_numu_numode) {
            // hack becasue I have been inconsistent in whether SPecRuns are
            // considered a different 'detector' or a different 'beam mode'
            std::string hname =
                (IsCAFAnaFormat
                     ? input_dir_i + ND_detector_tag + "_" + species_tag + "_" +
                           ND_SpecHCRun_detector_tag + "_" + beam_mode_tag
                     : input_dir_i + ND_detector_tag + "_" + beam_mode_tag +
                           "_" + ND_SpecHCRun_detector_tag + "_" + species_tag);
            NDSpecHCRunTweaks.back().back() =
                GetHistogram_uptr<TH1>(input_file, hname, true);
            if (!NDSpecHCRunTweaks.back().back()) {
              std::cout << "[WARN]: Failed to load tweak: " << hname << " from "
                        << input_file << std::endl;
            }
            if ((NDSpecHCRunIs2D == -1) && NDSpecHCRunTweaks.back().back()) {
              NDSpecHCRunIs2D =
                  NDSpecHCRunTweaks.back().back()->GetDimension() - 1;
            }
            NHistsLoaded += bool(NDSpecHCRunTweaks.back().back());
          } else {
            FDTweaks.back().back() =
                GetHistogram_uptr<TH1>(input_file, hname, true);
            if (!FDTweaks.back().back()) {
              std::cout << "[WARN]: Failed to load tweak: " << hname << " from "
                        << input_file << std::endl;
            }
            NHistsLoaded += bool(FDTweaks.back().back());
          }
          nucf += 1;
        }
      }
    }
    std::cout << "[EffectiveFluxParameters]: Loaded " << NHistsLoaded
              << " tweak histograms for parameter " << param_names[p_it]
              << std::endl;

    if (!NHistsLoaded) {
      NDTweaks.pop_back();
      NDSpecHCRunTweaks.pop_back();
      FDTweaks.pop_back();
    }
  }
}
int EffectiveFluxUncertaintyHelper::GetNuConfig(int nu_pdg, bool IsND,
                                                bool IsNuMode,
                                                bool isSpecHCRun) {

  int nucf;

  switch (nu_pdg) {
  case 14: {
    nucf =
        IsNuMode
            ? (isSpecHCRun ? kND_SpecHCRun_numu_numode : kND_numu_numode)
            : (isSpecHCRun ? kND_SpecHCRun_numu_nubarmode : kND_numu_nubarmode);
    break;
  }
  case -14: {
    nucf = IsNuMode ? (isSpecHCRun ? kND_SpecHCRun_numubar_numode
                                   : kND_numubar_numode)
                    : (isSpecHCRun ? kND_SpecHCRun_numubar_nubarmode
                                   : kND_numubar_nubarmode);
    break;
  }
  case 12: {
    nucf = IsNuMode ? (isSpecHCRun ? kND_SpecHCRun_nue_numode : kND_nue_numode)
                    : (isSpecHCRun ? kND_SpecHCRun_nue_nubarmode
                                   : kND_nue_nubarmode);
    break;
  }
  case -12: {
    nucf = IsNuMode
               ? (isSpecHCRun ? kND_SpecHCRun_nuebar_numode : kND_nuebar_numode)
               : (isSpecHCRun ? kND_SpecHCRun_nuebar_nubarmode
                              : kND_nuebar_nubarmode);
    break;
  }
  }

  return nucf + (IsND ? 0 : 16);
}

int EffectiveFluxUncertaintyHelper::GetNuConfig(int nu_pdg, double enu_GeV,
                                                double off_axis_pos_m,
                                                size_t param_id, bool IsND,
                                                bool IsNuMode,
                                                bool isSpecHCRun) {
  if (GetBin(nu_pdg, enu_GeV, off_axis_pos_m, param_id, IsND, IsNuMode,
             isSpecHCRun) == kInvalidBin) {
    return kUnhandled;
  }
  return GetNuConfig(nu_pdg, enu_GeV, off_axis_pos_m, 0, IsND, IsNuMode,
                     isSpecHCRun);
}

int EffectiveFluxUncertaintyHelper::GetBin(int nu_pdg, double enu_GeV,
                                           double off_axis_pos_m,
                                           size_t param_id, bool IsND,
                                           bool IsNuMode, bool isSpecHCRun) {
  int nucf = GetNuConfig(nu_pdg, IsND, IsNuMode, isSpecHCRun);

  if (nucf < kND_SpecHCRun_numu_numode) { // Is ND
    if (!NDTweaks[param_id][nucf]) {
      return kInvalidBin;
    } else {
      TH2JaggedF *jag;
      if ((jag = dynamic_cast<TH2JaggedF *>(NDTweaks[param_id][nucf].get()))) {
        Int_t gbin = jag->FindFixBin(enu_GeV, off_axis_pos_m);
        return jag->IsFlowBin(gbin) ? kInvalidBin : gbin;
      } else {
        Int_t xbin = NDTweaks[param_id][nucf]->GetXaxis()->FindFixBin(enu_GeV);
        if ((xbin == 0) ||
            (xbin == NDTweaks[param_id][nucf]->GetXaxis()->GetNbins() + 1)) {
          return kInvalidBin;
        }
        if (NDIs2D == 1) {
          Int_t ybin =
              NDTweaks[param_id][nucf]->GetYaxis()->FindFixBin(off_axis_pos_m);
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
  } else if (nucf < kFD_numu_numode) {
    if (!NDSpecHCRunTweaks[param_id][nucf]) {
      return kInvalidBin;
    } else {
      TH2JaggedF *jag;
      if ((jag = dynamic_cast<TH2JaggedF *>(
               NDSpecHCRunTweaks[param_id][nucf].get()))) {
        Int_t gbin = jag->FindFixBin(enu_GeV, off_axis_pos_m);
        return jag->IsFlowBin(gbin) ? kInvalidBin : gbin;
      } else {
        Int_t xbin =
            NDSpecHCRunTweaks[param_id][nucf]->GetXaxis()->FindFixBin(enu_GeV);
        if ((xbin == 0) ||
            (xbin ==
             NDSpecHCRunTweaks[param_id][nucf]->GetXaxis()->GetNbins() + 1)) {
          return kInvalidBin;
        }
        if (NDIs2D == 1) {
          Int_t ybin =
              NDSpecHCRunTweaks[param_id][nucf]->GetYaxis()->FindFixBin(
                  off_axis_pos_m);
          if ((ybin == 0) ||
              (ybin ==
               NDSpecHCRunTweaks[param_id][nucf]->GetYaxis()->GetNbins() + 1)) {
            return kInvalidBin;
          }
          return NDSpecHCRunTweaks[param_id][nucf]->GetBin(xbin, ybin);
        } else {
          return NDSpecHCRunTweaks[param_id][nucf]->GetBin(xbin);
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

  if (nucf < kND_SpecHCRun_numu_numode) {
    return 1 + param_val * (NDTweaks[param_id][nucf]->GetBinContent(bin));
  } else if (nucf < kFD_numu_numode) {
    return 1 +
           param_val * (NDSpecHCRunTweaks[param_id][nucf]->GetBinContent(bin));
  } else {
    return 1 + param_val * (FDTweaks[param_id][nucf]->GetBinContent(bin));
  }
}
