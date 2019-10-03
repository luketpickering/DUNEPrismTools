#include "dk2nu_TreeReader.hxx"

#include "GetUsage.hxx"
#include "HistogramBuilder.hxx"
#include "PhysicsUtility.hxx"
#include "ROOTUtility.hxx"

#include "TH2Jagged.h"

#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

inline double CosTheta(TVector3 const &v1, TVector3 const &v2) {
  return v1.Unit().Dot(v2.Unit());
}

std::tuple<double, double, double> GetNuWeight(DK2NuReader &dk2nuRdr,
                                               TVector3 const &DetPoint) {
  static const double detRadius = 100.0; // in cm
  static double const mumass = 0.105658389;

  double parent_mass = dk2nuRdr.GetParentMass();
  if (parent_mass == 0xdeadbeef) {
    return std::make_tuple(0xdeadbeef, 0xdeadbeef, 0xdeadbeef);
  }

  TLorentzVector parent_4mom;
  parent_4mom.SetXYZM(dk2nuRdr.decay_pdpx, dk2nuRdr.decay_pdpy,
                      dk2nuRdr.decay_pdpz, parent_mass);

  double enuzr = dk2nuRdr.decay_necm;

  TVector3 nuRay((DetPoint[0] - dk2nuRdr.decay_vx),
                 (DetPoint[1] - dk2nuRdr.decay_vy),
                 (DetPoint[2] - dk2nuRdr.decay_vz));

  if (!std::isnormal(nuRay.Mag())) {
    std::cout << "[ERROR]: Abnormal nuray, something is up." << std::endl;
    abort();
  }

  double parent_4mom_Gamma = parent_4mom.Gamma();
  double emrat =
      1.0 / (parent_4mom_Gamma *
             (1 - parent_4mom.Beta() * CosTheta(parent_4mom.Vect(), nuRay)));

  double nu_energy = emrat * enuzr;

  double sangDet = (1.0 - cos(atan(detRadius / nuRay.Mag()))) / 2.0;

  if (sangDet != sangDet) {
    std::cerr
        << "[ERROR]: Failed to calculate the solid angle element for nuRay: { "
        << nuRay[0] << ", " << nuRay[1] << ", " << nuRay[2] << " }."
        << std::endl;
    throw;
  }

  double nu_wght = sangDet * emrat * emrat;

  // done for all except polarized muon
  // in which case need to modify weight
  if (abs(dk2nuRdr.decay_ptype) == 13) {
    // boost new neutrino to mu decay cm
    TVector3 betaVect = parent_4mom.Vect();
    betaVect[0] /= parent_4mom.E();
    betaVect[1] /= parent_4mom.E();
    betaVect[2] /= parent_4mom.E();

    TVector3 nu_3mom = nuRay.Unit() * nu_energy;

    double partial = parent_4mom_Gamma * betaVect.Dot(nu_3mom);
    partial = nu_energy - partial / (parent_4mom_Gamma + 1.);

    TLorentzVector p_dcm_nu;
    for (int i = 0; i < 3; i++) {
      p_dcm_nu[i] = nu_3mom[i] - betaVect[i] * parent_4mom_Gamma * partial;
    }
    p_dcm_nu[3] = 0;
    for (int i = 0; i < 3; i++) {
      p_dcm_nu[3] += p_dcm_nu[i] * p_dcm_nu[i];
    }
    p_dcm_nu[3] = sqrt(p_dcm_nu[3]);

    // boost parent of mu to mu production cm
    parent_4mom_Gamma = dk2nuRdr.decay_ppenergy / parent_mass;
    TVector3 parentParentBetaVect(dk2nuRdr.decay_ppdxdz * dk2nuRdr.decay_pppz,
                                  dk2nuRdr.decay_ppdydz * dk2nuRdr.decay_pppz,
                                  dk2nuRdr.decay_pppz);
    parentParentBetaVect[0] /= dk2nuRdr.decay_ppenergy;
    parentParentBetaVect[1] /= dk2nuRdr.decay_ppenergy;
    parentParentBetaVect[2] /= dk2nuRdr.decay_ppenergy;

    double partial_mp =
        parent_4mom_Gamma * (betaVect[0] * dk2nuRdr.decay_muparpx +
                             betaVect[1] * dk2nuRdr.decay_muparpy +
                             betaVect[2] * dk2nuRdr.decay_muparpz);
    partial_mp = dk2nuRdr.decay_mupare - partial_mp / (parent_4mom_Gamma + 1.);

    TLorentzVector p_pcm_mp;
    p_pcm_mp[0] =
        dk2nuRdr.decay_muparpx - betaVect[0] * parent_4mom_Gamma * partial_mp;
    p_pcm_mp[1] =
        dk2nuRdr.decay_muparpy - betaVect[1] * parent_4mom_Gamma * partial_mp;
    p_pcm_mp[2] =
        dk2nuRdr.decay_muparpz - betaVect[2] * parent_4mom_Gamma * partial_mp;
    p_pcm_mp[3] = 0;
    for (int i = 0; i < 3; i++) {
      p_pcm_mp[3] += p_pcm_mp[i] * p_pcm_mp[i];
    }
    p_pcm_mp[3] = sqrt(p_pcm_mp[3]);

    double wt_ratio = 1.;
    // have to check p_pcm_mp
    // it can be 0 if mupar.=0. (I guess muons created in target??)
    if (p_pcm_mp[3] != 0.) {
      // calc new decay angle w.r.t. (anti)spin direction
      double costh = CosTheta(p_dcm_nu.Vect(), p_pcm_mp.Vect());

      // calc relative weight due to angle difference
      if (abs(dk2nuRdr.decay_ntype) == 12) {
        wt_ratio = 1. - costh;
      } else if (abs(dk2nuRdr.decay_ntype) == 14) {
        double xnu = 2. * enuzr / mumass;
        wt_ratio =
            ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
      } else {
        std::cout << "eventRates:: Bad neutrino type = " << dk2nuRdr.decay_ntype
                  << std::endl;
        throw;
      }
    }
    nu_wght *= wt_ratio;
  }

  if (!std::isnormal(nu_energy) || !std::isnormal(nu_wght)) {
    std::cout << "[ERROR]: Calculated bad nustats (E = " << nu_energy
              << ",  W = " << nu_wght << ")." << std::endl;
    abort();
  }
  return std::make_tuple(nu_energy, nuRay.Theta(), nu_wght);
}

enum OffAxisStepUnits { kPostion_m = 0, kmrad, kdegrees };
std::vector<int> NuPDGTargets = {-12, 12, -14, 14};
std::vector<std::string> NuPDGTargetNames = {"nuebar", "nue", "numubar",
                                             "numu"};

struct Config {

  struct Binning {
    std::vector<fhicl::ParameterSet> energy;
    std::vector<bool> is_off_axis;
    std::vector<bool> is_jagged;
    std::vector<std::vector<double>> off_axis;
    OffAxisStepUnits off_axis_type;
    std::string determine_binning_file;
    bool isdef_determine_binning_file_;
    size_t determine_binning_N;
    bool isdef_determine_binning_N_;
    double determine_binning_width;
    bool isdef_determine_binning_width_;
    double determine_binning_max_width;
    bool isdef_determine_binning_max_width_;
  };
  Binning binning;

  struct Flux_Window {
    double detector_half_height_cm;
    double z_from_target_cm;
    std::pair<double, double> z_limits_cm;
    double z_rel_min_cm, z_extent_cm;
  };
  Flux_Window flux_window;

  struct Input {
    size_t max_decay_parents;
    bool limit_decay_parent_use;
    bool use_dk2nu_lite;
    bool use_dk2nu_ppfx;
    bool isdef_dk2nu_ppfx_;
    bool use_dk2nu_ppfx_allweights;
    bool isdef_dk2nu_ppfx_allweights_;
    size_t number_ppfx_universes;
    bool isdef_number_ppfx_universes_;
    double ppfx_weightcap;

    std::string input_descriptor;
  };
  Input input;

  struct Output {
    int only_nu_species_pdg;
    bool isdef_only_nu_species_pdg_;

    bool separate_by_hadron_species;
    bool use_THF;
    bool use_reciprocal_energy;
    bool ignore_prediction_integral;
    std::string write_nu_ray_tree_to;

    std::string output_filename;
  };
  Output output;

  Config() {
    // Only need to set defaults for items that can be overriden by the CLI
    input.use_dk2nu_ppfx = false;
    input.isdef_dk2nu_ppfx_ = true;
    input.use_dk2nu_ppfx_allweights = false;
    input.isdef_dk2nu_ppfx_allweights_ = true;
    input.number_ppfx_universes = 0;
    input.isdef_number_ppfx_universes_ = true;

    output.only_nu_species_pdg = 0;
    output.isdef_only_nu_species_pdg_ = true;

    binning.determine_binning_file = "";
    binning.isdef_determine_binning_file_ = true;
    binning.determine_binning_N = 1000;
    binning.isdef_determine_binning_N_ = true;
    binning.determine_binning_width = 0.025;
    binning.isdef_determine_binning_width_ = true;
    binning.determine_binning_max_width = 2;
    binning.isdef_determine_binning_max_width_ = true;
  }

  void SetFromFHiCL(fhicl::ParameterSet const &ps) {

    fhicl::ParameterSet binning_ps =
        ps.get<fhicl::ParameterSet>("binning", fhicl::ParameterSet());

    if (!binning.isdef_determine_binning_file_) {
      binning.determine_binning_file =
          binning_ps.get<std::string>("optimized_binning_file", "");
    }
    if (!binning.isdef_determine_binning_N_) {
      binning.determine_binning_N =
          binning_ps.get<size_t>("optimized_binning_N", 1000);
    }
    if (!binning.isdef_determine_binning_width_) {
      binning.determine_binning_width =
          binning_ps.get<double>("optimized_binning_min_width", 0.25);
    }
    if (!binning.isdef_determine_binning_max_width_) {
      binning.determine_binning_max_width =
          binning_ps.get<double>("optimized_binning_max_width", 2);
    }

    binning.off_axis.resize(NuPDGTargetNames.size());
    std::fill_n(std::back_inserter(binning.is_off_axis), 4, false);
    std::fill_n(std::back_inserter(binning.is_jagged), 4, false);
    for (size_t nu_it = 0; nu_it < NuPDGTargetNames.size(); ++nu_it) {
      std::string const &nuname = NuPDGTargetNames[nu_it];

      fhicl::ParameterSet def;
      def.put<std::string>("off_axis", "0_0:1");
      def.put<std::string>("energy", "0_20:0.2");

      binning.energy.push_back(
          binning_ps.get<fhicl::ParameterSet>(nuname, def));

      bool has_off_axis = (binning.energy.back().get<std::string>(
                               "off_axis", "0_0:1") != "0_0:1");
      bool is_jagged = binning.energy.back().has_key("xattr");

      if (has_off_axis || is_jagged) {
        binning.is_off_axis[nu_it] = true;
        if (has_off_axis) {
          binning.off_axis[nu_it] =
              BuildBinEdges(binning.energy.back().get<std::string>("off_axis"));
        } else { // IsJagged
          binning.is_jagged[nu_it] = true;
          binning.off_axis[nu_it] =
              ParseJaggedUniformAxis(binning.energy.back());
        }
      } else {
        binning.off_axis[nu_it] = std::vector<double>{-2.5, 2.5};
      }
    }

    std::string off_axis_step =
        binning_ps.get<std::string>("type", "position_m");
    if (off_axis_step == "mrad") {
      binning.off_axis_type = kmrad;
    } else if (off_axis_step == "degrees") {
      binning.off_axis_type = kdegrees;
    } else if (off_axis_step == "position_m") {
      binning.off_axis_type = kPostion_m;
    } else {
      std::cout
          << "[ERROR]: Read \"binning.off_axis.type: " << off_axis_step
          << "\", but only understand one of: [mrad, degrees, <position_m>]."
          << std::endl;
      exit(1);
    }

    fhicl::ParameterSet flux_window_ps =
        ps.get<fhicl::ParameterSet>("flux_window", fhicl::ParameterSet());

    flux_window.detector_half_height_cm =
        flux_window_ps.get<double>("height_m") * 1.0E2 / 2.0;

    if (flux_window_ps.has_key("z_from_target_m")) {
      flux_window.z_from_target_cm =
          flux_window_ps.get<double>("z_from_target_m") * 1.0E2;
    } else if (flux_window_ps.has_key("z_from_target_km")) {
      flux_window.z_from_target_cm =
          flux_window_ps.get<double>("z_from_target_km") * 1.0E5;
    } else {
      std::cout
          << "[ERROR]: Expected to find key: \"flux_window.z_from_target_m\" "
             "or \"flux_window.z_from_target_km\"."
          << std::endl;
      throw;
    }

    if (flux_window_ps.has_key("width_m")) {

      double half_width_m = flux_window_ps.get<double>("width_m") / 2.0;
      for (size_t nu_it = 0; nu_it < NuPDGTargetNames.size(); ++nu_it) {
        if (binning.is_off_axis[nu_it]) {
          std::cout << "[WARN]: Found flux_window.width_m but "
                       "an off axis binning was provided for nu species: "
                    << NuPDGTargetNames[nu_it]
                    << ", the window width will be ignored." << std::endl;
        } else {
          binning.off_axis[nu_it] =
              std::vector<double>{-half_width_m, half_width_m};
        }
      }
    }
    for (size_t nu_it = 0; nu_it < NuPDGTargetNames.size(); ++nu_it) {
      if (!binning.off_axis[nu_it].size()) {
        std::cout << "[ERROR]: No off axis binning or window width found for "
                     "species: "
                  << NuPDGTargetNames[nu_it] << std::endl;
        throw;
      }
    }

    // Default to using a plane rather than a volume
    flux_window.z_limits_cm =
        flux_window_ps.get<std::pair<double, double>>("z_limits_cm", {0, 0});

    flux_window.z_rel_min_cm = flux_window.z_limits_cm.first;
    flux_window.z_extent_cm =
        flux_window.z_limits_cm.second - flux_window.z_limits_cm.first;

    fhicl::ParameterSet input_ps =
        ps.get<fhicl::ParameterSet>("input", fhicl::ParameterSet());

    input.max_decay_parents = input_ps.get<size_t>(
        "max_decay_parents", std::numeric_limits<size_t>::max());
    input.limit_decay_parent_use =
        input_ps.get<bool>("limit_decay_parent_use", false);

    input.use_dk2nu_lite = input_ps.get<bool>("use_dk2nu_lite", true);

    if (input.isdef_dk2nu_ppfx_) {
      input.use_dk2nu_ppfx = input_ps.get<bool>("use_dk2nu_ppfx", false);
    }
    if (input.isdef_dk2nu_ppfx_allweights_) {
      input.use_dk2nu_ppfx_allweights =
          input_ps.get<bool>("use_dk2nu_ppfx_allweights", false);
    }

    if (input.isdef_number_ppfx_universes_) {
      input.number_ppfx_universes =
          input_ps.get<size_t>("number_ppfx_universes", 100);
    }

    input.ppfx_weightcap = input_ps.get<double>("ppfx_weightcap", 10);

    fhicl::ParameterSet output_ps =
        ps.get<fhicl::ParameterSet>("output", fhicl::ParameterSet());

    if (output.isdef_only_nu_species_pdg_) {
      output.only_nu_species_pdg = output_ps.get<int>("only_nu_species_pdg", 0);
    }

    output.separate_by_hadron_species =
        output_ps.get<bool>("separate_by_hadron_species", false);
    output.use_THF = output_ps.get<bool>("use_THF", false);
    output.use_reciprocal_energy =
        output_ps.get<bool>("use_reciprocal_energy", false);
    output.ignore_prediction_integral =
        output_ps.get<bool>("ignore_prediction_integral", false);
    output.write_nu_ray_tree_to =
        output_ps.get<std::string>("write_nu_ray_tree_to", "");
  }
};

Config config;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "flux_tools") << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  bool got_fhicl_config;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      config.input.input_descriptor = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      config.output.output_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "--PPFX") {
      config.input.use_dk2nu_ppfx = true;
      config.input.isdef_dk2nu_ppfx_ = false;
    } else if (std::string(argv[opt]) == "--PPFX-Components") {
      config.input.use_dk2nu_ppfx_allweights = true;
      config.input.isdef_dk2nu_ppfx_allweights_ = false;
    } else if (std::string(argv[opt]) == "--NPPFX-Universes") {
      config.input.number_ppfx_universes = str2T<size_t>(argv[++opt]);
      config.input.isdef_number_ppfx_universes_ = false;
    } else if (std::string(argv[opt]) == "--only-pdg") {
      config.output.only_nu_species_pdg = str2T<int>(argv[++opt]);
      config.output.isdef_only_nu_species_pdg_ = false;
    } else if (std::string(argv[opt]) == "--optimized-binning-file") {
      config.binning.determine_binning_file = argv[++opt];
      config.binning.isdef_determine_binning_file_ = false;
    } else if (std::string(argv[opt]) == "--optimized-binning-N") {
      config.binning.determine_binning_N = str2T<size_t>(argv[++opt]);
      config.binning.isdef_determine_binning_N_ = false;
    } else if (std::string(argv[opt]) == "--optimized-binning-min-width") {
      config.binning.determine_binning_width = str2T<double>(argv[++opt]);
      config.binning.isdef_determine_binning_width_ = false;
    } else if (std::string(argv[opt]) == "--optimized-binning-max-width") {
      config.binning.determine_binning_max_width = str2T<double>(argv[++opt]);
      config.binning.isdef_determine_binning_max_width_ = false;
    } else if (std::string(argv[opt]) == "--fhicl") {
      config.SetFromFHiCL(fhicl::make_ParameterSet(argv[++opt]));
      got_fhicl_config = true;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
  if (!got_fhicl_config) {
    std::cout << "[ERROR]: Expected to find fhicl config file specifier."
              << std::endl;
    SayUsage(argv);
    exit(1);
  }
}

constexpr double rad2deg = 90.0 / asin(1);

std::pair<double, TVector3>
GetRandomFluxWindowPosition(TRandom3 &rnjesus, double win_min, double win_max) {
  double OffAxisCenter = (win_max + win_min) / 2.0;
  double OffAxisHalfRange = (win_max - win_min) / 2.0;

  double zpos_cm = config.flux_window.z_from_target_cm +
                   config.flux_window.z_rel_min_cm +
                   config.flux_window.z_extent_cm * rnjesus.Uniform();

  TVector3 rndDetPos(0,
                     (2.0 * rnjesus.Uniform() - 1.0) *
                         config.flux_window.detector_half_height_cm,
                     zpos_cm);

  double RandomOffAxisPos =
      OffAxisCenter + (2.0 * rnjesus.Uniform() - 1.0) * OffAxisHalfRange;
  switch (config.binning.off_axis_type) {
  case kPostion_m: {
    rndDetPos[0] = RandomOffAxisPos * 100.0; // to cm
    break;
  }
  case kmrad: {
    rndDetPos[0] =
        tan(RandomOffAxisPos * 1E-3) * config.flux_window.z_from_target_cm;
    break;
  }
  case kdegrees: {
    static const double deg2rad = asin(1) / 90.0;
    rndDetPos[0] =
        tan(RandomOffAxisPos * deg2rad) * config.flux_window.z_from_target_cm;
    break;
  }
  }

  if (!std::isnormal(rndDetPos.Mag())) {
    std::cout << "[ERROR]: Abnormal det pos, something is up." << std::endl;
    abort();
  }

  return std::make_pair(RandomOffAxisPos, rndDetPos);
}

void AllInOneGo(DK2NuReader &dk2nuRdr, double TotalPOT) {

  TRandom3 rnjesus;

  TFile *nuray_file = nullptr;
  std::vector<double> nuray_weights;
  std::vector<double> nuray_energy;
  std::vector<double> nuray_theta;
  TTree *nuray_tree = nullptr;
  if (config.output.write_nu_ray_tree_to.length()) {
    nuray_file = CheckOpenFile(config.output.write_nu_ray_tree_to, "RECREATE");
    nuray_tree = new TTree("nuray", "");
    nuray_tree->Branch("flux_ray_weight", &nuray_weights);
    nuray_tree->Branch("flux_ray_energy", &nuray_energy);
    nuray_tree->Branch("flux_ray_theta", &nuray_theta);
  }

  size_t NDecayParents =
      std::min(config.input.max_decay_parents, size_t(dk2nuRdr.GetEntries()));
  size_t NDecayParentsUsed = 0;

  double TotalPOTGuess =
      TotalPOT * (double(NDecayParents) / double(dk2nuRdr.GetEntries()));
  std::cout << "Only using the first " << NDecayParents << " events out of "
            << dk2nuRdr.GetEntries() << ", scaling POT to ~" << TotalPOTGuess
            << std::endl;
  std::cout << "Reading " << NDecayParents << " Dk2Nu entries." << std::endl;

  bool determine_binning = config.binning.determine_binning_file.size();

  TFile *outfile = nullptr;
  if (!determine_binning) {
    outfile = CheckOpenFile(config.output.output_filename, "RECREATE");
  }

  if (config.output.only_nu_species_pdg) {
    std::vector<int>::iterator pdg_it =
        std::find(NuPDGTargets.begin(), NuPDGTargets.end(),
                  config.output.only_nu_species_pdg);
    if (pdg_it == NuPDGTargets.end()) {
      std::cout << "[ERROR]: config.output.only_nu_species_pdg = "
                << config.output.only_nu_species_pdg << " is invalid."
                << std::endl;
      throw;
    }
    size_t idx = std::distance(NuPDGTargets.begin(), pdg_it);
    config.binning.energy =
        std::vector<fhicl::ParameterSet>{config.binning.energy[idx]};
    config.binning.is_off_axis =
        std::vector<bool>{config.binning.is_off_axis[idx]};
    config.binning.off_axis =
        std::vector<std::vector<double>>{config.binning.off_axis[idx]};

    NuPDGTargets = std::vector<int>{NuPDGTargets[idx]};
  } else if (determine_binning) {
    std::cout << "[ERROR]: When determining binning, expect only single "
                 "neutrino species, please specify output.only_nu_species_pdg "
                 "= <pdg> in configuration."
              << std::endl;
    abort();
  }

  std::vector<std::vector<std::vector<TH1 *>>> Hists;
  std::vector<std::vector<std::vector<TH2 *>>> Hists_2D;

  size_t NPPFXU = 1;
  if (config.input.use_dk2nu_ppfx) {
    if (config.input.use_dk2nu_ppfx_allweights) {
      NPPFXU =
          (config.input.number_ppfx_universes * DK2NuReader::kNPPFXAllWeights) +
          2;
    } else {
      NPPFXU = config.input.number_ppfx_universes + 2;
    }
  }
  bool use_PPFX = (NPPFXU > 1);
  if (determine_binning) {
    NPPFXU = 1;
  }

  if (!determine_binning) {
    Hists.resize(NPPFXU);
    Hists_2D.resize(NPPFXU);
  }

  for (size_t ppfx_univ_it = 0; (!determine_binning) && (ppfx_univ_it < NPPFXU);
       ++ppfx_univ_it) {

    Hists[ppfx_univ_it].resize(NuPDGTargets.size());
    Hists_2D[ppfx_univ_it].resize(NuPDGTargets.size());

    for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
      std::stringstream hist_name("");
      hist_name << "LBNF_" << GetSpeciesName(NuPDGTargets[nuPDG_it]) << "_flux";

      Int_t IsOffAxis = config.binning.is_off_axis[nuPDG_it];

      if (use_PPFX) {
        hist_name << GetPPFXHistName(ppfx_univ_it,
                                     config.input.number_ppfx_universes);
      }

      std::string div_hist_name = hist_name.str() + "_div";

      std::string hist_title =
          config.output.use_reciprocal_energy
              ? ";1/#it{E}_{#nu} (GeV^{-1});#Phi_{#nu} (cm^{-2} per POT per 1 "
                "GeV^{-1})"
              : ";#it{E}_{#nu} (GeV);#Phi_{#nu} (cm^{-2} per POT per 1 GeV)";

      std::string div_hist_title =
          config.output.use_reciprocal_energy
              ? ";1/#it{E}_{#nu} (GeV^{-1});Incoming Neutrino #theta "
                "(mrad);#Phi_{#nu} (cm^{-2} per POT per 1 GeV^{-1})"
              : ";#it{E}_{#nu} (GeV);Incoming Neutrino #theta "
                "(mrad);#Phi_{#nu} (cm^{-2} per POT per 1 GeV)";
      if (IsOffAxis) {
        switch (config.binning.off_axis_type) {
        case kPostion_m: {
          hist_title = config.output.use_reciprocal_energy
                           ? ";1/#it{E}_{#nu} (GeV^{-1});Off-axis postion "
                             "(m);#Phi_{#nu} (cm^{-2} per POT per 1 GeV^{-1})"
                           : ";#it{E}_{#nu} (GeV);Off-axis postion "
                             "(m);#Phi_{#nu} (cm^{-2} per POT per 1 GeV)";
          div_hist_title =
              config.output.use_reciprocal_energy
                  ? ";1/#it{E}_{#nu} (GeV^{-1});Incoming Neutrino #theta "
                    "(mrad);Off-axis postion (m);#Phi_{#nu} (cm^{-2} per "
                    "POT per 1 GeV^{-1})"
                  : ";#it{E}_{#nu} (GeV);Incoming Neutrino #theta "
                    "(mrad);Off-axis postion (m);#Phi_{#nu} (cm^{-2} per "
                    "POT per 1 GeV)";
          break;
        }
        case kmrad: {
          hist_title =
              config.output.use_reciprocal_energy
                  ? ";1/#it{E}_{#nu} (GeV^{-1});Off-axis angle "
                    "(mrad);#Phi_{#nu} (cm^{-2} per POT per 1 GeV^{-1})"
                  : ";#it{E}_{#nu} (GeV);Off-axis angle (mrad);#Phi_{#nu} "
                    "(cm^{-2} per POT per 1 GeV)";
          div_hist_title = ";#it{E}_{#nu} (GeV);Incoming Neutrino #theta "
                           "(mrad);Off-axis angle (mrad);#Phi_{#nu} "
                           "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        case kdegrees: {
          hist_title =
              config.output.use_reciprocal_energy
                  ? "1/;#it{E}_{#nu} (GeV^{-1});Off-axis angle "
                    "(degrees);#Phi_{#nu} (cm^{-2} per POT per 1 GeV^{-1})"
                  : ";#it{E}_{#nu} (GeV);Off-axis angle "
                    "(degrees);#Phi_{#nu} (cm^{-2} per POT per 1 GeV)";
          div_hist_title =
              config.output.use_reciprocal_energy
                  ? "1/;#it{E}_{#nu} (GeV^{-1});Incoming Neutrino #theta "
                    "(mrad);Off-axis angle (degrees);#Phi_{#nu} (cm^{-2} "
                    "per POT per 1 GeV^{-1})"
                  : ";#it{E}_{#nu} (GeV);Incoming Neutrino #theta "
                    "(mrad);Off-axis angle (degrees);#Phi_{#nu} (cm^{-2} "
                    "per POT per 1 GeV)";
          break;
        }
        }
      }

      std::vector<std::string> parent_suffixes = {"", "pi", "k", "k0", "mu"};
      for (size_t ps_it = 0; ps_it < parent_suffixes.size(); ++ps_it) {
        std::string parent_suffix = parent_suffixes[ps_it];
        if (parent_suffix.size()) {
          if (!config.output.separate_by_hadron_species) {
            continue;
          }
          parent_suffix = std::string("_") + parent_suffix;
        }

        if (!ppfx_univ_it) {
          if (config.output.use_THF) {
            Hists[ppfx_univ_it][nuPDG_it].push_back(
                BuildHistogram<float>(config.binning.energy[nuPDG_it],
                                      hist_name.str() + parent_suffix,
                                      hist_title, "energy", "off_axis"));
          } else {
            Hists[ppfx_univ_it][nuPDG_it].push_back(
                BuildHistogram<double>(config.binning.energy[nuPDG_it],
                                       hist_name.str() + parent_suffix,
                                       hist_title, "energy", "off_axis"));
          }
        } else { // We can just clone if we already have one for this species.
          Hists[ppfx_univ_it][nuPDG_it].push_back(
              static_cast<TH1 *>(Hists[0][nuPDG_it][ps_it]->Clone(
                  (hist_name.str() + parent_suffix).c_str())));
        }

        if (IsOffAxis) {
          Hists_2D[ppfx_univ_it][nuPDG_it].push_back(
              static_cast<TH2 *>(Hists[ppfx_univ_it][nuPDG_it].back()));
        }
      }
    }
  } // End loop over making many histos

  std::vector<size_t> NNuPDGTargets;
  for (auto i : NuPDGTargets) {
    (void)i;
    NNuPDGTargets.push_back(0);
  }

  std::vector<std::tuple<float, float, double>> binning_opt_points;

  size_t updateStep = (NDecayParents / 10) ? NDecayParents / 10 : 1;
  for (size_t nu_it = 0; nu_it < NDecayParents; ++nu_it) {
    if (!(nu_it % updateStep)) {
      std::cout << "--" << nu_it << "/" << NDecayParents << std::endl;
    }

    dk2nuRdr.GetEntry(nu_it);

    int nuPDG_it =
        std::distance(NuPDGTargets.begin(),
                      std::find(NuPDGTargets.begin(), NuPDGTargets.end(),
                                dk2nuRdr.decay_ntype));

    Int_t NOffAxisBins = Int_t(config.binning.off_axis[nuPDG_it].size() - 1);
    Int_t IsOffAxis = config.binning.is_off_axis[nuPDG_it];

    nuray_weights.clear();
    nuray_energy.clear();
    nuray_theta.clear();
    // Set up space for ray wegihts
    if (config.output.only_nu_species_pdg &&
        (config.output.only_nu_species_pdg != dk2nuRdr.decay_ntype)) {
      if (nuray_tree) {
        for (Int_t ang_it = 0; ang_it < NOffAxisBins; ++ang_it) {
          nuray_weights.push_back(0);
          nuray_energy.push_back(0);
          nuray_theta.push_back(0);
          // If we aren't re-using the parents then we have placed this neutrino
          // randomly in the 2D range and should now move to the next one.
          if (config.input.limit_decay_parent_use) {
            break;
          }
        }
        nuray_tree->Fill();
      }
      // Also have to increment here to get POT scaling correct!
      NDecayParentsUsed++;
      continue;
    }

    double wF = (dk2nuRdr.decay_nimpwt / TMath::Pi());

    for (Int_t ang_it = 0; ang_it < NOffAxisBins; ++ang_it) {

      // If we are not re-using the decay parents, then this is placed randomly
      // over the whole off-axis range.
      double win_min = config.input.limit_decay_parent_use
                           ? config.binning.off_axis[nuPDG_it].front()
                           : config.binning.off_axis[nuPDG_it][ang_it];
      double win_max = config.input.limit_decay_parent_use
                           ? config.binning.off_axis[nuPDG_it].back()
                           : config.binning.off_axis[nuPDG_it][ang_it + 1];
      std::pair<double, TVector3> det_point =
          GetRandomFluxWindowPosition(rnjesus, win_min, win_max);

      std::tuple<double, double, double> nuStats =
          GetNuWeight(dk2nuRdr, det_point.second);

      if (determine_binning) {
        binning_opt_points.emplace_back(float(det_point.first),
                                        float(std::get<0>(nuStats)),
                                        std::get<2>(nuStats) * wF);
        continue;
      }

      if (config.output.use_reciprocal_energy) {
        std::get<0>(nuStats) = 1.0 / std::get<0>(nuStats);
      }

      double w = std::get<2>(nuStats) * wF;

      if (nuPDG_it == 4) {
        std::cout << "Warning, couldn't find plot index for NuPDG: "
                  << dk2nuRdr.decay_ntype << std::endl;
        exit(1);
      }

      if ((!std::isnormal(std::get<0>(nuStats)) || (!std::isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << ")." << std::endl;
        throw;
      }

      size_t extra_hist_index = 0;
      if (config.output.separate_by_hadron_species) {
        switch (dk2nuRdr.decay_ptype) {
        case 211:
        case -211: {
          extra_hist_index = 1;
          break;
        }
        case 321:
        case -321: {
          extra_hist_index = 2;
          break;
        }
        case 311:
        case 310:
        case 130: {
          extra_hist_index = 3;
          break;
        }
        case 13:
        case -13: {
          extra_hist_index = 4;
          break;
        }
        default: {
          std::cout << "[ERROR]: Unexpected decay_ptype: "
                    << dk2nuRdr.decay_ptype << std::endl;
          exit(1);
        }
        }
      }

      if (ang_it == 0) {
        NNuPDGTargets[nuPDG_it]++;
      }

      Int_t JagBin = config.binning.is_jagged[nuPDG_it]
                         ? Hists_2D[0][nuPDG_it][0]->FindFixBin(
                               std::get<0>(nuStats), det_point.first)
                         : 0; // Haven't checked
      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {

        double ppfx_w = 1;
        if (use_PPFX) {
          ppfx_w = std::min(config.input.ppfx_weightcap,
                            GetPPFXWeight(ppfx_univ_it,
                                          config.input.number_ppfx_universes,
                                          dk2nuRdr));
        }

        double tw = w * ppfx_w;

        // AddBinContent with known bin fails to track errors, must fill each
        // time.
        if (IsOffAxis) {
          if (config.binning.is_jagged[nuPDG_it]) { // but with TH2Jagged we can
                                                    // fill a known bin and keep
                                                    // the stats correct
            if (config.output.use_THF) {
              static_cast<TH2JaggedF *>(Hists_2D[ppfx_univ_it][nuPDG_it][0])
                  ->FillKnownBin(JagBin, tw);
            } else {
              static_cast<TH2JaggedD *>(Hists_2D[ppfx_univ_it][nuPDG_it][0])
                  ->FillKnownBin(JagBin, tw);
            }
          } else {
            Hists_2D[ppfx_univ_it][nuPDG_it][0]->Fill(std::get<0>(nuStats),
                                                      det_point.first, tw);
          }
        } else {
          Hists[ppfx_univ_it][nuPDG_it][0]->Fill(std::get<0>(nuStats), tw);
        }
        if (config.output.separate_by_hadron_species) {
          if (IsOffAxis) {
            if (config.binning
                    .is_jagged[nuPDG_it]) { // but with TH2Jagged we
                                            // can fill a known bin and
                                            // keep the stats correct
              if (config.output.use_THF) {
                static_cast<TH2JaggedF *>(
                    Hists_2D[ppfx_univ_it][nuPDG_it][extra_hist_index])
                    ->FillKnownBin(JagBin, tw);
              } else {
                static_cast<TH2JaggedD *>(
                    Hists_2D[ppfx_univ_it][nuPDG_it][extra_hist_index])
                    ->FillKnownBin(JagBin, tw);
              }
            } else {
              Hists_2D[ppfx_univ_it][nuPDG_it][extra_hist_index]->Fill(
                  std::get<0>(nuStats), det_point.first, tw);
            }
          } else {
            Hists[ppfx_univ_it][nuPDG_it][extra_hist_index]->Fill(
                std::get<0>(nuStats), tw);
          }
        }

      } // End loop over PPFX universes

      // If we aren't re-using the parents then we have placed this neutrino
      // randomly in the 2D range and should now move to the next one.
      if (config.input.limit_decay_parent_use) {
        break;
      }

      nuray_weights.push_back(w);
      nuray_energy.push_back(std::get<0>(nuStats));
      nuray_theta.push_back(std::get<1>(nuStats) * 1E3);

    } // End loop over angle bins
    if (nuray_tree) {
      nuray_tree->Fill();
    }
    NDecayParentsUsed++;
  } // End loop over decay parents

  if (determine_binning) {
    std::vector<TAxis> axes = DetermineOptimalNonUniformBinning(
        TAxis(config.binning.off_axis[0].size() - 1,
              config.binning.off_axis[0].data()),
        binning_opt_points, config.binning.determine_binning_N, 0, 20,
        config.binning.determine_binning_width,
        config.binning.determine_binning_max_width);

    std::ofstream fo(config.binning.determine_binning_file.c_str());

    fo << "binning: [ " << std::endl;
    for (size_t ax_it = 0; ax_it < axes.size(); ++ax_it) {
      fo << "  {\n  OA_m: \"" << config.binning.off_axis[0][ax_it] << ", "
         << config.binning.off_axis[0][ax_it + 1] << "\"" << std::endl;
      fo << "  E_GeV: \"";
      for (Int_t b_it = 0; b_it < axes[ax_it].GetNbins(); ++b_it) {
        fo << axes[ax_it].GetBinLowEdge(b_it + 1) << ",";
      }
      fo << axes[ax_it].GetBinUpEdge(axes[ax_it].GetNbins()) << "\" },"
         << std::endl;
    }
    fo << "  ]" << std::endl;

    return;
  }

  double TotalPOTUsed =
      TotalPOT * (double(NDecayParentsUsed) / double(dk2nuRdr.GetEntries()));
  std::cout << "Only used the first " << NDecayParentsUsed << " events out of "
            << dk2nuRdr.GetEntries() << ", scaling POT to " << TotalPOTUsed
            << std::endl;

  double MToCMAndPOTUsedScale = 1.0E-4 / TotalPOTUsed;

  for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
    std::cout << "[INFO]: Found " << NNuPDGTargets[nuPDG_it]
              << " Neutrinos with PDG: " << NuPDGTargets[nuPDG_it] << std::endl;

    if (!config.output.ignore_prediction_integral) {
      double integ = Hists[0][nuPDG_it][0]->Integral();
      if (!std::isnormal(integ)) {
        std::cerr << "[ERROR]: Flux for PDG: " << NuPDGTargets[nuPDG_it]
                  << " has bad integral (" << integ << ")" << std::endl;
        throw;
      }
    }

    Int_t IsOffAxis = config.binning.is_off_axis[nuPDG_it];
    // Scale to /cm^2 and per GeV (xbin width)
    if (IsOffAxis) {
      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {
        std::vector<TH2 *> &univ_hists = Hists_2D[ppfx_univ_it][nuPDG_it];
        for (size_t hist_index = 0;
             hist_index < (config.output.separate_by_hadron_species ? 5 : 1);
             ++hist_index) {
          TH2 *hist = univ_hists[hist_index];
          hist->Scale(MToCMAndPOTUsedScale, "width");
          // To accomodate TH2Jagged it is easier to scale by cell area and then
          // loop through Y and scale up by y;
          for (Int_t ybin_it = 0; ybin_it < hist->GetYaxis()->GetNbins();
               ++ybin_it) {
            double ybw = hist->GetYaxis()->GetBinWidth(ybin_it + 1);

            TAxis const *xax;
            if (config.output.use_THF) {
              xax = config.binning.is_jagged[nuPDG_it]
                        ? dynamic_cast<TH2JaggedF *>(hist)->GetNonUniformAxis(
                              ybin_it + 1)
                        : hist->GetXaxis();
            } else {
              xax = config.binning.is_jagged[nuPDG_it]
                        ? dynamic_cast<TH2JaggedD *>(hist)->GetNonUniformAxis(
                              ybin_it + 1)
                        : hist->GetXaxis();
            }

            for (Int_t xbin_it = 0; xbin_it < xax->GetNbins(); ++xbin_it) {

              Int_t gbin = hist->GetBin(xbin_it + 1, ybin_it + 1);

              double bc = hist->GetBinContent(gbin);
              double be = hist->GetBinError(gbin);
              hist->SetBinContent(gbin, bc * ybw);
              hist->SetBinError(gbin, be * ybw);
            } // Loop over x bins
          }   // Loop over y bins
          outfile->WriteTObject(
              Hists[ppfx_univ_it][nuPDG_it][hist_index],
              Hists[ppfx_univ_it][nuPDG_it][hist_index]->GetName());
          Hists[ppfx_univ_it][nuPDG_it][hist_index]->SetDirectory(nullptr);
          delete Hists[ppfx_univ_it][nuPDG_it][hist_index];
        } // Loop over decay parents
      }   // Loop over universes
    } else {
      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {
        for (size_t hist_index = 0;
             hist_index < (config.output.separate_by_hadron_species ? 5 : 1);
             ++hist_index) {
          Hists[ppfx_univ_it][nuPDG_it][hist_index]->Scale(MToCMAndPOTUsedScale,
                                                           "width");
          outfile->WriteTObject(
              Hists[ppfx_univ_it][nuPDG_it][hist_index],
              Hists[ppfx_univ_it][nuPDG_it][hist_index]->GetName());
          Hists[ppfx_univ_it][nuPDG_it][hist_index]->SetDirectory(nullptr);
          delete Hists[ppfx_univ_it][nuPDG_it][hist_index];
        }
      }
    }
  } // Loop over neutrino species

  if (nuray_file) {
    nuray_file->Write();
    nuray_file->Close();
  }

  outfile->Write();
  outfile->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  DK2NuReader *dk2nuRdr = new DK2NuReader(
      config.input.use_dk2nu_lite ? "dk2nuTree_lite" : "dk2nuTree",
      config.input.input_descriptor, config.input.use_dk2nu_lite);

  if (config.input.use_dk2nu_ppfx) {
    dk2nuRdr->SetPPFXBranchAddresses(config.input.number_ppfx_universes);
  }

  if (!dk2nuRdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }

  DKMetaReader *dkmRdr = new DKMetaReader(
      config.input.use_dk2nu_lite ? "dkmetaTree_lite" : "dkmetaTree",
      config.input.input_descriptor, config.input.use_dk2nu_lite);

  int metaNEntries = dkmRdr->GetEntries();

  double TotalPOT = 0;
  for (int i = 0; i < metaNEntries; ++i) {
    dkmRdr->GetEntry(i);
    TotalPOT += dkmRdr->pots;
  }

  std::cout << "Total POT: " << TotalPOT << std::endl;
  AllInOneGo(*dk2nuRdr, TotalPOT);
  dk2nuRdr->DumpMaxWeights();
}
