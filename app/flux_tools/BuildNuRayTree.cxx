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

#ifdef USE_DK2NU
#include "dk2nu/tree/calcLocationWeights.h"
#endif

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

// If we have dk2nu support, use that to calculate the decay weight
#ifdef USE_DK2NU
  bsim::Dk2Nu const &dk2nu = dk2nuRdr.GetDk2Nu();

  double nu_energy, nu_wght;
  if (bsim::calcEnuWgt(dk2nu.decay, DetPoint, nu_energy, nu_wght)) {
    std::cout << "[ERROR]: bsim::calcLocationWeight failed." << std::endl;
    throw;
  }

  TVector3 nuRay((DetPoint[0] - dk2nu.decay.vx), (DetPoint[1] - dk2nu.decay.vy),
                 (DetPoint[2] - dk2nu.decay.vz));

  return std::make_tuple(nu_energy, nuRay.Theta(), nu_wght);

#else

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
#endif
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

    bool use_dk2nu_lite;
    bool limit_decay_parent_use;
    bool isdef_limit_decay_parent_use_;

    std::string input_descriptor;
  };
  Input input;

  struct Output {
    int only_nu_species_pdg;
    bool isdef_only_nu_species_pdg_;

    bool separate_by_hadron_species;
    bool use_THF;
    bool ignore_prediction_integral;

    std::string output_filename;
  };
  Output output;

  Config() {
    // Only need to set defaults for items that can be overriden by the CLI

    output.only_nu_species_pdg = 0;
    output.isdef_only_nu_species_pdg_ = true;
    input.isdef_limit_decay_parent_use_ = true;
  }

  void SetFromFHiCL(fhicl::ParameterSet const &ps) {

    fhicl::ParameterSet binning_ps =
        ps.get<fhicl::ParameterSet>("binning", fhicl::ParameterSet());

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
    input.use_dk2nu_lite = input_ps.get<bool>("use_dk2nu_lite", true);
    if (input.isdef_limit_decay_parent_use_) {
      input.limit_decay_parent_use =
          input_ps.get<bool>("limit_decay_parent_use", false);
    }

    fhicl::ParameterSet output_ps =
        ps.get<fhicl::ParameterSet>("output", fhicl::ParameterSet());

    if (output.isdef_only_nu_species_pdg_) {
      output.only_nu_species_pdg = output_ps.get<int>("only_nu_species_pdg", 0);
    }
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
    } else if (std::string(argv[opt]) == "--only-pdg") {
      config.output.only_nu_species_pdg = str2T<int>(argv[++opt]);
      config.output.isdef_only_nu_species_pdg_ = false;
    } else if (std::string(argv[opt]) == "-N") {
      config.input.max_decay_parents = str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "--limit-decay-parent") {
      config.input.limit_decay_parent_use = true;
      config.input.isdef_limit_decay_parent_use_ = false;

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

struct NuRay {
  float phase_space_weight;
  float imp_weight;
  short pdg;
  float four_mom[4];
  float ray_end_point[3];

  short parent_pdg;
  float parent_position[3];
};

void AllInOneGo(DK2NuReader &dk2nuRdr, double TotalPOT) {

  TRandom3 rnjesus(0);

  TFile *nuray_file = CheckOpenFile(config.output.output_filename, "RECREATE");
  TTree *nuray_tree = new TTree("nurays", "");

  NuRay nr;

  nuray_tree->Branch("phase_space_weight", &nr.phase_space_weight,
                     "phase_space_weight/F");
  nuray_tree->Branch("imp_weight", &nr.imp_weight, "imp_weight/F");
  nuray_tree->Branch("pdg", &nr.pdg, "pdg/S");
  nuray_tree->Branch("four_mom", &nr.four_mom, "four_mom[4]/F");
  nuray_tree->Branch("ray_end_point", &nr.ray_end_point, "ray_end_point[3]/F");
  nuray_tree->Branch("parent_pdg", &nr.parent_pdg, "parent_pdg/S");
  nuray_tree->Branch("parent_position", &nr.parent_position,
                     "parent_position[3]/F");

  size_t NDecayParents =
      std::min(config.input.max_decay_parents, size_t(dk2nuRdr.GetEntries()));

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
    NuPDGTargets = std::vector<int>{
        NuPDGTargets[*pdg_it],
    };
  }

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

    if (config.output.only_nu_species_pdg &&
        (config.output.only_nu_species_pdg != dk2nuRdr.decay_ntype)) {
      continue;
    }

    double wF = dk2nuRdr.decay_nimpwt;

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

      double w = std::get<2>(nuStats) * (wF / TMath::Pi());

      if ((!std::isnormal(std::get<0>(nuStats)) || (!std::isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << ")." << std::endl;
        abort();
      }

      TVector3 nuRay((det_point.second[0] - dk2nuRdr.decay_vx),
                     (det_point.second[1] - dk2nuRdr.decay_vy),
                     (det_point.second[2] - dk2nuRdr.decay_vz));

      nr.phase_space_weight = w;
      nr.imp_weight = wF;
      nr.pdg = dk2nuRdr.decay_ntype;
      nr.ray_end_point[0] = det_point.second.X();
      nr.ray_end_point[1] = det_point.second.Y();
      nr.ray_end_point[2] = det_point.second.Z();
      TVector3 nuMom = nuRay.Unit() * std::get<0>(nuStats);
      nr.four_mom[0] = nuMom.X();
      nr.four_mom[1] = nuMom.Y();
      nr.four_mom[2] = nuMom.Z();
      nr.four_mom[3] = std::get<0>(nuStats);
      nr.parent_pdg = dk2nuRdr.decay_ptype;
      nr.parent_position[0] = dk2nuRdr.decay_vx;
      nr.parent_position[1] = dk2nuRdr.decay_vy;
      nr.parent_position[2] = dk2nuRdr.decay_vz;

      nuray_tree->Fill();
      // If we aren't re-using the parents then we have placed this neutrino
      // randomly in the 2D range and should now move to the next one.
      if (config.input.limit_decay_parent_use) {
        break;
      }
    } // End loop over angle bins
  }   // End loop over decay parents

  nuray_file->Write();
  nuray_file->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  DK2NuReader *dk2nuRdr = new DK2NuReader(
      config.input.use_dk2nu_lite ? "dk2nuTree_lite" : "dk2nuTree",
      config.input.input_descriptor, config.input.use_dk2nu_lite);

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
