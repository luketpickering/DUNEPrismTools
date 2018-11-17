#include "dk2nu_TreeReader.hxx"

#include "GetUsage.hxx"
#include "PhysicsUtility.hxx"
#include "ROOTUtility.hxx"

#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
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
    std::cout << "[ERROR]: Calculated bad nustats." << std::endl;
    exit(1);
  }
  return std::make_tuple(nu_energy, nuRay.Theta(), nu_wght);
}

std::vector<std::vector<double>> OffAxisSteps;
enum OffAxisStepUnits { kPostion_m = 0, kmrad, kdegrees };
OffAxisStepUnits stepType;

std::string inpDir = ".";
bool DoExtra = false;

double detector_half_height = 0;
bool ReUseParents = true;
bool DK2NULite = false;
bool ReadDK2NULitePPFX = false;
UInt_t NPPFX_Universes = 100;
bool UseTHF = false;

double ZDist = 57400;
int NMaxNeutrinos = -1;

std::vector<std::vector<double>> EnergyBinning;

std::string outputFile;

std::string runPlanCfg, runPlanName = "";

int OnlySpecies = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "flux_tools") << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?" || std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-m") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
                << " " << argv[opt + 1] << std::endl;
      OffAxisSteps.push_back(BuildBinEdges(argv[++opt]));
      for (size_t i = 0; i < 3; ++i) {
        OffAxisSteps.push_back(OffAxisSteps.back());
      }
      stepType = kmrad;
    } else if (std::string(argv[opt]) == "-d") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
                << " " << argv[opt + 1] << std::endl;
      OffAxisSteps.push_back(BuildBinEdges(argv[++opt]));
      for (size_t i = 0; i < 3; ++i) {
        OffAxisSteps.push_back(OffAxisSteps.back());
      }
      stepType = kdegrees;
    } else if (std::string(argv[opt]) == "-x") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
                << " " << argv[opt + 1] << std::endl;
      OffAxisSteps.push_back(BuildBinEdges(argv[++opt]));
      for (size_t i = 0; i < 3; ++i) {
        OffAxisSteps.push_back(OffAxisSteps.back());
      }
      stepType = kPostion_m;
    } else if (std::string(argv[opt]) == "-e") {
      DoExtra = true;
    } else if (std::string(argv[opt]) == "-b") {
      EnergyBinning.clear();
      std::vector<double> binningDescriptor =
          ParseToVect<double>(argv[++opt], ",");
      if (binningDescriptor.size() != 3) {
        std::cout << "[ERROR]: Recieved " << binningDescriptor.size()
                  << " entrys for -b, expected 3." << std::endl;
        exit(1);
      }
      int NBins = int(binningDescriptor[0]);
      double BLow = binningDescriptor[1];
      double BUp = binningDescriptor[2];
      double bwidth = (BUp - BLow) / double(NBins);
      std::vector<double> ebins;
      ebins.push_back(BLow);
      for (Int_t i = 0; i < NBins; ++i) {
        ebins.push_back(ebins.back() + bwidth);
      }
      for (size_t i = 0; i < 4; ++i) {
        EnergyBinning.push_back(ebins);
      }
    } else if (std::string(argv[opt]) == "-vb") {
      EnergyBinning.clear();
      EnergyBinning.push_back(BuildBinEdges(argv[++opt]));
      for (size_t i = 0; i < 3; ++i) {
        EnergyBinning.push_back(EnergyBinning.back());
      }
    } else if (std::string(argv[opt]) == "-h") {
      detector_half_height = str2T<double>(argv[++opt]) / 2.0;
    } else if (std::string(argv[opt]) == "-n") {
      NMaxNeutrinos = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-z") {
      ZDist = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-P") {
      ReUseParents = false;
    } else if (std::string(argv[opt]) == "-S") {
      OnlySpecies = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-L") {
      DK2NULite = true;
    } else if (std::string(argv[opt]) == "--PPFX") {
      ReadDK2NULitePPFX = true;
    } else if (std::string(argv[opt]) == "--NPPFXU") {
      NPPFX_Universes = str2T<UInt_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "--fhicl") {

      fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[++opt]);

      //  std::vector<int> NuPDGTargets = {-12, 12, -14, 14};
      for (std::string const &nuname : {"nuebar", "nue", "numubar", "numu"}) {
        OffAxisSteps.push_back(BuildBinEdges(
            ps.get<std::string>(nuname + "_off_axis_binning", "0_0:1")));
        EnergyBinning.push_back(BuildBinEdges(
            ps.get<std::string>(nuname + "_energy_binning", "0_20:0.2")));
      }

      std::string off_axis_step =
          ps.get<std::string>("off_axis_step_type", "position_m");
      if (off_axis_step == "mrad") {
        stepType = kmrad;
      } else if (off_axis_step == "degrees") {
        stepType = kdegrees;
      } else if (off_axis_step == "position_m") {
        stepType = kPostion_m;
      } else {
        std::cout
            << "[ERROR]: Read \"off_axis_step: " << off_axis_step
            << "\", but only understand one of: [mrad, degrees, <position_m>]."
            << std::endl;
        exit(1);
      }

      DoExtra = ps.get<bool>("make_extra_plots", false);
      detector_half_height =
          ps.get<double>("flux_window_height_m") * 1.0E2 / 2.0;
      NMaxNeutrinos = ps.get<int>("max_decay_parents", -1);
      ZDist = ps.get<double>("flux_window_z_from_target_m") * 1.0E2;
      ReUseParents = !ps.get<bool>("limit_decay_parent_use", false);
      OnlySpecies = ps.get<int>("only_nu_species_pdg", 0);
      DK2NULite = ps.get<bool>("use_dk2nu_lite", true);
      ReadDK2NULitePPFX = ps.get<bool>("use_dk2nu_ppfx", true);
      NPPFX_Universes = ps.get<int>("number_ppfx_universes", 100);
      UseTHF = ps.get<bool>("use_THF", false);

    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

constexpr double rad2deg = 90.0 / asin(1);

std::pair<double, TVector3>
GetRandomFluxWindowPosition(TRandom3 &rnjesus,
                            std::vector<double> const &oasteps,
                            Int_t FluxWindow = -1) {
  double OffAxisCenter =
      ((FluxWindow == -1) ? (oasteps.back() + oasteps.front())
                          : (oasteps[FluxWindow + 1] + oasteps[FluxWindow])) /
      2.0;
  double OffAxisHalfRange =
      ((FluxWindow == -1) ? (oasteps.back() - oasteps.front())
                          : (oasteps[FluxWindow + 1] - oasteps[FluxWindow])) /
      2.0;

  TVector3 rndDetPos(0, (2.0 * rnjesus.Uniform() - 1.0) * detector_half_height,
                     ZDist);

  double RandomOffAxisPos =
      OffAxisCenter + (2.0 * rnjesus.Uniform() - 1.0) * OffAxisHalfRange;
  switch (stepType) {
  case kPostion_m: {
    rndDetPos[0] = RandomOffAxisPos * 100.0; // to cm
    break;
  }
  case kmrad: {
    rndDetPos[0] = tan(RandomOffAxisPos * 1E-3) * ZDist;
    break;
  }
  case kdegrees: {
    static const double deg2rad = asin(1) / 90.0;
    rndDetPos[0] = tan(RandomOffAxisPos * deg2rad) * ZDist;
    break;
  }
  }

  return std::make_pair(RandomOffAxisPos, rndDetPos);
}

void AllInOneGo(DK2NuReader &dk2nuRdr, double TotalPOT) {
  TRandom3 rnjesus;

  if (!OffAxisSteps.size()) {
    std::cout << "[ERROR]: No off-axis positions specified (Try `-x 0`)."
              << std::endl;
    throw;
  }
  size_t NNeutrinos = (NMaxNeutrinos == -1)
                          ? dk2nuRdr.GetEntries()
                          : std::min(NMaxNeutrinos, int(dk2nuRdr.GetEntries()));

  TotalPOT = TotalPOT * (double(NNeutrinos) / double(dk2nuRdr.GetEntries()));
  std::cout << "Only using the first " << NNeutrinos << " events out of "
            << dk2nuRdr.GetEntries() << ", scaling POT to " << TotalPOT
            << std::endl;
  std::cout << "Reding " << NNeutrinos << " Dk2Nu entries." << std::endl;

  TFile *outfile = CheckOpenFile(outputFile.c_str(), "RECREATE");

  std::vector<int> NuPDGTargets = {-12, 12, -14, 14};
  if (OnlySpecies) {
    std::vector<int>::iterator pdg_it =
        std::find(NuPDGTargets.begin(), NuPDGTargets.end(), OnlySpecies);
    if (pdg_it == NuPDGTargets.end()) {
      std::cout << "[ERROR]: OnlySpecies = " << OnlySpecies << " is invalid."
                << std::endl;
      throw;
    }
    size_t idx = std::distance(NuPDGTargets.begin(), pdg_it);
    EnergyBinning = std::vector<std::vector<double>>{EnergyBinning[idx]};
    OffAxisSteps = std::vector<std::vector<double>>{OffAxisSteps[idx]};
    NuPDGTargets = std::vector<int>{NuPDGTargets[idx]};
  }

  std::vector<std::vector<std::vector<TH1 *>>> Hists;
  std::vector<std::vector<std::vector<TH2 *>>> Hists_2D;
  size_t NPPFXU = (ReadDK2NULitePPFX ? NPPFX_Universes + 2 : 1);
  bool use_PPFX = (NPPFXU > 1);
  Hists.resize(NPPFXU);
  Hists_2D.resize(NPPFXU);

  for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {

    Hists[ppfx_univ_it].resize(NuPDGTargets.size());
    Hists_2D[ppfx_univ_it].resize(NuPDGTargets.size());

    for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
      std::stringstream hist_name("");
      hist_name << "LBNF_" << GetSpeciesName(NuPDGTargets[nuPDG_it]) << "_flux";

      Int_t NOffAxisBins = Int_t(OffAxisSteps[nuPDG_it].size()) - 1;

      if (use_PPFX) {
        if (ppfx_univ_it == 0) {
          hist_name << "_Nom";

        } else if (ppfx_univ_it == 1) {
          hist_name << "_CV";

        } else {
          hist_name << "_univ_" << (ppfx_univ_it - 2);
        }
      }

      std::string hist_title = ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                               "(cm^{-2} per POT per 1 GeV)";
      if (NOffAxisBins > 1) {
        switch (stepType) {
        case kPostion_m: {
          hist_title = ";#it{E}_{#nu} (GeV);Off-axis postion (m);#Phi_{#nu} "
                       "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        case kmrad: {
          hist_title = ";#it{E}_{#nu} (GeV);Off-axis angle (mrad);#Phi_{#nu} "
                       "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        case kdegrees: {
          hist_title =
              ";#it{E}_{#nu} (GeV);Off-axis angle (degrees);#Phi_{#nu} "
              "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        }
      }

      for (std::string parent_suffix : {"", "pi", "k", "k0", "mu"}) {
        if (parent_suffix.size()) {
          if (!DoExtra) {
            continue;
          }
          parent_suffix = std::string("_") + parent_suffix;
        }

        if (UseTHF) {
          Hists[ppfx_univ_it][nuPDG_it].push_back(
              (NOffAxisBins < 2)
                  ? new TH1F((hist_name.str() + parent_suffix).c_str(),
                             hist_title.c_str(),
                             (EnergyBinning[nuPDG_it].size() - 1),
                             EnergyBinning[nuPDG_it].data())
                  : static_cast<TH1 *>(
                        new TH2F((hist_name.str() + parent_suffix).c_str(),
                                 hist_title.c_str(),
                                 (EnergyBinning[nuPDG_it].size() - 1),
                                 EnergyBinning[nuPDG_it].data(),
                                 (OffAxisSteps[nuPDG_it].size() - 1),
                                 OffAxisSteps[nuPDG_it].data())));
        } else {
          Hists[ppfx_univ_it][nuPDG_it].push_back(
              (NOffAxisBins < 2)
                  ? new TH1D((hist_name.str() + parent_suffix).c_str(),
                             hist_title.c_str(),
                             (EnergyBinning[nuPDG_it].size() - 1),
                             EnergyBinning[nuPDG_it].data())
                  : static_cast<TH1 *>(
                        new TH2D((hist_name.str() + parent_suffix).c_str(),
                                 hist_title.c_str(),
                                 (EnergyBinning[nuPDG_it].size() - 1),
                                 EnergyBinning[nuPDG_it].data(),
                                 (OffAxisSteps[nuPDG_it].size() - 1),
                                 OffAxisSteps[nuPDG_it].data())));
        }

        if (NOffAxisBins > 1) {
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

  size_t updateStep = (NNeutrinos / 10) ? NNeutrinos / 10 : 1;
  for (size_t nu_it = 0; nu_it < NNeutrinos; ++nu_it) {
    if (!(nu_it % updateStep)) {
      std::cout << "--" << nu_it << "/" << NNeutrinos << std::endl;
    }

    dk2nuRdr.GetEntry(nu_it);

    if (OnlySpecies && (OnlySpecies != dk2nuRdr.decay_ntype)) {
      continue;
    }

    int nuPDG_it =
        std::distance(NuPDGTargets.begin(),
                      std::find(NuPDGTargets.begin(), NuPDGTargets.end(),
                                dk2nuRdr.decay_ntype));
    Int_t NOffAxisBins = Int_t(OffAxisSteps[nuPDG_it].size()) - 1;

    double wF = (dk2nuRdr.decay_nimpwt / TMath::Pi()) * (1.0 / TotalPOT);

    // If there are no off axis bins, then we want to just build flux on-axis,
    // this loop needs to happen once.
    for (Int_t ang_it = 0; (!NOffAxisBins) || (ang_it < NOffAxisBins);
         ++ang_it) {

      // If we are not re-using the decay parents, then this is placed randomly
      // over the whole off-axis range.
      std::pair<double, TVector3> det_point = GetRandomFluxWindowPosition(
          rnjesus, OffAxisSteps[nuPDG_it], ReUseParents ? ang_it : -1);

      std::tuple<double, double, double> nuStats =
          GetNuWeight(dk2nuRdr, det_point.second);

      double w = std::get<2>(nuStats) * wF;

      if (nuPDG_it == 4) {
        std::cout << "Warning, couldn't find plot index for NuPDG: "
                  << dk2nuRdr.decay_ntype << std::endl;
        exit(1);
      }

      if ((!std::isnormal(std::get<0>(nuStats)) || (!std::isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << "*" << TotalPOT << ")."
                  << std::endl;
        throw;
      }

      size_t extra_hist_index = 0;
      if (DoExtra) {
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

      Int_t bin_it;
      if (NOffAxisBins > 1) {
        Int_t bin_it_x =
            Hists_2D.front()[nuPDG_it].front()->GetXaxis()->FindFixBin(
                std::get<0>(nuStats));
        Int_t bin_it_y =
            Hists_2D.front()[nuPDG_it].front()->GetYaxis()->FindFixBin(
                det_point.first);
        bin_it = Hists_2D.front()[nuPDG_it].front()->GetBin(bin_it_x, bin_it_y);
      } else {
        bin_it = Hists.front()[nuPDG_it].front()->GetXaxis()->FindFixBin(
            std::get<0>(nuStats));
      }

      if ((ang_it == 0) || !ReUseParents) {
        NNuPDGTargets[nuPDG_it]++;
      }

      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {

        double ppfx_w = 1;
        if (use_PPFX) {
          if (ppfx_univ_it == 1) {
            ppfx_w = dk2nuRdr.ppfx_cvwgt;
          } else if (ppfx_univ_it > 1) {
            ppfx_w = dk2nuRdr.ppfx_vwgt_tot[ppfx_univ_it - 2];
          }
        }

        // Should, but doesn't seem to work when filling TH2Ds
        Hists[ppfx_univ_it][nuPDG_it][0]->AddBinContent(bin_it, w * ppfx_w);
        if (DoExtra) {
          Hists[ppfx_univ_it][nuPDG_it][extra_hist_index]->AddBinContent(
              bin_it, w * ppfx_w);
        }

      } // End loop over PPFX universes

      // If we aren't re-using the parents then we have placed this neutrino
      // randomly in the 2D range and should now move to the next one.
      if (!ReUseParents) {
        break;
      }

      // If we only have a single position, then we must use this to break out
      // of the loop.
      if (!NOffAxisBins) {
        break;
      }
    } // End loop over angle bins
  }   // End loop over decay parents

  for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
    std::cout << "[INFO]: Found " << NNuPDGTargets[nuPDG_it]
              << " Neutrinos with PDG: " << NuPDGTargets[nuPDG_it] << std::endl;

    double integ = Hists[0][nuPDG_it][0]->Integral();
    if (!std::isnormal(integ)) {
      std::cerr << "[ERROR]: Flux for PDG: " << NuPDGTargets[nuPDG_it]
                << " has bad integral (" << integ << ")" << std::endl;
      throw;
    }

    Int_t NOffAxisBins = Int_t(OffAxisSteps[nuPDG_it].size()) - 1;
    // Scale to /cm^2 and per GeV (xbin width)
    if (NOffAxisBins > 1) {
      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {
        std::vector<TH2 *> &univ_hists = Hists_2D[ppfx_univ_it][nuPDG_it];
        for (size_t hist_index = 0; hist_index < (DoExtra ? 5 : 1);
             ++hist_index) {
          TH2 *hist = univ_hists[hist_index];
          for (Int_t xbin_it = 0; xbin_it < hist->GetXaxis()->GetNbins();
               ++xbin_it) {
            double scale = 1E-4 / hist->GetXaxis()->GetBinWidth(xbin_it + 1);
            for (Int_t ybin_it = 0; ybin_it < hist->GetYaxis()->GetNbins();
                 ++ybin_it) {

              double bc = hist->GetBinContent(xbin_it + 1, ybin_it + 1);
              double be = hist->GetBinError(xbin_it + 1, ybin_it + 1);
              hist->SetBinContent(xbin_it + 1, ybin_it + 1, bc * scale);
              hist->SetBinError(xbin_it + 1, ybin_it + 1, be * scale);
            } // Loop over y bins
          }   // Loop over x bins
        }     // Loop over decay parents
      }       // Loop over universes
    } else {
      for (size_t ppfx_univ_it = 0; ppfx_univ_it < NPPFXU; ++ppfx_univ_it) {
        for (size_t hist_index = 0; hist_index < (DoExtra ? 5 : 1);
             ++hist_index) {
          Hists[ppfx_univ_it][nuPDG_it][hist_index]->Scale(1E-4, "width");
        }
      }
    }
  } // Loop over neutrino species

  outfile->Write();
  outfile->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!EnergyBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "40,0,10"};
    handleOpts(argc_dum, argv_dum);
  }

  if (!OffAxisSteps.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-x", "0"};
    handleOpts(argc_dum, argv_dum);
  }

  DK2NuReader *dk2nuRdr = new DK2NuReader(
      DK2NULite ? "dk2nuTree_lite" : "dk2nuTree", inpDir, DK2NULite);

  if (ReadDK2NULitePPFX) {
    dk2nuRdr->SetPPFXBranchAddresses(NPPFX_Universes);
  }

  if (!dk2nuRdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }

  DKMetaReader *dkmRdr = new DKMetaReader(
      DK2NULite ? "dkmetaTree_lite" : "dkmetaTree", inpDir, DK2NULite);

  int metaNEntries = dkmRdr->GetEntries();

  double TotalPOT = 0;
  for (int i = 0; i < metaNEntries; ++i) {
    dkmRdr->GetEntry(i);
    TotalPOT += dkmRdr->pots;
  }

  std::cout << "Total POT: " << TotalPOT << std::endl;
  AllInOneGo(*dk2nuRdr, TotalPOT);
}
