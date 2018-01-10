#include "DetectorStop.hxx"
#include "NuDataReader.hxx"
#include "Utils.hxx"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

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
  static const double detRadius = 100.0;  // in cm
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

  if (!isnormal(nu_energy) || !isnormal(nu_wght)) {
    std::cout << "[ERROR]: Calculated bad nustats." << std::endl;
    exit(1);
  }
  return std::make_tuple(nu_energy, nuRay.Theta(), nu_wght);
}

std::vector<double> mrads;
std::vector<double> latsteps;
std::string inpDir = ".";
bool DoExtra = false;

int NBins = 50;
double BLow = 0;
double BUp = 10;
double detector_width = 0;
double detector_height = 0;
bool ReUseParents = true;
bool DK2NULite = false;

double ZDist = 57400;
int NMaxNeutrinos = -1;

std::vector<double> varBin;
bool VariableBinning = false;

bool FillStopNeutrinoTree = false;
bool DoDivergence = false;

std::string outputFile;

std::string runPlanCfg, runPlanName = "";

int OnlySpecies = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << "\t-i /path/to/DUNE/nudata/files : Can include wildcards "
               "(remeber to quote \n"
               "\t                                to avoid shell expansion.)   "
               "           \n"
               "\n"
               "\t-o output.root                : File to write fluxes to.     "
               "           \n"
               "\n"
               "\t-a 1,2,3[,...]                : Specific mrad off-axis "
               "fluxes to        \n"
               "\t                                calculate.                   "
               "           \n"
               "\n"
               "\t-d <step_deg>,<max_deg>       : Calculate fluxes from deg = "
               "0 to       \n"
               "\t                                deg = <max_deg> in steps of "
               "           \n"
               "\t                                deg = <step_deg>.           "
               "           \n"
               "\n"
               "\t-x <step_x>,<max_x>           : Calculate fluxes from lateral"
               " offset =  \n"
               "\t                                x = <max_x> in steps of "
               "           \n"
               "\t                                x = <step_x>. x in cm.       "
               "           \n"
               "\n"
               "\t-r <RunPlan.XML>              : An XML file specifying a run "
               "plan  \n"
               "\t                                to build fluxes for. See     "
               "      \n"
               "\t                                documentation for XML "
               "structure.   \n"
               "\n"
               "\t-e                            : Build fluxes for specific "
               "neutrino decay\n"
               "\t                                parents.                     "
               "           \n"
               "\n"
               "\t-b <NBins>,<Low>,<High>       : Use uniform binning for flux "
               "histograms.\n"
               "\n"
               "\t-vb <bin0low>,<bin1low>_<binXUp>:step,..,<binYlow>,<binYup> "
               ": Use       \n"
               "\t                                variable binning specified "
               "by bin edges \n"
               "\t                                and step ranges.             "
               "           \n"
               "\n"
               "\t-w <DetectorWidth=0           : Width of detector plane (cm)."
               "\n\n"
               "\t-h <DetectorHeight=0>         : Height of detector plane (cm)"
               ".\n"
               "\n"
               "\t-n <NMaxNeutrinos>            : Only loop over -n nus.    \n"
               "\n"
               "\t-z <ZDist>                    : Detector ZDist (cm).       \n"
               "\n"
               "\t-P                            : Only use each decaying parent"
               " once.       \n"
               "\t-F                            : Fill tree of neutrino "
               "energies that can be used to \n"
               "\t                                rebin flux.\n"
               "\t-S <species PDG>              : Only fill information for "
               "neutrinos of a given species.\n"
               "\t-D                            : Calculate flux divergences.\n"
               "\t-L                            : Expect dk2nulite inputs.\n"
               "\n"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-a") {
      mrads = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-e") {
      DoExtra = true;
    } else if (std::string(argv[opt]) == "-b") {
      std::vector<double> binning = ParseToVect<double>(argv[++opt], ",");
      if (binning.size() != 3) {
        std::cout << "[ERROR]: Recieved " << binning.size()
                  << " entrys for -b, expected 3." << std::endl;
        exit(1);
      }
      NBins = int(binning[0]);
      BLow = binning[1];
      BUp = binning[2];
    } else if (std::string(argv[opt]) == "-vb") {
      std::vector<std::string> vbDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");
      varBin.clear();
      for (size_t vbd_it = 0; vbd_it < vbDescriptors.size(); ++vbd_it) {
        AppendVect(varBin, BuildDoubleList(vbDescriptors[vbd_it]));
      }

      for (size_t bin_it = 1; bin_it < varBin.size(); ++bin_it) {
        if (varBin[bin_it] == varBin[bin_it - 1]) {
          std::cout << "[INFO]: Removing duplciate bin low edge " << bin_it
                    << " low edge: " << varBin[bin_it] << std::endl;
          varBin.erase(varBin.begin() + bin_it);
        }
      }

      for (size_t bin_it = 1; bin_it < varBin.size(); ++bin_it) {
        if (varBin[bin_it] < varBin[bin_it - 1]) {
          std::cout << "[ERROR]: Bin " << bin_it
                    << " low edge: " << varBin[bin_it]
                    << " is smaller than bin " << (bin_it - 1)
                    << " low edge: " << varBin[bin_it - 1] << std::endl;
          exit(1);
        }
      }
      VariableBinning = true;
    } else if (std::string(argv[opt]) == "-d") {
      std::vector<double> degstep = ParseToVect<double>(argv[++opt], ",");
      if (degstep.size() != 2) {
        std::cout << "[ERROR]: Recieved " << degstep.size()
                  << " entrys for -d, expected 2." << std::endl;
        exit(1);
      }
      mrads.push_back(0);
      static const double deg2rad = asin(1) / 90.0;
      double curr = degstep[0];
      while (curr < degstep[1]) {
        mrads.push_back(curr * deg2rad * 1E3);
        curr += degstep[0];
      }

    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
    } else if (std::string(argv[opt]) == "-x") {
      std::vector<std::string> latDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");
      latsteps.clear();
      for (size_t vbd_it = 0; vbd_it < latDescriptors.size(); ++vbd_it) {
        AppendVect(latsteps, BuildDoubleList(latDescriptors[vbd_it]));
      }
    } else if (std::string(argv[opt]) == "-w") {
      detector_width = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-h") {
      detector_height = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-n") {
      NMaxNeutrinos = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-z") {
      ZDist = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-P") {
      ReUseParents = false;
    } else if (std::string(argv[opt]) == "-F") {
      FillStopNeutrinoTree = true;
    } else if (std::string(argv[opt]) == "-S") {
      OnlySpecies = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-D") {
      DoDivergence = true;
    } else if (std::string(argv[opt]) == "-L") {
      DK2NULite = true;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

constexpr double rad2deg = 90.0 / asin(1);

void AllInOneGo(DK2NuReader &dk2nuRdr, double TotalPOT) {
  TRandom3 rnjesus;

  size_t NAngs = mrads.size();
  size_t nNus = (NMaxNeutrinos == -1)
                    ? dk2nuRdr.GetEntries()
                    : std::min(NMaxNeutrinos, int(dk2nuRdr.GetEntries()));

  TotalPOT = TotalPOT * (double(nNus) / double(dk2nuRdr.GetEntries()));
  std::cout << "Only using the first " << nNus << " events out of "
            << dk2nuRdr.GetEntries() << ", scaling POT to " << TotalPOT
            << std::endl;
  std::cout << "Reding " << nNus << " Dk2Nu entries." << std::endl;

  TFile *outfile = new TFile(outputFile.c_str(), "RECREATE");

  std::vector<std::vector<std::vector<TH1D *> > > Hists;
  Hists.resize(NAngs);

  std::vector<int> NuPDGTargets = {-12, 12, -14, 14};
  if (OnlySpecies) {
    NuPDGTargets = std::vector<int>{OnlySpecies};
  }
  std::vector<std::string> NuNames = {"nueb", "nue", "numub", "numu"};

  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    Hists.push_back(std::vector<std::vector<TH1D *> >());
    for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
      Hists[ang_it].push_back(std::vector<TH1D *>());
    }
  }
  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
      std::stringstream ss("");
      if (!latsteps.size()) {
        ss << "LBNF_" << NuNames[nuPDG_it] << "_mrad_" << mrads[ang_it];
      } else {
        ss << "LBNF_" << NuNames[nuPDG_it] << "_lateral_displace_"
           << latsteps[ang_it] << "_cm";
      }

      Hists[ang_it][nuPDG_it].push_back(
          VariableBinning
              ? new TH1D(
                    ss.str().c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    varBin.size() - 1, varBin.data())
              : new TH1D(
                    ss.str().c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    NBins, BLow, BUp));

      if (DoExtra) {
        Hists[ang_it][nuPDG_it].push_back(
            VariableBinning ? new TH1D((ss.str() + "_pi").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       varBin.size() - 1, varBin.data())
                            : new TH1D((ss.str() + "_pi").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       NBins, BLow, BUp));
        Hists[ang_it][nuPDG_it].push_back(
            VariableBinning ? new TH1D((ss.str() + "_k").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       varBin.size() - 1, varBin.data())
                            : new TH1D((ss.str() + "_k").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       NBins, BLow, BUp));
        Hists[ang_it][nuPDG_it].push_back(
            VariableBinning ? new TH1D((ss.str() + "_k0").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       varBin.size() - 1, varBin.data())
                            : new TH1D((ss.str() + "_k0").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       NBins, BLow, BUp));
        Hists[ang_it][nuPDG_it].push_back(
            VariableBinning ? new TH1D((ss.str() + "_mu").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       varBin.size() - 1, varBin.data())
                            : new TH1D((ss.str() + "_mu").c_str(),
                                       ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                       "(GeV^{-1}cm^{-2} per POT)",
                                       NBins, BLow, BUp));
      }
    }
  }

  std::vector<TVector3> detPoses;
  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    detPoses.push_back(TVector3(ZDist * tan(mrads[ang_it] * 1.0E-3), 0, ZDist));

    std::cout << "[INFO]: Building flux at: "
              << (mrads[ang_it] * 1.0E-3 * rad2deg) << " degrees ("
              << mrads[ang_it] << " mrads) Offset = "
              << ((ZDist * 1.0E-2) * tan(mrads[ang_it] * 1.0E-3)) << " m."
              << std::endl;
  }
  size_t updateStep = nNus / 10 ? nNus / 10 : 1;
  for (size_t nu_it = 0; nu_it < nNus; ++nu_it) {
    if (!(nu_it % updateStep)) {
      std::cout << "--" << nu_it << "/" << nNus << std::endl;
    }

    dk2nuRdr.GetEntry(nu_it);

    if (OnlySpecies && (OnlySpecies != dk2nuRdr.decay_ntype)) {
      continue;
    }

    double wF = (dk2nuRdr.decay_nimpwt / TMath::Pi()) * (1.0 / TotalPOT);
    for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
      TVector3 det_point = detPoses[ang_it];

      if (detector_width != 0) {
        det_point[0] += (2.0 * rnjesus.Uniform() - 1.0) * detector_width;
      }

      if (detector_height != 0) {
        det_point[1] += (2.0 * rnjesus.Uniform() - 1.0) * detector_height;
      }

      std::tuple<double, double, double> nuStats =
          GetNuWeight(dk2nuRdr, det_point);

      double w = std::get<2>(nuStats) * wF;

      int nuPDG_it =
          std::distance(NuPDGTargets.begin(),
                        std::find(NuPDGTargets.begin(), NuPDGTargets.end(),
                                  dk2nuRdr.decay_ntype));
      if (nuPDG_it == 4) {
        std::cout << "Warning, couldn't find plot index for NuPDG: "
                  << dk2nuRdr.decay_ntype << std::endl;
        exit(1);
      }

      Hists[ang_it][nuPDG_it][0]->Fill(std::get<0>(nuStats), w);
      std::cout << "[VERBOSE]: " << ang_it << ", " << nuPDG_it << ", "
                << std::get<0>(nuStats) << ", " << w << std::endl;

      if ((!isnormal(std::get<0>(nuStats)) || (!isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << "*" << TotalPOT << ")."
                  << std::endl;
        throw;
      }

      if (DoExtra) {
        if ((dk2nuRdr.decay_ptype == 211) || (dk2nuRdr.decay_ptype == -211)) {
          Hists[ang_it][nuPDG_it][1]->Fill(std::get<0>(nuStats), w);
        } else if ((dk2nuRdr.decay_ptype == 321) ||
                   (dk2nuRdr.decay_ptype == -321)) {
          Hists[ang_it][nuPDG_it][2]->Fill(std::get<0>(nuStats), w);
        } else if ((dk2nuRdr.decay_ptype == 311) ||
                   (dk2nuRdr.decay_ptype == 310) ||
                   (dk2nuRdr.decay_ptype == 130)) {
          Hists[ang_it][nuPDG_it][3]->Fill(std::get<0>(nuStats), w);
        } else if ((dk2nuRdr.decay_ptype == 13) ||
                   (dk2nuRdr.decay_ptype == -13)) {
          Hists[ang_it][nuPDG_it][4]->Fill(std::get<0>(nuStats), w);
        }
      }
    }
  }

  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
      double integ = Hists[ang_it][nuPDG_it][0]->Integral();
      if (!std::isnormal(integ)) {
        std::cerr << "[ERROR]: Flux @ angular point: " << ang_it
                  << " for PDG: " << NuPDGTargets[nuPDG_it]
                  << " has bad integral." << std::endl;
        throw;
      }

      Hists[ang_it][nuPDG_it][0]->Scale(1E-4, "width");
      if (DoExtra) {
        Hists[ang_it][nuPDG_it][1]->Scale(1E-4, "width");
        Hists[ang_it][nuPDG_it][2]->Scale(1E-4, "width");
        Hists[ang_it][nuPDG_it][3]->Scale(1E-4, "width");
        Hists[ang_it][nuPDG_it][4]->Scale(1E-4, "width");
      }
    }
  }

  outfile->Write();
  outfile->Close();
}

void CalculateFluxesForRunPlan(DK2NuReader &dk2nuRdr, double TotalPOT,
                               std::vector<DetectorStop> &detStops) {
  TRandom3 rnjesus;
  size_t nNus = (NMaxNeutrinos == -1)
                    ? dk2nuRdr.GetEntries()
                    : std::min(NMaxNeutrinos, int(dk2nuRdr.GetEntries()));
  TotalPOT = TotalPOT * (double(nNus) / double(dk2nuRdr.GetEntries()));
  std::cout << "Only using the first " << nNus << " events out of "
            << dk2nuRdr.GetEntries() << ", scaling totla POT to " << TotalPOT
            << std::endl;
  std::cout << "Reading " << nNus << " Dk2Nu entries." << std::endl;

  TFile *outfile = new TFile(outputFile.c_str(), "RECREATE");

  std::vector<double> MeasurementsOffsets;
  std::vector<double> MeasurementHalfWidths;
  std::vector<double> MeasurementHalfHeights;
  std::vector<TVector3> detPositions;
  std::vector<std::pair<size_t, size_t> > SliceToStopSliceMap;

  std::vector<std::map<int, TH1D *> > Fluxes;
  std::vector<std::map<int, TH2D *> > Divergences;

  std::stringstream ss("");

  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    DetectorStop &ds = detStops[ds_it];

    for (size_t ms_it = 0; ms_it < ds.GetNMeasurementSlices(); ++ms_it) {
      MeasurementsOffsets.push_back(ds.GetAbsoluteOffsetOfSlice(ms_it));

      std::cout << "[INFO]: Building flux for measurement slice: Stop " << ds_it
                << ", slice: " << ms_it
                << " at absolute offset: " << MeasurementsOffsets.back()
                << " m." << std::endl;

      MeasurementHalfWidths.push_back(ds.MeasurementRegionWidth / 2.0);
      MeasurementHalfHeights.push_back(ds.DetectorFiducialHeight / 2.0);

      detPositions.push_back(
          TVector3(MeasurementsOffsets.back() * 100.0, 0, ZDist));
      SliceToStopSliceMap.push_back(std::make_pair(ds_it, ms_it));

      std::map<int, TH1D *> fluxhistos;
      std::map<int, TH2D *> divhistos;

      for (int species : {-14, -12, 12, 14}) {
        if (OnlySpecies && (OnlySpecies != species)) {
          continue;
        }

        ss.str("");
        ss << GetSpeciesName(species) << "_flux_stop_" << ds_it << "_slice_"
           << ms_it << "_" << MeasurementsOffsets.back() << "_m";

        fluxhistos[species] = VariableBinning
                                  ? new TH1D(ss.str().c_str(),
                                             ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                             "(GeV^{-1}cm^{-2} per POT)",
                                             varBin.size() - 1, varBin.data())
                                  : new TH1D(ss.str().c_str(),
                                             ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
                                             "(GeV^{-1}cm^{-2} per POT)",
                                             NBins, BLow, BUp);
        fluxhistos[species]->SetDirectory(nullptr);

        if (DoDivergence) {
          ss.str("");
          ss << GetSpeciesName(species) << "_divergence_stop_" << ds_it
             << "_slice_" << ms_it << "_" << MeasurementsOffsets.back() << "_m";
          divhistos[species] =
              new TH2D(ss.str().c_str(),
                       ";#it{E}_{#nu} (GeV);#theta_{#nu} (degrees);#Phi_{#nu} "
                       "per POT",
                       80, BLow, BUp, 240, 0, 6);
          divhistos[species]->SetDirectory(nullptr);
        }
      }

      Fluxes.push_back(fluxhistos);
      if (DoDivergence) {
        Divergences.push_back(divhistos);
      }
    }
  }

  double FluxNormPOT = TotalPOT;
  if (!ReUseParents) {
    FluxNormPOT /= Fluxes.size();
    std::cout
        << "[INFO]: Each decay parent may only be used once -- POT per slice: "
        << FluxNormPOT << std::endl;
  }

  size_t NSlices = detPositions.size();
  size_t updateStep = nNus / 10 ? nNus / 10 : 1;
  size_t slice_it = 0;
  for (size_t nu_it = 0; nu_it < nNus; ++nu_it) {
    if (!(nu_it % updateStep)) {
      std::cout << "--" << nu_it << "/" << nNus << std::endl;
    }

    dk2nuRdr.GetEntry(nu_it);

    if (OnlySpecies && (OnlySpecies != dk2nuRdr.decay_ntype)) {
      continue;
    }

    double wF = (dk2nuRdr.decay_nimpwt / TMath::Pi()) * (1.0 / FluxNormPOT);

    if (ReUseParents) {
      for (slice_it = 0; slice_it < NSlices; ++slice_it) {
        TVector3 det_point = detPositions[slice_it];

        det_point[0] += (2.0 * rnjesus.Uniform() - 1.0) *
                        MeasurementHalfWidths[slice_it] * 100.0;

        det_point[1] += (2.0 * rnjesus.Uniform() - 1.0) *
                        MeasurementHalfHeights[slice_it] * 100.0;

        std::tuple<double, double, double> nuStats =
            GetNuWeight(dk2nuRdr, det_point);

        double w = std::get<2>(nuStats) * wF;

        if ((std::get<0>(nuStats) != std::get<0>(nuStats)) || (w != w)) {
          std::cout << std::get<0>(nuStats) << ", " << w << std::endl;
          throw;
        }

        Fluxes[slice_it][dk2nuRdr.decay_ntype]->Fill(std::get<0>(nuStats), w);
        if (DoDivergence) {
          Divergences[slice_it][dk2nuRdr.decay_ntype]->Fill(
              std::get<0>(nuStats),
              std::get<1>(nuStats) * (180.0 / TMath::Pi()), w);
        }

        if (FillStopNeutrinoTree) {
          std::pair<size_t, size_t> dsms = SliceToStopSliceMap[slice_it];

          detStops[dsms.first].FillNeutrino(dsms.second, dk2nuRdr.decay_ntype,
                                            std::get<0>(nuStats), w);
        }
      }
    } else {
      slice_it = (slice_it + 1) % NSlices;

      TVector3 det_point = detPositions[slice_it];

      det_point[0] += (2.0 * rnjesus.Uniform() - 1.0) *
                      MeasurementHalfWidths[slice_it] * 100.0;

      det_point[1] += (2.0 * rnjesus.Uniform() - 1.0) *
                      MeasurementHalfHeights[slice_it] * 100.0;

      std::tuple<double, double, double> nuStats =
          GetNuWeight(dk2nuRdr, det_point);

      double w = std::get<2>(nuStats) * wF;

      if ((!isnormal(std::get<0>(nuStats)) || (!isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << "*" << TotalPOT << ")."
                  << std::endl;
        throw;
      }

      Fluxes[slice_it][dk2nuRdr.decay_ntype]->Fill(std::get<0>(nuStats), w);
      if (DoDivergence) {
        Divergences[slice_it][dk2nuRdr.decay_ntype]->Fill(
            std::get<0>(nuStats), std::get<1>(nuStats) * (180.0 / TMath::Pi()),
            w);
      }

      if (FillStopNeutrinoTree) {
        std::pair<size_t, size_t> dsms = SliceToStopSliceMap[slice_it];

        detStops[dsms.first].FillNeutrino(dsms.second, dk2nuRdr.decay_ntype,
                                          std::get<0>(nuStats), w);
      }
    }
  }

  slice_it = 0;
  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    DetectorStop &ds = detStops[ds_it];
    for (size_t ms_it = 0; ms_it < ds.GetNMeasurementSlices(); ++ms_it) {
      for (int species : {-14, -12, 12, 14}) {
        if (OnlySpecies && (OnlySpecies != species)) {
          continue;
        }
        Fluxes[slice_it][species]->Scale(1E-4, "width");
        ds.AddSliceFlux(ms_it, species, Fluxes[slice_it][species]);
        if (DoDivergence) {
          ds.AddSliceDivergence(ms_it, species, Divergences[slice_it][species]);
        }
      }
      slice_it++;
    }
    ds.Write();
  }

  std::map<int, TGraph> FluxNorms = {{-14, TGraph(Fluxes.size())},
                                     {14, TGraph(Fluxes.size())},
                                     {-12, TGraph(Fluxes.size())},
                                     {12, TGraph(Fluxes.size())}};

  slice_it = 0;
  for (size_t ds_it = 0; ds_it < detStops.size(); ++ds_it) {
    DetectorStop &ds = detStops[ds_it];
    for (size_t ms_it = 0; ms_it < ds.GetNMeasurementSlices(); ++ms_it) {
      for (int species : {-14, -12, 12, 14}) {
        if (OnlySpecies && (OnlySpecies != species)) {
          continue;
        }
        double sliceOffset = ds.GetAbsoluteOffsetOfSlice(ms_it);
        double theta_rad = atan(sliceOffset / ZDist);
        FluxNorms[species].SetPoint(
            slice_it, theta_rad * rad2deg,
            ds.GetFluxForSpecies(ms_it, species)->Integral("width"));
      }
      slice_it++;
    }
  }

  for (int species : {-14, -12, 12, 14}) {
    if (OnlySpecies && (OnlySpecies != species)) {
      continue;
    }
    FluxNorms[species].GetXaxis()->SetTitle("Off-axis angle (Degrees)");
    FluxNorms[species].GetYaxis()->SetTitle(
        "#Phi_{#nu}^{Total} per m^{2} per POT");
    ss.str("");
    ss << GetSpeciesName(species) << "_norm";
    FluxNorms[species].Write(ss.str().c_str(), TObject::kOverwrite);
  }

  outfile->Write();
  outfile->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if (!mrads.size() && latsteps.size()) {
    for (auto ls : latsteps) {
      double theta = atan(ls / ZDist);
      mrads.push_back(theta * 1.0E3);

      std::cout << "[INFO]: Building flux [" << mrads.size()
                << "] at lateral position: " << ls
                << " (mrad = " << mrads.back() << ")." << std::endl;
    }
  }

  if (!mrads.size()) {
    mrads.push_back(0);
  }

  DK2NuReader *dk2nuRdr = new DK2NuReader(
      DK2NULite ? "dk2nuTree_lite" : "dk2nuTree", inpDir, DK2NULite);

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

  if (runPlanCfg.length()) {
    std::vector<DetectorStop> detStops =
        ReadDetectorStopConfig(runPlanCfg, runPlanName);
    CalculateFluxesForRunPlan(*dk2nuRdr, TotalPOT, detStops);
  } else {
    AllInOneGo(*dk2nuRdr, TotalPOT);
  }
}
