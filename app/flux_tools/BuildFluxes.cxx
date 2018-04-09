#include "dk2nu_TreeReader.hxx"

#include "ROOTUtility.hxx"
#include "PhysicsUtility.hxx"
#include "GetUsage.hxx"

#include "TChain.h"
#include "TFile.h"
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

  if (!std::isnormal(nu_energy) || !std::isnormal(nu_wght)) {
    std::cout << "[ERROR]: Calculated bad nustats." << std::endl;
    exit(1);
  }
  return std::make_tuple(nu_energy, nuRay.Theta(), nu_wght);
}

std::vector<double> OffAxisSteps;
enum OffAxisStepUnits { kPostion_m = 0, kmrad, kdegrees };
OffAxisStepUnits stepType;

std::string inpDir = ".";
bool DoExtra = false;

double detector_half_height = 0;
bool ReUseParents = true;
bool DK2NULite = false;

double ZDist = 57400;
int NMaxNeutrinos = -1;

std::vector<double> EnergyBinning;

// bool DoDivergence = false;

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
    if (std::string(argv[opt]) == "-?"||
               std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-m") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
        << " " << argv[opt+1] << std::endl;
      OffAxisSteps = BuildBinEdges(argv[++opt]);
      stepType = kmrad;
    } else if (std::string(argv[opt]) == "-d") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
        << " " << argv[opt+1] << std::endl;
      OffAxisSteps = BuildBinEdges(argv[++opt]);
      stepType = kdegrees;
    }  else if (std::string(argv[opt]) == "-x") {
      std::cout << "\t--Handling off-axis position specifier: " << argv[opt]
        << " " << argv[opt+1] << std::endl;
      OffAxisSteps = BuildBinEdges(argv[++opt]);
      stepType = kPostion_m;
    } else if (std::string(argv[opt]) == "-e") {
      DoExtra = true;
    } else if (std::string(argv[opt]) == "-b") {
      std::vector<double> binningDescriptor = ParseToVect<double>(argv[++opt], ",");
      if (binningDescriptor.size() != 3) {
        std::cout << "[ERROR]: Recieved " << binningDescriptor.size()
                  << " entrys for -b, expected 3." << std::endl;
        exit(1);
      }
      int NBins = int(binningDescriptor[0]);
      double BLow = binningDescriptor[1];
      double BUp = binningDescriptor[2];
      double bwidth = (BUp - BLow)/double(NBins);
      EnergyBinning.push_back(BLow);
      for(Int_t i = 0; i < NBins; ++i){
        EnergyBinning.push_back(EnergyBinning.back() + bwidth);
      }
    } else if (std::string(argv[opt]) == "-vb") {
      EnergyBinning = BuildBinEdges(argv[++opt]);
    } else if (std::string(argv[opt]) == "-h") {
      detector_half_height = str2T<double>(argv[++opt])/2.0;
    } else if (std::string(argv[opt]) == "-n") {
      NMaxNeutrinos = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-z") {
      ZDist = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-P") {
      ReUseParents = false;
    } else if (std::string(argv[opt]) == "-S") {
      OnlySpecies = str2T<int>(argv[++opt]);
    // } else if (std::string(argv[opt]) == "-D") {
    //   DoDivergence = true;
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

std::pair<double, TVector3> GetRandomFluxWindowPosition(TRandom3 &rnjesus,
  Int_t FluxWindow=-1){
  double OffAxisCenter = ((FluxWindow == -1) ?
    (OffAxisSteps.back()+OffAxisSteps.front()) :
    (OffAxisSteps[FluxWindow+1] + OffAxisSteps[FluxWindow])) / 2.0;
  double OffAxisHalfRange = ((FluxWindow == -1) ?
    (OffAxisSteps.back()-OffAxisSteps.front()) :
    (OffAxisSteps[FluxWindow+1] - OffAxisSteps[FluxWindow])) / 2.0;

    // std::cout << "[INFO]: Flux window[" << FluxWindow << "/" << OffAxisSteps.size() << "] = {"
    // << OffAxisCenter << " \\pm " << OffAxisHalfRange << "}" << std::endl;

    TVector3 rndDetPos(0, (2.0 * rnjesus.Uniform() - 1.0) * detector_half_height, ZDist);

    double RandomOffAxisPos = OffAxisCenter + (2.0 * rnjesus.Uniform() - 1.0) * OffAxisHalfRange;
    switch(stepType){
      case kPostion_m:{
        rndDetPos[0] = RandomOffAxisPos*100.0; // to cm
        break;
      }
      case kmrad:{
        rndDetPos[0] = tan(RandomOffAxisPos*1E-3) * ZDist;
        break;
      }
      case kdegrees:{
        static const double deg2rad = asin(1) / 90.0;
        rndDetPos[0] = tan(RandomOffAxisPos*deg2rad) * ZDist;
        break;
      }
    }

    // std::cout << "\t Threw = { " << rndDetPos[0] << ", "
    //   << rndDetPos[1] << ", " << rndDetPos[2] << "}." << std::endl;

    return std::make_pair(RandomOffAxisPos, rndDetPos);
}

void AllInOneGo(DK2NuReader &dk2nuRdr, double TotalPOT) {
  TRandom3 rnjesus;

  Int_t NOffAxisBins = Int_t(OffAxisSteps.size()) - 1;
  if(NOffAxisBins == -1){
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
    NuPDGTargets = std::vector<int>{OnlySpecies};
  }

  std::vector< std::vector<TH1 *> > Hists;
  std::vector< std::vector<TH2D *> > Hists_2D;
  Hists.resize(NuPDGTargets.size());
  Hists_2D.resize(NuPDGTargets.size());

  for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
    std::stringstream hist_name("");
    hist_name << "LBNF_" << GetSpeciesName(NuPDGTargets[nuPDG_it]) << "_flux";

    std::string hist_title = ";#it{E}_{#nu} (GeV);#Phi_{#nu} "
    "(cm^{-2} per POT per 1 GeV)";
    if(NOffAxisBins > 1){
      switch(stepType){
        case kPostion_m:{
          hist_title = ";#it{E}_{#nu} (GeV);Off-axis postion (m);#Phi_{#nu} "
          "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        case kmrad:{
          hist_title = ";#it{E}_{#nu} (GeV);Off-axis angle (mrad);#Phi_{#nu} "
          "(cm^{-2} per POT per 1 GeV)";
          break;
        }
        case kdegrees:{
          hist_title = ";#it{E}_{#nu} (GeV);Off-axis angle (degrees);#Phi_{#nu} "
          "(cm^{-2} per POT per 1 GeV)";
          break;
        }
      }
    }

    Hists[nuPDG_it].push_back(
        (NOffAxisBins < 2) ? new TH1D(hist_name.str().c_str(),
                  hist_title.c_str(),
                  (EnergyBinning.size() - 1), EnergyBinning.data())
                  : static_cast<TH1 *>(new TH2D(hist_name.str().c_str(),
                  hist_title.c_str(),
                  (EnergyBinning.size() - 1), EnergyBinning.data(),
                  (OffAxisSteps.size() - 1), OffAxisSteps.data()))
                );
    if(NOffAxisBins > 1){
      Hists_2D[nuPDG_it].push_back(static_cast<TH2D *>(Hists[nuPDG_it].back()));
    }

    if (DoExtra) {
      Hists[nuPDG_it].push_back(
          (NOffAxisBins < 2) ? new TH1D((hist_name.str() + "_pi").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data())
                    : static_cast<TH1 *>(new TH2D((hist_name.str() + "_pi").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data(),
                    (OffAxisSteps.size() - 1), OffAxisSteps.data()))
                  );
      if(NOffAxisBins > 1){
        Hists_2D[nuPDG_it].push_back(static_cast<TH2D *>(Hists[nuPDG_it].back()));
      }
      Hists[nuPDG_it].push_back(
          (NOffAxisBins < 2) ? new TH1D((hist_name.str() + "_k").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data())
                    : static_cast<TH1 *>(new TH2D((hist_name.str() + "_k").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data(),
                    (OffAxisSteps.size() - 1), OffAxisSteps.data()))
                  );
      if(NOffAxisBins > 1){
        Hists_2D[nuPDG_it].push_back(static_cast<TH2D *>(Hists[nuPDG_it].back()));
      }
      Hists[nuPDG_it].push_back(
          (NOffAxisBins < 2) ? new TH1D((hist_name.str() + "_k0").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data())
                    : static_cast<TH1 *>(new TH2D((hist_name.str() + "_k0").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data(),
                    (OffAxisSteps.size() - 1), OffAxisSteps.data()))
                  );
      if(NOffAxisBins > 1){
        Hists_2D[nuPDG_it].push_back(static_cast<TH2D *>(Hists[nuPDG_it].back()));
      }
      Hists[nuPDG_it].push_back(
          (NOffAxisBins < 2) ? new TH1D((hist_name.str() + "_mu").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data())
                    : static_cast<TH1 *>(new TH2D((hist_name.str() + "_mu").c_str(),
                    hist_title.c_str(),
                    (EnergyBinning.size() - 1), EnergyBinning.data(),
                    (OffAxisSteps.size() - 1), OffAxisSteps.data()))
                  );
      if(NOffAxisBins > 1){
        Hists_2D[nuPDG_it].push_back(static_cast<TH2D *>(Hists[nuPDG_it].back()));
      }
    }
  }

  size_t updateStep = NNeutrinos / 10 ? NNeutrinos / 10 : 1;
  for (size_t nu_it = 0; nu_it < NNeutrinos; ++nu_it) {
    if (!(nu_it % updateStep)) {
      std::cout << "--" << nu_it << "/" << NNeutrinos << std::endl;
    }

    dk2nuRdr.GetEntry(nu_it);

    if (OnlySpecies && (OnlySpecies != dk2nuRdr.decay_ntype)) {
      continue;
    }

    double wF = (dk2nuRdr.decay_nimpwt / TMath::Pi()) * (1.0 / TotalPOT);

    for (Int_t ang_it = 0; (!NOffAxisBins) || (ang_it < NOffAxisBins);
      ++ang_it) {

      //If we are not re-using the decay parents, then this is placed randomly
      //over the whole off-axis range.
      std::pair<double,TVector3> det_point =
        GetRandomFluxWindowPosition(rnjesus,ReUseParents?ang_it:-1);

      std::tuple<double, double, double> nuStats =
          GetNuWeight(dk2nuRdr, det_point.second);

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

      if(NOffAxisBins > 1){
        Hists_2D[nuPDG_it][0]->Fill(std::get<0>(nuStats), det_point.first, w);
      } else {
        Hists[nuPDG_it][0]->Fill(std::get<0>(nuStats), w);
      }

      if ((!std::isnormal(std::get<0>(nuStats)) || (!std::isnormal(w)))) {
        std::cout << std::get<0>(nuStats) << ", " << w << "("
                  << std::get<2>(nuStats) << "*" << dk2nuRdr.decay_nimpwt
                  << "/(" << TMath::Pi() << "*" << TotalPOT << ")."
                  << std::endl;
        throw;
      }

      if (DoExtra) {
        if ((dk2nuRdr.decay_ptype == 211) || (dk2nuRdr.decay_ptype == -211)) {
          if(NOffAxisBins > 1){
            Hists_2D[nuPDG_it][1]->Fill(std::get<0>(nuStats), det_point.first, w);
          } else {
            Hists[nuPDG_it][1]->Fill(std::get<0>(nuStats), w);
          }
        } else if ((dk2nuRdr.decay_ptype == 321) ||
                   (dk2nuRdr.decay_ptype == -321)) {
         if(NOffAxisBins > 1){
           Hists_2D[nuPDG_it][2]->Fill(std::get<0>(nuStats), det_point.first, w);
         } else {
           Hists[nuPDG_it][2]->Fill(std::get<0>(nuStats), w);
         }
        } else if ((dk2nuRdr.decay_ptype == 311) ||
                   (dk2nuRdr.decay_ptype == 310) ||
                   (dk2nuRdr.decay_ptype == 130)) {
         if(NOffAxisBins > 1){
           Hists_2D[nuPDG_it][3]->Fill(std::get<0>(nuStats), det_point.first, w);
         } else {
           Hists[nuPDG_it][3]->Fill(std::get<0>(nuStats), w);
         }
        } else if ((dk2nuRdr.decay_ptype == 13) ||
                   (dk2nuRdr.decay_ptype == -13)) {
         if(NOffAxisBins > 1){
           Hists_2D[nuPDG_it][4]->Fill(std::get<0>(nuStats), det_point.first, w);
         } else {
           Hists[nuPDG_it][4]->Fill(std::get<0>(nuStats), w);
         }
        }
      }
      //If we aren't re-using the parents then we have placed this neutrino
      //randomly in the 2D range and should now move to the next one.
      if(!ReUseParents){
        break;
      }

      //If we only have a single position, then we must use this to break out of
      //the loop.
      if(!NOffAxisBins){
        break;
      }
    }
  }

  for (size_t nuPDG_it = 0; nuPDG_it < NuPDGTargets.size(); ++nuPDG_it) {
    double integ = Hists[nuPDG_it][0]->Integral();
    if (!std::isnormal(integ)) {
      std::cerr << "[ERROR]: Flux for PDG: " << NuPDGTargets[nuPDG_it]
                << " has bad integral." << std::endl;
      throw;
    }

    // Scale to /cm^2 and per GeV (xbin width)
    if(NOffAxisBins > 1){
      for(Int_t xbin_it = 0; xbin_it < Hists_2D[nuPDG_it][0]->GetXaxis()->GetNbins();
        ++ xbin_it){
        for(Int_t ybin_it = 0;
          ybin_it < Hists_2D[nuPDG_it][0]->GetYaxis()->GetNbins(); ++ ybin_it){
          Hists_2D[nuPDG_it][0]->SetBinContent(xbin_it +1 ,ybin_it + 1,
            Hists_2D[nuPDG_it][0]->GetBinContent(xbin_it +1 ,ybin_it + 1) *
            1E-4 * Hists_2D[nuPDG_it][0]->GetXaxis()->GetBinWidth(xbin_it+1));
          if (DoExtra) {
            Hists_2D[nuPDG_it][1]->SetBinContent(xbin_it +1 ,ybin_it + 1,
              Hists_2D[nuPDG_it][1]->GetBinContent(xbin_it +1 ,ybin_it + 1) *
              1E-4 * Hists_2D[nuPDG_it][1]->GetXaxis()->GetBinWidth(xbin_it+1));
            Hists_2D[nuPDG_it][2]->SetBinContent(xbin_it +1 ,ybin_it + 1,
              Hists_2D[nuPDG_it][2]->GetBinContent(xbin_it +1 ,ybin_it + 1) *
              1E-4 * Hists_2D[nuPDG_it][2]->GetXaxis()->GetBinWidth(xbin_it+1));
            Hists_2D[nuPDG_it][3]->SetBinContent(xbin_it +1 ,ybin_it + 1,
              Hists_2D[nuPDG_it][3]->GetBinContent(xbin_it +1 ,ybin_it + 1) *
              1E-4 * Hists_2D[nuPDG_it][3]->GetXaxis()->GetBinWidth(xbin_it+1));
            Hists_2D[nuPDG_it][4]->SetBinContent(xbin_it +1 ,ybin_it + 1,
              Hists_2D[nuPDG_it][4]->GetBinContent(xbin_it +1 ,ybin_it + 1) *
              1E-4 * Hists_2D[nuPDG_it][4]->GetXaxis()->GetBinWidth(xbin_it+1));
          }
        }
      }
    } else {
      Hists[nuPDG_it][0]->Scale(1E-4,"width");
      if (DoExtra) {
        Hists[nuPDG_it][1]->Scale(1E-4,"width");
        Hists[nuPDG_it][2]->Scale(1E-4,"width");
        Hists[nuPDG_it][3]->Scale(1E-4,"width");
        Hists[nuPDG_it][4]->Scale(1E-4,"width");
      }
    }
  }

  outfile->Write();
  outfile->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  if(!EnergyBinning.size()){
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "40,0,10"};
    handleOpts(argc_dum, argv_dum);
  }

  if(!OffAxisSteps.size()){
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-x", "0"};
    handleOpts(argc_dum, argv_dum);
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
  AllInOneGo(*dk2nuRdr, TotalPOT);
}
