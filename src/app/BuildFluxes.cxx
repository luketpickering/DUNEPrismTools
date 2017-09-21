#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

int NuPDGTarget = 14;

struct NuDataReader {
  Int_t Ntype;
  Double_t Nimpwt;

  Float_t Vx;
  Float_t Vy;
  Float_t Vz;

  Float_t pdPx;
  Float_t pdPy;
  Float_t pdPz;

  Float_t ppdxdz;
  Float_t ppdydz;
  Float_t pppz;
  Float_t ppenergy;

  Float_t muparpx;
  Float_t muparpy;
  Float_t muparpz;
  Float_t mupare;

  Int_t ptype;

  Float_t Necm;

  void SetBranchAddresses(TTree *t) {
    t->SetBranchAddress("Ntype", &Ntype);
    t->SetBranchAddress("Nimpwt", &Nimpwt);
    t->SetBranchAddress("Vx", &Vx);
    t->SetBranchAddress("Vy", &Vy);
    t->SetBranchAddress("Vz", &Vz);
    t->SetBranchAddress("pdPx", &pdPx);
    t->SetBranchAddress("pdPy", &pdPy);
    t->SetBranchAddress("pdPz", &pdPz);
    t->SetBranchAddress("ppdxdz", &ppdxdz);
    t->SetBranchAddress("ppdydz", &ppdydz);
    t->SetBranchAddress("pppz", &pppz);
    t->SetBranchAddress("ppenergy", &ppenergy);
    t->SetBranchAddress("muparpx", &muparpx);
    t->SetBranchAddress("muparpy", &muparpy);
    t->SetBranchAddress("muparpz", &muparpz);
    t->SetBranchAddress("mupare", &mupare);
    t->SetBranchAddress("ptype", &ptype);
    t->SetBranchAddress("Necm", &Necm);
  }

  Int_t GetPDG() {
    switch (Ntype) {
      case 52:
        return -12;
      case 53:
        return 12;

      case 55:
        return -14;
      case 56:
        return 14;

      default:
        return 0;
    }
  }

  double GetParentMass() {
    static double const pimass = 0.13957;  // in GeV
    static double const kmass = 0.49368;
    static double const k0mass = 0.49767;
    static double const mumass = 0.105658389;

    if ((ptype == 8) || (ptype == 9)) {
      return pimass;
    } else if ((ptype == 11) || (ptype == 12)) {
      return kmass;
    } else if (ptype == 10) {
      return k0mass;
    } else if ((ptype == 5) || (ptype == 6)) {
      return mumass;
    } else {
      return 0xdeadbeef;
    }
  }
};

double CosTheta(TVector3 const &v1, TVector3 const &v2) {
  return v1.Unit().Dot(v2.Unit());
}

std::pair<double, double> GetNuWeight(NuDataReader &ndr,
                                      TVector3 const &DetPoint) {
  static const double detRadius = 100.0;  // in cm
  static double const mumass = 0.105658389;

  double parent_mass = ndr.GetParentMass();
  if (parent_mass == 0xdeadbeef) {
    return std::make_pair(0xdeadbeef, 0xdeadbeef);
  }

  TLorentzVector parent_4mom;
  parent_4mom.SetXYZM(ndr.pdPx, ndr.pdPy, ndr.pdPz, parent_mass);

  double parent_energy = parent_4mom.E();

  double enuzr = ndr.Necm;

  TVector3 VToDetVect((DetPoint[0] - ndr.Vx), (DetPoint[1] - ndr.Vy),
                      (DetPoint[2] - ndr.Vz));

  double parent_4mom_Gamma = parent_4mom.Gamma();
  double emrat =
      1.0 /
      (parent_4mom_Gamma *
       (1 - parent_4mom.Beta() * CosTheta(parent_4mom.Vect(), VToDetVect)));

  double nu_energy = emrat * enuzr;

  double sangDet = (detRadius * detRadius / VToDetVect.Mag2() / 4.);

  double nu_wght = sangDet * emrat * emrat;

  // done for all except polarized muon
  // in which case need to modify weight
  if ((ndr.ptype == 5) || (ndr.ptype == 6)) {
    // boost new neutrino to mu decay cm
    TVector3 betaVect = parent_4mom.Vect();
    betaVect[0] /= parent_4mom.E();
    betaVect[1] /= parent_4mom.E();
    betaVect[2] /= parent_4mom.E();

    TVector3 nu_3mom = VToDetVect.Unit() * nu_energy;

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
    parent_4mom_Gamma = ndr.ppenergy / parent_mass;
    TVector3 parentParentBetaVect(ndr.ppdxdz * ndr.pppz, ndr.ppdydz * ndr.pppz,
                                  ndr.pppz);
    parentParentBetaVect[0] /= ndr.ppenergy;
    parentParentBetaVect[1] /= ndr.ppenergy;
    parentParentBetaVect[2] /= ndr.ppenergy;

    double partial_mp = parent_4mom_Gamma *
                        (betaVect[0] * ndr.muparpx + betaVect[1] * ndr.muparpy +
                         betaVect[2] * ndr.muparpz);
    partial_mp = ndr.mupare - partial_mp / (parent_4mom_Gamma + 1.);

    TLorentzVector p_pcm_mp;
    p_pcm_mp[0] = ndr.muparpx - betaVect[0] * parent_4mom_Gamma * partial_mp;
    p_pcm_mp[1] = ndr.muparpy - betaVect[1] * parent_4mom_Gamma * partial_mp;
    p_pcm_mp[2] = ndr.muparpz - betaVect[2] * parent_4mom_Gamma * partial_mp;
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
      if (abs(ndr.GetPDG()) == NuPDGTarget) {
        wt_ratio = 1. - costh;
      } else if (abs(ndr.GetPDG()) == NuPDGTarget) {
        double xnu = 2. * enuzr / mumass;
        wt_ratio =
            ((3. - 2. * xnu) - (1. - 2. * xnu) * costh) / (3. - 2. * xnu);
      } else {
        std::cout << "eventRates:: Bad neutrino type = " << ndr.Ntype
                  << std::endl;
        throw;
      }
    }
    nu_wght *= wt_ratio;
  }

  return std::make_pair(nu_energy, nu_wght);
}

double str2d(std::string str) {
  std::istringstream stream(str);
  double d;
  stream >> d;
  return d;
}

int str2i(std::string str) {
  std::istringstream stream(str);
  int d;
  stream >> d;
  return d;
}

std::vector<double> ParseToDbl(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<double> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    std::istringstream stream(temp_string);
    double entry;
    stream >> entry;

    vals.push_back(entry);
  }
  return vals;
}

std::vector<std::string> ParseToStr(std::string str, const char *del) {
  std::istringstream stream(str);
  std::string temp_string;
  std::vector<std::string> vals;

  while (std::getline(stream >> std::ws, temp_string, *del)) {
    if (temp_string.empty()) continue;
    vals.push_back(temp_string);
  }

  return vals;
}

std::vector<double> mrads;
std::string inpDir = ".";
bool DoExtra = false;

int NBins = 50;
double BLow = 0;
double BUp = 10;

std::vector<double> varBin;
bool VariableBinning = false;

bool Rebin = false;
int rNBins = 50;
double rBLow = 0;
double rBUp = 10;

bool DontSaveCalcBin = false;
bool DontSmoothInterpolation = false;

std::string outputFile;

bool DoBinOptimiz = false;

int NBO = 5E3;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << " -i /path/to/DUNE/nudata/files -a "
               "1,2,3,4,5,<other mrad fluxes to calc> -o [<Output ROOT file>] "
               "[-e (Outputs flux from "
               "specific parent species.)] [-b <NBins>,<FluxLow>,<FluxHigh>] "
               "[-d <step_deg>,<max_deg>] [-vb "
               "<bin0low>,<bin1low>_<binXUp>:step,..,<binYlow>,<binYup>] -r "
               "<RebinNBins>,<RebinFluxLow>,<RebinFluxHigh> -BO <NNuPerBin> "
               "(BinOptimisation -- Warning, very slow) -C (keep calc fluxes) "
               "-I (Don't smooth interpolation) -u <NuPDG>"
            << std::endl;
}

void AppendDblVect(std::vector<double> &target,
                   std::vector<double> const &toApp) {
  for (size_t i = 0; i < toApp.size(); ++i) {
    target.push_back(toApp[i]);
  }
}

// Converts "5_10:1" into a vector containing: 5,6,7,8,9,10
std::vector<double> BuildDoubleList(std::string const &str) {
  std::vector<std::string> steps = ParseToStr(str, ":");
  if (steps.size() != 2) {
    return ParseToDbl(str, ",");
  }
  double step = str2d(steps[1]);

  std::vector<double> range = ParseToDbl(steps[0], "_");
  if (!steps.size() == 2) {
    std::cout
        << "[ERROR]: When attempting to parse bin range descriptor: \" " << str
        << "\", couldn't determine range. Expect form: <bin1low>_<binXUp>:step"
        << std::endl;
    exit(1);
  }

  int nsteps = (range[1] - range[0]) / step;

  std::vector<double> rtn;
  for (int step_it = 0; step_it <= nsteps; ++step_it) {
    rtn.push_back(range[0] + step * step_it);
  }
  return rtn;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-i") {
      inpDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-a") {
      mrads = ParseToDbl(argv[++opt], ",");
      mrads.insert(mrads.begin(), 0);
    } else if (std::string(argv[opt]) == "-C") {
      DontSaveCalcBin = true;
    } else if (std::string(argv[opt]) == "-I") {
      DontSmoothInterpolation = true;
    } else if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-e") {
      DoExtra = true;
    } else if (std::string(argv[opt]) == "-b") {
      std::vector<double> binning = ParseToDbl(argv[++opt], ",");
      if (binning.size() != 3) {
        std::cout << "[ERROR]: Recieved " << binning.size()
                  << " entrys for -b, expected 3." << std::endl;
        exit(1);
      }
      NBins = int(binning[0]);
      BLow = binning[1];
      BUp = binning[2];
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<double> binning = ParseToDbl(argv[++opt], ",");
      if (binning.size() != 3) {
        std::cout << "[ERROR]: Recieved " << binning.size()
                  << " entrys for -r, expected 3." << std::endl;
        exit(1);
      }
      rNBins = int(binning[0]);
      rBLow = binning[1];
      rBUp = binning[2];
      Rebin = true;
    } else if (std::string(argv[opt]) == "-vb") {
      std::vector<std::string> vbDescriptors = ParseToStr(argv[++opt], ",");
      varBin.clear();
      for (size_t vbd_it = 0; vbd_it < vbDescriptors.size(); ++vbd_it) {
        AppendDblVect(varBin, BuildDoubleList(vbDescriptors[vbd_it]));
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
      std::vector<double> radStep = ParseToDbl(argv[++opt], ",");
      if (radStep.size() != 2) {
        std::cout << "[ERROR]: Recieved " << radStep.size()
                  << " entrys for -d, expected 2." << std::endl;
        exit(1);
      }
      mrads.push_back(0);
      static const double deg2rad = asin(1) / 90.0;
      double curr = radStep[0];
      while (curr < radStep[1]) {
        mrads.push_back(curr * deg2rad * 1E3);
        curr += radStep[0];
      }

    } else if (std::string(argv[opt]) == "-BO") {
      NBO = str2i(argv[++opt]);
      DoBinOptimiz = true;
    } else if (std::string(argv[opt]) == "-u") {
      NuPDGTarget = str2i(argv[++opt]);
      if ((abs(NuPDGTarget) != 14) && (abs(NuPDGTarget) != 12)) {
        std::cout << "[ERROR]: -u option recieved " << argv[opt - 1]
                  << ". Expected one of {12,-12,14,-14}." << std::endl;
      }
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

void AllInOneGo(TTree *nuChain, double TotalPOT, double ZDist) {
  NuDataReader ndr;
  ndr.SetBranchAddresses(nuChain);

  size_t NAngs = mrads.size();
  size_t nNus = nuChain->GetEntries();
  std::cout << "Have " << nNus << " nudata entries." << std::endl;

  TFile *outfile = new TFile(outputFile.c_str(), "RECREATE");

  std::vector<std::vector<TH1D *> > Hists;
  Hists.resize(NAngs);
  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    std::stringstream ss("");
    ss << "LBNF_numu_mrad_" << mrads[ang_it];

    Hists[ang_it].push_back(
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
      Hists[ang_it].push_back(
          new TH1D((ss.str() + "_pdp").c_str(), "", 1010, -10, 1000));

      Hists[ang_it].push_back(
          VariableBinning
              ? new TH1D(
                    (ss.str() + "_pi").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    varBin.size() - 1, varBin.data())
              : new TH1D(
                    (ss.str() + "_pi").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    NBins, BLow, BUp));
      Hists[ang_it].push_back(
          VariableBinning
              ? new TH1D(
                    (ss.str() + "_k").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    varBin.size() - 1, varBin.data())
              : new TH1D(
                    (ss.str() + "_k").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    NBins, BLow, BUp));
      Hists[ang_it].push_back(
          VariableBinning
              ? new TH1D(
                    (ss.str() + "_k0").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    varBin.size() - 1, varBin.data())
              : new TH1D(
                    (ss.str() + "_k0").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    NBins, BLow, BUp));
      Hists[ang_it].push_back(
          VariableBinning
              ? new TH1D(
                    (ss.str() + "_mu").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    varBin.size() - 1, varBin.data())
              : new TH1D(
                    (ss.str() + "_mu").c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    NBins, BLow, BUp));
    }
  }

  std::vector<TVector3> detPoses;
  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    detPoses.push_back(TVector3(ZDist * mrads[ang_it] * 1E-3, 0, ZDist));

    static const double rad2deg = 90.0 / asin(1);

    std::cout << "[INFO]: Building flux at: "
              << (mrads[ang_it] * 1.0E-3 * rad2deg) << " degrees ("
              << mrads[ang_it] << " mrads) Offset = "
              << ((ZDist * 1.0E-2) * tan(mrads[ang_it] * 1.0E-3)) << " m."
              << std::endl;
  }
  for (size_t nu_it = 0; nu_it < nNus; ++nu_it) {
    nuChain->GetEntry(nu_it);

    if (ndr.GetPDG() != NuPDGTarget) {  // Only care about numu at the moment
      continue;
    }

    double wF = (ndr.Nimpwt / TMath::Pi()) * (1.0 / TotalPOT);
    for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
      std::pair<double, double> nuStats = GetNuWeight(ndr, detPoses[ang_it]);

      double w = nuStats.second * wF;
      Hists[ang_it][0]->Fill(nuStats.first, w);

      if ((nuStats.first != nuStats.first) || (w != w)) {
        std::cout << nuStats.first << ", " << w << std::endl;
        throw;
      }

      if (DoExtra) {
        Hists[ang_it][1]->Fill(ndr.Vz * 1E-2, w);

        if ((ndr.ptype == 8) || (ndr.ptype == 9)) {
          Hists[ang_it][2]->Fill(nuStats.first, w);
        } else if ((ndr.ptype == 11) || (ndr.ptype == 12)) {
          Hists[ang_it][3]->Fill(nuStats.first, w);
        } else if (ndr.ptype == 10) {
          Hists[ang_it][4]->Fill(nuStats.first, w);
        } else if ((ndr.ptype == 5) || (ndr.ptype == 6)) {
          Hists[ang_it][5]->Fill(nuStats.first, w);
        }
      }
    }
  }
  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    Hists[ang_it][0]->Scale(1E-4, "width");
    if (DoExtra) {
      Hists[ang_it][2]->Scale(1E-4, "width");
      Hists[ang_it][3]->Scale(1E-4, "width");
      Hists[ang_it][4]->Scale(1E-4, "width");
      Hists[ang_it][5]->Scale(1E-4, "width");
    }
  }

  if (Rebin) {
    std::vector<std::vector<TH1D *> > ReBinHists;

    for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
      ReBinHists.push_back(std::vector<TH1D *>());
      for (size_t h_it = 0; h_it < Hists[ang_it].size(); ++h_it) {
        if (h_it == 1) {  // skip PDP
          continue;
        }
        std::string name = Hists[ang_it][h_it]->GetName();
        Hists[ang_it][h_it]->SetName((name + "_CalcBinning").c_str());
        if (DontSaveCalcBin) {
          Hists[ang_it][h_it]->SetDirectory(NULL);
        }
        TGraph TheEvaluator(Hists[ang_it][h_it]->GetXaxis()->GetNbins());

        for (size_t bi_it = 1;
             bi_it < Hists[ang_it][h_it]->GetXaxis()->GetNbins() + 1; ++bi_it) {
          TheEvaluator.SetPoint(
              bi_it - 1, Hists[ang_it][h_it]->GetXaxis()->GetBinCenter(bi_it),
              Hists[ang_it][h_it]->GetBinContent(bi_it));
        }

        ReBinHists[ang_it].push_back(
            new TH1D(name.c_str(),
                     ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                     rNBins, rBLow, rBUp));

        for (size_t bin_it = 1;
             bin_it < ReBinHists[ang_it].back()->GetXaxis()->GetNbins() + 1;
             ++bin_it) {
          double eval = TheEvaluator.Eval(
              ReBinHists[ang_it].back()->GetXaxis()->GetBinCenter(bin_it), 0,
              "S");
          ReBinHists[ang_it].back()->SetBinContent(bin_it, eval);
          ReBinHists[ang_it].back()->SetBinError(bin_it, eval * 0.01);
        }
      }
    }
  }

  outfile->Write();
  outfile->Close();
}

// Only want to sort on pair.first
bool SortNumuByE(std::pair<double, double> const &l,
                 std::pair<double, double> const &r) {
  return l.first < r.first;
}

void BinOptimisation(TTree *nuChain, double TotalPOT, double ZDist) {
  NuDataReader ndr;
  ndr.SetBranchAddresses(nuChain);

  size_t NAngs = mrads.size();
  size_t nNus = nuChain->GetEntries();
  std::cout << "Have " << nNus << " nudata entries." << std::endl;

  TFile *outfile = new TFile(outputFile.c_str(), "RECREATE");

  for (size_t ang_it = 0; ang_it < NAngs; ++ang_it) {
    TVector3 detPos(ZDist * mrads[ang_it] * 1E-3, 0, ZDist);

    static const double rad2deg = 90.0 / asin(1);

    std::cout << "[INFO]: Building flux at: "
              << (mrads[ang_it] * 1.0E-3 * rad2deg) << " degrees ("
              << mrads[ang_it] << " mrads) Offset = "
              << ((ZDist * 1.0E-2) * tan(mrads[ang_it] * 1.0E-3)) << " m."
              << std::endl;

    TH1D *Hist;
    std::vector<std::pair<double, double> > numus;
    for (size_t nu_it = 0; nu_it < nNus; ++nu_it) {
      nuChain->GetEntry(nu_it);

      if (ndr.GetPDG() != NuPDGTarget) {  // Only care about numu at the moment
        continue;
      }

      double wF = (ndr.Nimpwt / TMath::Pi()) * (1.0 / TotalPOT);
      numus.push_back(GetNuWeight(ndr, detPos));
      numus.back().second *= wF;
    }
    std::cout << "Sorting..." << std::endl;
    std::stable_sort(numus.begin(), numus.end(), SortNumuByE);
    std::cout << "Binning..." << std::endl;

    std::vector<double> binEdges;

    size_t InBin = 0;
    size_t NNumus = numus.size();

    if (!NNumus) {
      std::cout << "[ERROR]: Didn't find any neutrinos." << std::endl;
      exit(1);
    }
    binEdges.push_back(0);
    for (size_t nu_it = 0; nu_it < (NNumus - 1); ++nu_it) {
      if (InBin < NBO) {
        InBin++;
        continue;
      }
      binEdges.push_back((numus[nu_it + 1].first + numus[nu_it].first) / 2.0);
      // std::cout << "Bin [" << binEdges.size()
      //           << "] low edge: " << binEdges.back() << std::endl;
      InBin = 0;
    }
    binEdges.push_back(numus[NNumus - 1].first +
                       (numus[NNumus - 1].first - numus[NNumus - 2].first) /
                           2.0);

    std::stringstream ss("");
    ss << "LBNF_numu_mrad_" << mrads[ang_it];
    Hist = new TH1D(ss.str().c_str(),
                    ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                    binEdges.size() - 1, binEdges.data());

    for (size_t nu_it = 0; nu_it < (NNumus - 1); ++nu_it) {
      Hist->Fill(numus[nu_it].first, numus[nu_it].second);
    }

    Hist->Scale(1E-4, "width");

    if (Rebin) {
      std::string name = Hist->GetName();
      Hist->SetName((name + "_CalcBinning").c_str());
      if (DontSaveCalcBin) {
        Hist->SetDirectory(NULL);
      }
      TH1D *ReBinHist =
          new TH1D(name.c_str(),
                   ";#it{E}_{#nu} (GeV);#Phi_{#nu} (GeV^{-1}cm^{-2} per POT)",
                   rNBins, rBLow, rBUp);

      TGraph TheEvaluator(Hist->GetXaxis()->GetNbins());

      for (size_t bi_it = 1; bi_it < Hist->GetXaxis()->GetNbins() + 1;
           ++bi_it) {
        TheEvaluator.SetPoint(bi_it - 1, Hist->GetXaxis()->GetBinCenter(bi_it),
                              Hist->GetBinContent(bi_it));
      }

      double Peak = std::numeric_limits<double>::min();
      bool foundPeak = false;
      double peval = 0;
      for (size_t bin_it = 1; bin_it < ReBinHist->GetXaxis()->GetNbins() + 1;
           ++bin_it) {
        double eval;
        if (DontSmoothInterpolation) {
          eval = Hist->Interpolate(ReBinHist->GetXaxis()->GetBinCenter(bin_it));
        } else {
          eval = TheEvaluator.Eval(ReBinHist->GetXaxis()->GetBinCenter(bin_it),
                                   0, "S");
          eval = eval < 0 ? 0 : eval;

          if (!foundPeak) {
            if (eval > Peak) {
              Peak = eval;
            } else {
              foundPeak = true;
            }
          } else {
            if (eval > peval) {
              eval = TheEvaluator.Eval(
                  ReBinHist->GetXaxis()->GetBinCenter(bin_it));
              if (eval > peval) {
                eval = peval;
              }
            }
          }
        }

        if (ReBinHist->GetXaxis()->GetBinCenter(bin_it) >
            Hist->GetXaxis()->GetBinCenter(Hist->GetXaxis()->GetNbins())) {
          eval = 0;
        }

        ReBinHist->SetBinContent(bin_it, eval);
        ReBinHist->SetBinError(bin_it,
                               eval * (sqrt(double(NBO)) / double(NBO)));
        peval = eval;
      }
    }
  }
  outfile->Write();
  outfile->Close();
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  double ZDist_ND = 57400;  // cm
  double ZDist_FD = 1287E5;  // cm

  double ZDist = ZDist_FD;

  double PotPerFile = 1E5;

  TChain *nuChain = new TChain("nudata");
  std::string loc = inpDir + "/g4lbne_v3r3p8_QGSP_BERT_DK2Nu_200kA*.root";
  int nFiles = nuChain->Add(loc.c_str());

  double TotalPOT = PotPerFile * nFiles;
  nuChain->SetMakeClass(true);

  DoBinOptimiz ? BinOptimisation(nuChain, TotalPOT, ZDist)
               : AllInOneGo(nuChain, TotalPOT, ZDist);
}
