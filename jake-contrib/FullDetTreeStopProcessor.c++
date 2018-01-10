#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "FullDetTreeReader.h"

#include "TTree.h"

#include <utility>

struct DetBox {
  double XOffset;
  double XWidth_fv;
  double YWidth_fv;
  double ZWidth_fv;

  double XWidth_det;
  double YWidth_det;
  double ZWidth_det;

  size_t X_fv[2];
  size_t X_veto_left[2];
  size_t X_veto_right[2];
};

struct EDep {
  double vtx[3];

  double LepDep_det;
  double ProtonDep_FV;
  double ProtonDep_veto;
  double NeutronDep_FV;
  double NeutronDep_veto;
  double PiCDep_FV;
  double PiCDep_veto;
  double Pi0Dep_FV;
  double Pi0Dep_veto;
  double OtherDep_FV;
  double OtherDep_veto;

  double Enu;
  int NuPDG;
  int LepPDG;
  double TotalDep_FV;
  double TotalDep_veto;
};

std::vector<DetectorStop> DetectorStops;
std::vector<std::tuple<DetBox, EDep, TTree *> > Detectors;

std::string inpDir, outputFile;
std::string runPlanCfg, runPlanName = "";

double X_Veto_cm = 50, Y_Veto_cm = 50, Z_Veto_cm = 50;
double XShift_cm = 1800;
double XArLen_cm = 4000;
double XBinWidth_cm = 10;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <input dir>                : Input directory containing "
               "condensed arbox.py output.\n"
               "\t-r <RunPlan.XML>              : An XML file specifying a run "
               "plan  \n"
               "\t                                to build fluxes for. See     "
               "      \n"
               "\t                                documentation for XML "
               "structure.   \n"
               "\t-o <output.root>             : Output file name.\n"
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
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }

    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int BinSearch(std::vector<double> &BinEdges, double find) {
  if ((find < BinEdges.front()) || (find > BinEdges.back())) {
    std::cout << "[ERROR]: Failed to find " << find << " in binning ["
              << BinEdges.front() << ", " << BinEdges.back() << "]."
              << std::endl;
    throw;
  }

  int max = BinEdges.size() - 1;
  int min = 1;
  while (true) {
    if (min == max) {
      std::cout << "[ERROR]: Couldn't find " << find
                << " but should have been able to. Last check was " << min
                << std::endl;
      throw;
    }

    int bi = min + floor((max - min) / 2);

    std::cout << "[BS]: min = " << min << ", max = " << max
              << ", check = " << bi << " (" << BinEdges[bi - 1]
              << "), looking for = " << find << std::endl;

    bool eq = (fabs(find - BinEdges[bi - 1]) < 1E-4);
    if (eq) {
      std::cout << "[BS]\tFound equal in " << bi << " @ " << BinEdges[bi - 1]
                << std::endl;
      return bi - 1;
    }

    bool lt = (find < BinEdges[bi - 1]);
    if (lt) {
      std::cout << "[BS] " << find << " < " << BinEdges[bi - 1] << std::endl;
      max = bi;
      continue;
    }

    bool ib = (find < BinEdges[bi]);
    if (ib) {
      std::cout << "[BS]\tFound in bin " << bi << " @ " << BinEdges[bi - 1]
                << std::endl;
      return bi - 1;
    }

    std::cout << "[BS] " << find << " > " << BinEdges[bi] << std::endl;

    min = bi + 1;
    continue;
  }
}

double Jaccumulate(double (&arr)[400][3][3], size_t ind1, size_t ind2,
                   double init) {
  double sum = init;
  for (size_t i = ind1; i < ind2; ++i) {
    sum += arr[i][1][1];
  }
  return sum;
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);
  FullDetTreeReader *rdr = new FullDetTreeReader("fullDetTree", inpDir);

  if (!rdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }
  if (!runPlanCfg.length()) {
    std::cout << "[ERROR]: Found no run plan configuration file." << std::endl;
  }

  TFile *outfile = CheckOpenFile(outputFile, "RECREATE");

  // Read in det stops
  DetectorStops = ReadDetectorStopConfig(runPlanCfg, runPlanName);

  std::vector<double> XBins;
  XBins.push_back(XShift_cm - XArLen_cm / 2.0);
  for (size_t bi_it = 1; bi_it < 40; ++bi_it) {
    XBins.push_back(XBins.back() + XBinWidth_cm);
  }
  XBins.push_back(XBins.back() + XBinWidth_cm);

  // translate to detboxes
  for (auto const &ds : DetectorStops) {
    DetBox db;

    db.XOffset = ds.LateralOffset * 100.0;
    db.XWidth_fv = ds.DetectorFiducialWidth * 100.0;
    db.YWidth_fv = ds.DetectorFiducialHeight * 100.0;
    db.ZWidth_fv = ds.DetectorFiducialDepth * 100.0;
    db.XWidth_det = db.XWidth_fv + 2 * X_Veto_cm;
    db.YWidth_det = db.YWidth_fv + 2 * Y_Veto_cm;
    db.ZWidth_det = db.ZWidth_fv + 2 * Z_Veto_cm;

    double FVLow = db.XOffset - db.XWidth_fv / 2.0;
    double FVHigh = db.XOffset + db.XWidth_fv / 2.0;
    double Detlow = db.XOffset - db.XWidth_det / 2.0;
    double DetHigh = db.XOffset + db.XWidth_det / 2.0;

    db.X_fv[0] = BinSearch(XBins, FVLow);
    db.X_fv[1] = BinSearch(XBins, FVHigh);

    db.X_veto_left[0] = BinSearch(XBins, Detlow);
    db.X_veto_left[1] = BinSearch(XBins, Detlow + X_Veto_cm);

    db.X_veto_right[0] = BinSearch(XBins, DetHigh - X_Veto_cm);
    db.X_veto_right[1] = BinSearch(XBins, DetHigh);

    TTree *stoptree = new TTree(
        (std::string("EDep_Stop") + to_str(db.XOffset) + "_m").c_str(), "");

    Detectors.push_back(std::make_tuple(db, EDep(), stoptree));

    stoptree->Branch("vtx", &std::get<1>(Detectors.back()).vtx, "vtx[3]/D");
    stoptree->Branch("LepDep_det", &std::get<1>(Detectors.back()).LepDep_det,
                     "LepDep_det/D");
    stoptree->Branch("ProtonDep_FV",
                     &std::get<1>(Detectors.back()).ProtonDep_FV,
                     "ProtonDep_FV/D");
    stoptree->Branch("ProtonDep_veto",
                     &std::get<1>(Detectors.back()).ProtonDep_veto,
                     "ProtonDep_veto/D");
    stoptree->Branch("NeutronDep_FV",
                     &std::get<1>(Detectors.back()).NeutronDep_FV,
                     "NeutronDep_FV/D");
    stoptree->Branch("NeutronDep_veto",
                     &std::get<1>(Detectors.back()).NeutronDep_veto,
                     "NeutronDep_veto/D");
    stoptree->Branch("PiCDep_FV", &std::get<1>(Detectors.back()).PiCDep_FV,
                     "PiCDep_FV/D");
    stoptree->Branch("PiCDep_veto", &std::get<1>(Detectors.back()).PiCDep_veto,
                     "PiCDep_veto/D");
    stoptree->Branch("Pi0Dep_FV", &std::get<1>(Detectors.back()).Pi0Dep_FV,
                     "Pi0Dep_FV/D");
    stoptree->Branch("Pi0Dep_veto", &std::get<1>(Detectors.back()).Pi0Dep_veto,
                     "Pi0Dep_veto/D");
    stoptree->Branch("OtherDep_FV", &std::get<1>(Detectors.back()).OtherDep_FV,
                     "OtherDep_FV/D");
    stoptree->Branch("OtherDep_veto",
                     &std::get<1>(Detectors.back()).OtherDep_veto,
                     "OtherDep_veto/D");

    stoptree->Branch("Enu", &std::get<1>(Detectors.back()).Enu, "Enu/D");
    stoptree->Branch("NuPDG", &std::get<1>(Detectors.back()).NuPDG, "NuPDG/I");
    stoptree->Branch("LepPDG", &std::get<1>(Detectors.back()).LepPDG,
                     "LepPDG/I");

    stoptree->Branch("TotalDep_FV", &std::get<1>(Detectors.back()).TotalDep_FV,
                     "TotalDep_FV/D");
    stoptree->Branch("TotalDep_veto",
                     &std::get<1>(Detectors.back()).TotalDep_veto,
                     "TotalDep_veto/D");
  }

  for (size_t e_it = 0; e_it < rdr->GetEntries(); ++e_it) {
    rdr->GetEntry(e_it);
    // DetBox + Reader -> EDep -> Fill out tree

    for (std::tuple<DetBox, EDep, TTree *> &db : Detectors) {
      std::get<1>(db).vtx[0] = rdr->vtx_X;
      std::get<1>(db).vtx[1] = rdr->vtx_Y;
      std::get<1>(db).vtx[2] = rdr->vtx_Z;

      std::get<1>(db).LepDep_det =
          Jaccumulate(rdr->eLepPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->eLepSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).ProtonDep_FV =
          Jaccumulate(rdr->eProtonPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).ProtonDep_veto =
          Jaccumulate(rdr->eProtonPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eProtonPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).NeutronDep_FV =
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->eNeutronSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).NeutronDep_veto =
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eNeutronSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->eNeutronSecondaryDep,
                      std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).PiCDep_FV =
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).PiCDep_veto =
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).Pi0Dep_FV =
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).Pi0Dep_veto =
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).OtherDep_FV =
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).OtherDep_veto =
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0) +
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0);

      std::get<1>(db).Enu = rdr->Enu;
      std::get<1>(db).NuPDG = rdr->nuPDG;
      std::get<1>(db).LepPDG = rdr->lepPDG;
      std::get<1>(db).TotalDep_FV =
          std::get<1>(db).ProtonDep_FV + std::get<1>(db).NeutronDep_FV +
          std::get<1>(db).PiCDep_FV + std::get<1>(db).Pi0Dep_FV +
          std::get<1>(db).OtherDep_FV;

      std::get<1>(db).TotalDep_veto =
          std::get<1>(db).ProtonDep_veto + std::get<1>(db).NeutronDep_veto +
          std::get<1>(db).PiCDep_veto + std::get<1>(db).Pi0Dep_veto +
          std::get<1>(db).OtherDep_veto;

      std::get<2>(db)->Fill();
    }
  }

  for (std::tuple<DetBox, EDep, TTree *> &db : Detectors) {
    std::get<2>(db)->Write(std::get<2>(db)->GetName(), TObject::kOverwrite);
  }

  outfile->Write();
  outfile->Close();
}
