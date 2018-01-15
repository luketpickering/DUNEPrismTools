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

  double X_Range_fv[2];
  double Y_Range_fv[2];
  double Z_Range_fv[2];

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

  bool flagLepExitBack;
  bool flagLepExitFront;
  bool flagLepExitY;
  bool flagLepExitXLow;
  bool flagLepExitXHigh;

  double lepExitingPosX;
  double lepExitingPosY;
  double lepExitingPosZ;
  double lepExitingMomX;
  double lepExitingMomY;
  double lepExitingMomZ;

  double Q2True;
  double yTrue;
  double W_rest;

};

std::vector<DetectorStop> DetectorStops;
std::vector<std::tuple<DetBox, EDep, TTree *> > Detectors;

std::string inpDir, outputFile;
std::string runPlanCfg, runPlanName = "";

double X_Veto_cm = 50, Y_Veto_cm = 50, Z_Veto_cm = 50;
double XShift_cm = -1800;
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
  // Special case -- if your up-bin is at the top of the range then
  // give them the 'overflow' bin as this will be set as the uncount top
  // of the sum range.
  if ((fabs(find - BinEdges.back()) < 1E-4)) {
    return BinEdges.size() - 1;
  }

  if ((find < BinEdges.front()) || (find > BinEdges.back())) {
    std::cout << "[ERROR]: Failed to find " << find << " in binning ["
              << BinEdges.front() << ", " << BinEdges.back() << "]."
              << std::endl;
    throw;
  }

  int max = BinEdges.size() - 1;
  int min = 1;

  // std::cout << "[BS]: Looking for " << find << " in binning ["
  //           << BinEdges.front() << ", " << BinEdges.back() << "]." <<
  //           std::endl;

  while (true) {
    if (min == max) {
      std::cout << "[ERROR]: Couldn't find " << find
                << " but should have been able to. Last check was " << min
                << std::endl;
      throw;
    }

    int bi = min + floor((max - min) / 2);

    // std::cout << "[BS]: min = " << min << ", max = " << max
    //           << ", check = " << bi << " (" << BinEdges[bi - 1]
    //           << "), looking for = " << find << std::endl;

    bool eq = (fabs(find - BinEdges[bi - 1]) < 1E-4);
    if (eq) {
      // std::cout << "[BS]\tFound equal in " << bi << " @ " << BinEdges[bi - 1]
      //           << std::endl;
      return bi - 1;
    }

    bool lt = (find < BinEdges[bi - 1]);
    if (lt) {
      // std::cout << "[BS] " << find << " < " << BinEdges[bi - 1] << std::endl;
      max = bi;
      continue;
    }

    bool ib = (find < BinEdges[bi]);
    if (ib) {
      // std::cout << "[BS]\tFound in bin " << bi << " @ " << BinEdges[bi - 1]
      //           << std::endl;
      return bi - 1;
    }

    // std::cout << "[BS] " << find << " > " << BinEdges[bi] << std::endl;

    min = bi + 1;
    continue;
  }
}

double Jaccumulate(double (&arr)[400][3][3], size_t ind1, size_t ind2,
                   double init, bool IncludeFVYZ = true,
                   bool IncludeOOFVYZ = false) {
  double sum = init;
  for (size_t i = ind1; i < ind2; ++i) {
    if ((arr[i][1][1] != 0) && (!std::isnormal(arr[i][1][1]))) {
      std::cout << "[" << i << "][1][1]: bin non-normal = " << arr[i][1][1]
                << std::endl;
    }

    if (IncludeFVYZ) {
      sum += arr[i][1][1];
    }

    if (IncludeOOFVYZ) {
      sum += arr[i][0][0];
      sum += arr[i][0][1];
      sum += arr[i][0][2];

      sum += arr[i][1][0];
      sum += arr[i][1][2];

      sum += arr[i][2][0];
      sum += arr[i][2][1];
      sum += arr[i][2][2];

      if ((arr[i][0][0] != 0) && (!std::isnormal(arr[i][0][0]))) {
        std::cout << "[" << i << "][0][0]: bin non-normal = " << arr[i][0][0]
                  << std::endl;
      }
      if ((arr[i][0][1] != 0) && (!std::isnormal(arr[i][0][1]))) {
        std::cout << "[" << i << "][0][1]: bin non-normal = " << arr[i][0][1]
                  << std::endl;
      }
      if ((arr[i][0][2] != 0) && (!std::isnormal(arr[i][0][2]))) {
        std::cout << "[" << i << "][0][2]: bin non-normal = " << arr[i][0][2]
                  << std::endl;
      }

      if ((arr[i][1][0] != 0) && (!std::isnormal(arr[i][1][0]))) {
        std::cout << "[" << i << "][1][0]: bin non-normal = " << arr[i][1][0]
                  << std::endl;
      }
      if ((arr[i][1][2] != 0) && (!std::isnormal(arr[i][1][2]))) {
        std::cout << "[" << i << "][1][2]: bin non-normal = " << arr[i][1][2]
                  << std::endl;
      }

      if ((arr[i][2][0] != 0) && (!std::isnormal(arr[i][2][0]))) {
        std::cout << "[" << i << "][2][0]: bin non-normal = " << arr[i][2][0]
                  << std::endl;
      }
      if ((arr[i][2][1] != 0) && (!std::isnormal(arr[i][2][1]))) {
        std::cout << "[" << i << "][2][1]: bin non-normal = " << arr[i][2][1]
                  << std::endl;
      }
      if ((arr[i][2][2] != 0) && (!std::isnormal(arr[i][2][2]))) {
        std::cout << "[" << i << "][2][2]: bin non-normal = " << arr[i][2][2]
                  << std::endl;
      }
    }
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
  //DetectorStops.resize(1);

  std::vector<double> XBins;

  XBins.push_back(XShift_cm - XArLen_cm / 2.0);
  for (size_t bi_it = 1; bi_it < 400; ++bi_it) {
    XBins.push_back(XBins.back() + XBinWidth_cm);
  }
  XBins.push_back(XBins.back() + XBinWidth_cm);

  size_t NDets = DetectorStops.size();
  Detectors.reserve(NDets);

  // translate to detboxes
  for (size_t d_it = 0; d_it < NDets; ++d_it) {
    DetBox db;

    TTree *stoptree =
        new TTree((std::string("EDep_Stop") +
                   to_str(DetectorStops[d_it].LateralOffset) + "_m")
                      .c_str(),
                  "");

    Detectors.push_back(std::make_tuple(db, EDep(), stoptree));

    std::get<0>(Detectors[d_it]).XOffset =
        -DetectorStops[d_it].LateralOffset * 100.0;
    std::get<0>(Detectors[d_it]).XWidth_fv =
        DetectorStops[d_it].DetectorFiducialWidth * 100.0;
    std::get<0>(Detectors[d_it]).YWidth_fv =
        DetectorStops[d_it].DetectorFiducialHeight * 100.0;
    std::get<0>(Detectors[d_it]).ZWidth_fv =
        DetectorStops[d_it].DetectorFiducialDepth * 100.0;
    std::get<0>(Detectors[d_it]).XWidth_det =
        std::get<0>(Detectors[d_it]).XWidth_fv + 2 * X_Veto_cm;
    std::get<0>(Detectors[d_it]).YWidth_det =
        std::get<0>(Detectors[d_it]).YWidth_fv + 2 * Y_Veto_cm;
    std::get<0>(Detectors[d_it]).ZWidth_det =
        std::get<0>(Detectors[d_it]).ZWidth_fv + 2 * Z_Veto_cm;

    double Detlow = std::get<0>(Detectors[d_it]).XOffset -
                    std::get<0>(Detectors[d_it]).XWidth_det / 2.0;
    double DetHigh = std::get<0>(Detectors[d_it]).XOffset +
                     std::get<0>(Detectors[d_it]).XWidth_det / 2.0;

    std::get<0>(Detectors[d_it]).X_Range_fv[0] =
        std::get<0>(Detectors[d_it]).XOffset -
        std::get<0>(Detectors[d_it]).XWidth_fv / 2.0;
    std::get<0>(Detectors[d_it]).X_Range_fv[1] =
        std::get<0>(Detectors[d_it]).XOffset +
        std::get<0>(Detectors[d_it]).XWidth_fv / 2.0;

    std::get<0>(Detectors[d_it]).Y_Range_fv[0] =
        -std::get<0>(Detectors[d_it]).YWidth_fv / 2.0;
    std::get<0>(Detectors[d_it]).Y_Range_fv[1] =
        std::get<0>(Detectors[d_it]).YWidth_fv / 2.0;

    std::get<0>(Detectors[d_it]).Z_Range_fv[0] =
        -std::get<0>(Detectors[d_it]).ZWidth_fv / 2.0;
    std::get<0>(Detectors[d_it]).Z_Range_fv[1] =
        std::get<0>(Detectors[d_it]).ZWidth_fv / 2.0;

    std::cout << "[INFO]: Det: [" << Detlow << ", " << DetHigh << "], FV: ["
              << std::get<0>(Detectors[d_it]).X_Range_fv[0] << ", "
              << std::get<0>(Detectors[d_it]).X_Range_fv[1] << "]."
              << std::endl;

    std::get<0>(Detectors[d_it]).X_fv[0] =
        BinSearch(XBins, std::get<0>(Detectors[d_it]).X_Range_fv[0]);
    std::get<0>(Detectors[d_it]).X_fv[1] =
        BinSearch(XBins, std::get<0>(Detectors[d_it]).X_Range_fv[1]);

    std::get<0>(Detectors[d_it]).X_veto_left[0] = BinSearch(XBins, Detlow);
    std::get<0>(Detectors[d_it]).X_veto_left[1] =
        BinSearch(XBins, Detlow + X_Veto_cm);

    std::get<0>(Detectors[d_it]).X_veto_right[0] =
        BinSearch(XBins, DetHigh - X_Veto_cm);
    std::get<0>(Detectors[d_it]).X_veto_right[1] = BinSearch(XBins, DetHigh);

    stoptree->Branch("vtx", &std::get<1>(Detectors[d_it]).vtx, "vtx[3]/D");
    stoptree->Branch("LepDep_det", &std::get<1>(Detectors[d_it]).LepDep_det,
                     "LepDep_det/D");
    stoptree->Branch("ProtonDep_FV", &std::get<1>(Detectors[d_it]).ProtonDep_FV,
                     "ProtonDep_FV/D");
    stoptree->Branch("ProtonDep_veto",
                     &std::get<1>(Detectors[d_it]).ProtonDep_veto,
                     "ProtonDep_veto/D");
    stoptree->Branch("NeutronDep_FV",
                     &std::get<1>(Detectors[d_it]).NeutronDep_FV,
                     "NeutronDep_FV/D");
    stoptree->Branch("NeutronDep_veto",
                     &std::get<1>(Detectors[d_it]).NeutronDep_veto,
                     "NeutronDep_veto/D");
    stoptree->Branch("PiCDep_FV", &std::get<1>(Detectors[d_it]).PiCDep_FV,
                     "PiCDep_FV/D");
    stoptree->Branch("PiCDep_veto", &std::get<1>(Detectors[d_it]).PiCDep_veto,
                     "PiCDep_veto/D");
    stoptree->Branch("Pi0Dep_FV", &std::get<1>(Detectors[d_it]).Pi0Dep_FV,
                     "Pi0Dep_FV/D");
    stoptree->Branch("Pi0Dep_veto", &std::get<1>(Detectors[d_it]).Pi0Dep_veto,
                     "Pi0Dep_veto/D");
    stoptree->Branch("OtherDep_FV", &std::get<1>(Detectors[d_it]).OtherDep_FV,
                     "OtherDep_FV/D");
    stoptree->Branch("OtherDep_veto",
                     &std::get<1>(Detectors[d_it]).OtherDep_veto,
                     "OtherDep_veto/D");

    stoptree->Branch("Enu", &std::get<1>(Detectors[d_it]).Enu, "Enu/D");

    // std::get<1>(Detectors[d_it]).Enu = 1;
    // std::cout << "Enu: " << &std::get<1>(Detectors[d_it]).Enu << std::endl;

    stoptree->Branch("NuPDG", &std::get<1>(Detectors[d_it]).NuPDG, "NuPDG/I");
    stoptree->Branch("LepPDG", &std::get<1>(Detectors[d_it]).LepPDG,
                     "LepPDG/I");
    stoptree->Branch("yTrue", &std::get<1>(Detectors[d_it]).yTrue, "yTrue/D");
    stoptree->Branch("W_rest", &std::get<1>(Detectors[d_it]).W_rest, "W_rest/D");
    stoptree->Branch("Q2True", &std::get<1>(Detectors[d_it]).Q2True, "Q2True/D");
    stoptree->Branch("flagLepExitBack", &std::get<1>(Detectors[d_it]).flagLepExitBack,"flagLepExitBack/B");
    stoptree->Branch("flagLepExitFront", &std::get<1>(Detectors[d_it]).flagLepExitFront,"flagLepExitFront/B");
    stoptree->Branch("flagLepExitY", &std::get<1>(Detectors[d_it]).flagLepExitY,"flagLepExitY/B");
    stoptree->Branch("flagLepExitXLow", &std::get<1>(Detectors[d_it]).flagLepExitXLow,"flagLepExitXLow/B");
    stoptree->Branch("flagLepExitXHigh", &std::get<1>(Detectors[d_it]).flagLepExitXHigh,"flagLepExitXHigh/B");
    stoptree->Branch("lepExitingPosX",&std::get<1>(Detectors[d_it]).lepExitingPosX,"lepExitingPosX/D");
    stoptree->Branch("lepExitingPosY",&std::get<1>(Detectors[d_it]).lepExitingPosY,"lepExitingPosY/D");
    stoptree->Branch("lepExitingPosZ",&std::get<1>(Detectors[d_it]).lepExitingPosZ,"lepExitingPosZ/D");
    stoptree->Branch("lepExitingMomX",&std::get<1>(Detectors[d_it]).lepExitingMomX,"lepExitingMomX/D");
    stoptree->Branch("lepExitingMomY",&std::get<1>(Detectors[d_it]).lepExitingMomY,"lepExitingMomY/D");
    stoptree->Branch("lepExitingMomZ",&std::get<1>(Detectors[d_it]).lepExitingMomZ,"lepExitingMomZ/D");

    stoptree->Branch("TotalDep_FV", &std::get<1>(Detectors[d_it]).TotalDep_FV,
                     "TotalDep_FV/D");
    stoptree->Branch("TotalDep_veto",
                     &std::get<1>(Detectors[d_it]).TotalDep_veto,
                     "TotalDep_veto/D");
  }

  std::cout << "[INFO]: Reading " << rdr->GetEntries() << " input entries."
            << std::endl;

  size_t loud_every = rdr->GetEntries() / 10;

  for (size_t e_it = 0; e_it < rdr->GetEntries(); ++e_it) {
    rdr->GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << rdr->vtx_X << ", " << rdr->vtx_Y << ", " << rdr->vtx_Z
                << "}, Enu: " << rdr->Enu << " )" << std::flush;
    }

    // DetBox + Reader -> EDep -> Fill out tree

    for (size_t d_it = 0; d_it < NDets; ++d_it) {
      std::tuple<DetBox, EDep, TTree *> &db = Detectors[d_it];
      std::get<1>(db).vtx[0] = rdr->vtx_X;
      std::get<1>(db).vtx[1] = rdr->vtx_Y;
      std::get<1>(db).vtx[2] = rdr->vtx_Z;

      if ((rdr->vtx_X < std::get<0>(db).X_Range_fv[0]) ||
          (rdr->vtx_X > std::get<0>(db).X_Range_fv[1]) ||
          (rdr->vtx_Y < std::get<0>(db).Y_Range_fv[0]) ||
          (rdr->vtx_Y > std::get<0>(db).Y_Range_fv[1]) ||
          (rdr->vtx_Z < std::get<0>(db).Z_Range_fv[0]) ||
          (rdr->vtx_Z > std::get<0>(db).Z_Range_fv[1])) {
        continue;
      }

      //Checking if lepton exits
      double Detlow = std::get<0>(Detectors[d_it]).XOffset -
                      std::get<0>(Detectors[d_it]).XWidth_det / 2.0;
      double DetHigh = std::get<0>(Detectors[d_it]).XOffset +
                       std::get<0>(Detectors[d_it]).XWidth_det / 2.0;
      //Exits through the side. Should supercede exiting back/front/y
      if(rdr->lepExitingPosX < Detlow){
        std::get<1>(db).flagLepExitBack = false;
        std::get<1>(db).flagLepExitFront = false;
        std::get<1>(db).flagLepExitY = false;
        std::get<1>(db).flagLepExitXLow = true;
        std::get<1>(db).flagLepExitXHigh = false;
      }
      else if(rdr->lepExitingPosX > DetHigh){
        std::get<1>(db).flagLepExitBack = false;
        std::get<1>(db).flagLepExitFront = false;
        std::get<1>(db).flagLepExitY = false;
        std::get<1>(db).flagLepExitXHigh = true;        
        std::get<1>(db).flagLepExitXLow = false;        
      }
      else{
        if(std::get<1>(Detectors[d_it]).flagLepExitBack && std::get<1>(Detectors[d_it]).flagLepExitY){
          if(rdr->lepExitingPosZ - std::get<0>(Detectors[d_it]).ZWidth_det / 2.0 > 
             fabs(rdr->lepExitingPosY) - std::get<0>(Detectors[d_it]).YWidth_det / 2.0){
            std::get<1>(db).flagLepExitBack = true;
            std::get<1>(db).flagLepExitFront = rdr->flagLepExitFront;
            std::get<1>(db).flagLepExitY = false; 
          }
          else{
            std::get<1>(db).flagLepExitBack = false;
            std::get<1>(db).flagLepExitFront = rdr->flagLepExitFront;
            std::get<1>(db).flagLepExitY = true; 
          }
        }
        else if(std::get<1>(Detectors[d_it]).flagLepExitFront && std::get<1>(Detectors[d_it]).flagLepExitY){
          if(fabs(rdr->lepExitingPosZ) - std::get<0>(Detectors[d_it]).ZWidth_det / 2.0 > 
             fabs(rdr->lepExitingPosY) - std::get<0>(Detectors[d_it]).YWidth_det / 2.0){
            std::get<1>(db).flagLepExitFront = true;
            std::get<1>(db).flagLepExitBack = rdr->flagLepExitBack;
            std::get<1>(db).flagLepExitY = false; 
          }
          else{
            std::get<1>(db).flagLepExitFront = false;
            std::get<1>(db).flagLepExitBack = rdr->flagLepExitBack;
            std::get<1>(db).flagLepExitY = true; 
          }
        }
        else{
          std::get<1>(db).flagLepExitBack = rdr->flagLepExitBack;
          std::get<1>(db).flagLepExitFront = rdr->flagLepExitFront;
          std::get<1>(db).flagLepExitY = rdr->flagLepExitY;
        }
        std::get<1>(db).flagLepExitXHigh = false; 
        std::get<1>(db).flagLepExitXLow = false;                
      }
      std::get<1>(db).lepExitingMomX = rdr->lepExitingMomX;
      std::get<1>(db).lepExitingMomY = rdr->lepExitingMomY;
      std::get<1>(db).lepExitingMomZ = rdr->lepExitingMomZ;
      ////////End exiting lepton section
     
      //Truth info 
      std::get<1>(db).Q2True = rdr->Q2True;
      std::get<1>(db).yTrue = rdr->yTrue;
      std::get<1>(db).W_rest = rdr->W_rest;
      /////////
      
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
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eProtonPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eProtonPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true) +
          Jaccumulate(rdr->eProtonSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true);

      std::get<1>(db).NeutronDep_FV =
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->eNeutronSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).NeutronDep_veto =
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eNeutronSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eNeutronSecondaryDep,
                      std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eNeutronPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true) +
          Jaccumulate(rdr->eNeutronSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true);

      std::get<1>(db).PiCDep_FV =
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).PiCDep_veto =
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->ePiCPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true) +
          Jaccumulate(rdr->ePiCSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true);

      std::get<1>(db).Pi0Dep_FV =
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).Pi0Dep_veto =
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->ePi0PrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true) +
          Jaccumulate(rdr->ePi0SecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true);

      std::get<1>(db).OtherDep_FV =
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0);
      std::get<1>(db).OtherDep_veto =
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_veto_left[0],
                      std::get<0>(db).X_veto_left[1], 0, true, true) +
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_veto_right[0],
                      std::get<0>(db).X_veto_right[1], 0, true, true) +
          Jaccumulate(rdr->eOtherPrimaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true) +
          Jaccumulate(rdr->eOtherSecondaryDep, std::get<0>(db).X_fv[0],
                      std::get<0>(db).X_fv[1], 0, false, true);

      // std::cout << "\nEnu: " << &std::get<1>(Detectors[d_it]).Enu << ", " <<
      // std::get<1>(Detectors[d_it]).Enu << std::endl;

      std::get<1>(db).Enu = rdr->Enu;

      // std::cout << "\nEnu: " << &std::get<1>(Detectors[d_it]).Enu << ", " <<
      // std::get<1>(Detectors[d_it]).Enu << std::endl;

      std::get<1>(db).NuPDG = rdr->nuPDG;
      std::get<1>(db).LepPDG = rdr->lepPDG;
      std::get<1>(db).TotalDep_FV =
          std::get<1>(db).ProtonDep_FV + std::get<1>(db).NeutronDep_FV +
          std::get<1>(db).PiCDep_FV + std::get<1>(db).Pi0Dep_FV +
          std::get<1>(db).OtherDep_FV;

      std::get<1>(db).TotalDep_veto =
          std::get<1>(db).ProtonDep_veto + std::get<1>(db).NeutronDep_veto +
          std::get<1>(db).PiCDep_veto + /*std::get<1>(db).Pi0Dep_veto +*/
          std::get<1>(db).OtherDep_veto;

      if (((std::get<1>(db).TotalDep_veto != 0) &&
           (!std::isnormal(std::get<1>(db).TotalDep_veto))) ||
          ((std::get<1>(db).TotalDep_FV != 0) &&
           (!std::isnormal(std::get<1>(db).TotalDep_FV)))) {
        std::cout << "\n[INFO] XBins: " << std::get<0>(db).X_fv[0] << " -- "
                  << std::get<0>(db).X_fv[1]
                  << ", Veto left: " << std::get<0>(db).X_veto_left[0] << " -- "
                  << std::get<0>(db).X_veto_left[1]
                  << ", Veto right: " << std::get<0>(db).X_veto_right[0]
                  << " -- " << std::get<0>(db).X_veto_right[1] << std::endl;

        std::cout << "FV -- Total: " << std::get<1>(db).TotalDep_FV
                  << ", Proton: " << std::get<1>(db).ProtonDep_FV
                  << ", Neutron: " << std::get<1>(db).NeutronDep_FV
                  << ", PiC: " << std::get<1>(db).PiCDep_FV
                  << ", Pi0: " << std::get<1>(db).Pi0Dep_FV
                  << ", Other: " << std::get<1>(db).OtherDep_FV << std::endl;

        std::cout << "Veto -- Total: " << std::get<1>(db).TotalDep_veto
                  << ", Proton: " << std::get<1>(db).ProtonDep_veto
                  << ", Neutron: " << std::get<1>(db).NeutronDep_veto
                  << ", PiC: " << std::get<1>(db).PiCDep_veto
                  << ", Pi0: " << std::get<1>(db).Pi0Dep_veto
                  << ", Other: " << std::get<1>(db).OtherDep_veto << std::endl
                  << std::endl;
      }


      std::get<2>(db)->Fill();
    }

    // break;
  }

  std::cout << "\r[INFO]: Read " << rdr->GetEntries() << " entries."
            << std::endl;

  // for (std::tuple<DetBox, EDep, TTree *> &db : Detectors) {
  //   std::get<2>(db)->Write(std::get<2>(db)->GetName(), TObject::kOverwrite);
  // }

  outfile->Write();
  outfile->Close();
}
