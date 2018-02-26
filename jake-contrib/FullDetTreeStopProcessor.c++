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
  int stop;

  double vtx[3];
  TObjString * eventCode = new TObjString();
  double Enu;
  double yTrue;
  double Q2True;
  TLorentzVector * qTrue;
  double W_rest;

  int NuPDG;
  int LepPDG;

  int nPi0;
  int nPiC;
  int nProton;
  int nNeutron;
  int nGamma;

  double eLepTrue;
  double ePi0True;
  double ePiCTrue;
  double eProtonTrue;
  double eNeutronTrue;
  double eGammaTrue;
  double TotalNonlep_eTrue;

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

  double TotalNonlep_Dep_FV;
  double TotalNonlep_Dep_veto;

  bool flagLepExitBack;
  bool flagLepExitFront;
  bool flagLepExitYLow;
  bool flagLepExitYHigh;
  bool flagLepExit;
  bool flagLepExitXLow;
  bool flagLepExitXHigh;

  double lepExitingPosX;
  double lepExitingPosY;
  double lepExitingPosZ;
  double lepExitingMomX;
  double lepExitingMomY;
  double lepExitingMomZ;

  double eElectronShowerDepInside;
  double eElectronShowerDepOutside;
};

std::vector<DetectorStop> DetectorStops;
std::vector<DetBox> Detectors;

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
               "\t-o <output.root>              : Output file name.\n"
               "\t-S <shift_cm>                 : Detector XShift in cm "
               "(default = -1800).\n"
               "\t-A <LAr length cm>            : Full span of detector (FV + "
               "veto) range (default = 4000).\n"
               "\t-V <Veto width cm>            : Veto region in each "
               "dimension, on each side of the FV (default = 50).\n"
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

    } else if (std::string(argv[opt]) == "-S") {
      XShift_cm = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-A") {
      XArLen_cm = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-V") {
      X_Veto_cm = str2T<double>(argv[++opt]);
      Y_Veto_cm = X_Veto_cm;
      Z_Veto_cm = X_Veto_cm;
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

  while (true) {
    if (min == max) {
      std::cout << "[ERROR]: Couldn't find " << find
                << " but should have been able to. Last check was " << min
                << std::endl;
      throw;
    }

    int bi = min + floor((max - min) / 2);

    bool eq = (fabs(find - BinEdges[bi - 1]) < 1E-4);
    if (eq) {
      return bi - 1;
    }

    bool lt = (find < BinEdges[bi - 1]);
    if (lt) {
      max = bi;
      continue;
    }

    bool ib = (find < BinEdges[bi]);
    if (ib) {
      return bi - 1;
    }

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

  std::vector<double> XBins;

  XBins.push_back(XShift_cm - XArLen_cm / 2.0);
  for (size_t bi_it = 1; bi_it < 400; ++bi_it) {
    XBins.push_back(XBins.back() + XBinWidth_cm);
  }
  XBins.push_back(XBins.back() + XBinWidth_cm);

  size_t NDets = DetectorStops.size();
  Detectors.reserve(NDets);

  EDep OutputEDep;
  TTree *OutputTree = new TTree("EDeps", "");

  // translate to detboxes
  for (size_t d_it = 0; d_it < NDets; ++d_it) {
    DetBox db;

    db.XOffset = -DetectorStops[d_it].LateralOffset * 100.0;
    db.XWidth_fv = DetectorStops[d_it].DetectorFiducialWidth * 100.0;
    db.YWidth_fv = DetectorStops[d_it].DetectorFiducialHeight * 100.0;
    db.ZWidth_fv = DetectorStops[d_it].DetectorFiducialDepth * 100.0;
    db.XWidth_det = db.XWidth_fv + 2 * X_Veto_cm;
    db.YWidth_det = db.YWidth_fv + 2 * Y_Veto_cm;
    db.ZWidth_det = db.ZWidth_fv + 2 * Z_Veto_cm;

    double Detlow = db.XOffset - db.XWidth_det / 2.0;
    double DetHigh = db.XOffset + db.XWidth_det / 2.0;

    db.X_Range_fv[0] = db.XOffset - db.XWidth_fv / 2.0;
    db.X_Range_fv[1] = db.XOffset + db.XWidth_fv / 2.0;

    db.Y_Range_fv[0] = -db.YWidth_fv / 2.0;
    db.Y_Range_fv[1] = db.YWidth_fv / 2.0;

    db.Z_Range_fv[0] = -db.ZWidth_fv / 2.0;
    db.Z_Range_fv[1] = db.ZWidth_fv / 2.0;

    std::cout << "[INFO]: Det: [" << Detlow << ", " << DetHigh << "], FV: ["
              << db.X_Range_fv[0] << ", " << db.X_Range_fv[1] << "]."
              << std::endl;

    db.X_fv[0] = BinSearch(XBins, db.X_Range_fv[0]);
    db.X_fv[1] = BinSearch(XBins, db.X_Range_fv[1]);

    db.X_veto_left[0] = BinSearch(XBins, Detlow);
    db.X_veto_left[1] = BinSearch(XBins, Detlow + X_Veto_cm);

    db.X_veto_right[0] = BinSearch(XBins, DetHigh - X_Veto_cm);
    db.X_veto_right[1] = BinSearch(XBins, DetHigh);

    Detectors.push_back(db);
  }

  //Init qTrue 
  OutputEDep.qTrue = new TLorentzVector();

  OutputTree->Branch("stop", &OutputEDep.stop, "stop/I");

  OutputTree->Branch("vtx", &OutputEDep.vtx, "vtx[3]/D");
  OutputTree->Branch("eventCode", &OutputEDep.eventCode);

  OutputTree->Branch("Enu", &OutputEDep.Enu, "Enu/D");
  OutputTree->Branch("yTrue", &OutputEDep.yTrue, "yTrue/D");
  OutputTree->Branch("W_rest", &OutputEDep.W_rest, "W_rest/D");
  OutputTree->Branch("Q2True", &OutputEDep.Q2True, "Q2True/D");
  OutputTree->Branch("qTrue", &OutputEDep.qTrue);

  OutputTree->Branch("NuPDG", &OutputEDep.NuPDG, "NuPDG/I");
  OutputTree->Branch("LepPDG", &OutputEDep.LepPDG, "LepPDG/I");

  OutputTree->Branch("nPi0", &OutputEDep.nPi0, "nPi0/I");
  OutputTree->Branch("nPiC", &OutputEDep.nPiC, "nPiC/I");
  OutputTree->Branch("nProton", &OutputEDep.nProton, "nProton/I");
  OutputTree->Branch("nNeutron", &OutputEDep.nNeutron, "nNeutron/I");
  OutputTree->Branch("nGamma", &OutputEDep.nGamma, "nGamma/I");

  OutputTree->Branch("eLepTrue", &OutputEDep.eLepTrue, "eLepTrue/D");
  OutputTree->Branch("ePi0True", &OutputEDep.ePi0True, "ePi0True/D");
  OutputTree->Branch("ePiCTrue", &OutputEDep.ePiCTrue, "ePiCTrue/D");
  OutputTree->Branch("eProtonTrue", &OutputEDep.eProtonTrue, "eProtonTrue/D");
  OutputTree->Branch("eNeutronTrue", &OutputEDep.eNeutronTrue,
                     "eNeutronTrue/D");
  OutputTree->Branch("eGammaTrue", &OutputEDep.eGammaTrue, "eGammaTrue/D");

  OutputTree->Branch("TotalNonlep_eTrue", &OutputEDep.TotalNonlep_eTrue,
                     "TotalNonlep_eTrue/D");

  OutputTree->Branch("LepDep_det", &OutputEDep.LepDep_det, "LepDep_det/D");
  OutputTree->Branch("ProtonDep_FV", &OutputEDep.ProtonDep_FV,
                     "ProtonDep_FV/D");
  OutputTree->Branch("ProtonDep_veto", &OutputEDep.ProtonDep_veto,
                     "ProtonDep_veto/D");
  OutputTree->Branch("NeutronDep_FV", &OutputEDep.NeutronDep_FV,
                     "NeutronDep_FV/D");
  OutputTree->Branch("NeutronDep_veto", &OutputEDep.NeutronDep_veto,
                     "NeutronDep_veto/D");
  OutputTree->Branch("PiCDep_FV", &OutputEDep.PiCDep_FV, "PiCDep_FV/D");
  OutputTree->Branch("PiCDep_veto", &OutputEDep.PiCDep_veto, "PiCDep_veto/D");
  OutputTree->Branch("Pi0Dep_FV", &OutputEDep.Pi0Dep_FV, "Pi0Dep_FV/D");
  OutputTree->Branch("Pi0Dep_veto", &OutputEDep.Pi0Dep_veto, "Pi0Dep_veto/D");
  OutputTree->Branch("OtherDep_FV", &OutputEDep.OtherDep_FV, "OtherDep_FV/D");
  OutputTree->Branch("OtherDep_veto", &OutputEDep.OtherDep_veto,
                     "OtherDep_veto/D");

  OutputTree->Branch("TotalNonlep_Dep_FV", &OutputEDep.TotalNonlep_Dep_FV,
                     "TotalNonlep_Dep_FV/D");
  OutputTree->Branch("TotalNonlep_Dep_veto", &OutputEDep.TotalNonlep_Dep_veto,
                     "TotalNonlep_Dep_veto/D");

  OutputTree->Branch("flagLepExitBack", &OutputEDep.flagLepExitBack,
                     "flagLepExitBack/B");
  OutputTree->Branch("flagLepExitFront", &OutputEDep.flagLepExitFront,
                     "flagLepExitFront/B");
  OutputTree->Branch("flagLepExitYLow", &OutputEDep.flagLepExitYLow,
                     "flagLepExitYLow/B");
  OutputTree->Branch("flagLepExitYHigh", &OutputEDep.flagLepExitYHigh,
                     "flagLepExitYHigh/B");
  OutputTree->Branch("flagLepExit", &OutputEDep.flagLepExit, "flagLepExit/B");

  OutputTree->Branch("flagLepExitXLow", &OutputEDep.flagLepExitXLow,
                     "flagLepExitXLow/B");
  OutputTree->Branch("flagLepExitXHigh", &OutputEDep.flagLepExitXHigh,
                     "flagLepExitXHigh/B");
  OutputTree->Branch("lepExitingPosX", &OutputEDep.lepExitingPosX,
                     "lepExitingPosX/D");
  OutputTree->Branch("lepExitingPosY", &OutputEDep.lepExitingPosY,
                     "lepExitingPosY/D");
  OutputTree->Branch("lepExitingPosZ", &OutputEDep.lepExitingPosZ,
                     "lepExitingPosZ/D");
  OutputTree->Branch("lepExitingMomX", &OutputEDep.lepExitingMomX,
                     "lepExitingMomX/D");
  OutputTree->Branch("lepExitingMomY", &OutputEDep.lepExitingMomY,
                     "lepExitingMomY/D");
  OutputTree->Branch("lepExitingMomZ", &OutputEDep.lepExitingMomZ,
                     "lepExitingMomZ/D");
  OutputTree->Branch("eElectronShowerDepInside", 
                     &OutputEDep.eElectronShowerDepInside,
                     "eElectronShowerDepInside/D");
  OutputTree->Branch("eElectronShowerDepOutside", 
                     &OutputEDep.eElectronShowerDepOutside,
                     "eElectronShowerDepOutside/D");


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

    int stop = -1;
    for (size_t d_it = 0; d_it < NDets; ++d_it) {
      DetBox &db = Detectors[d_it];

      if ((rdr->vtx_X < db.X_Range_fv[0]) || (rdr->vtx_X > db.X_Range_fv[1]) ||
          (rdr->vtx_Y < db.Y_Range_fv[0]) || (rdr->vtx_Y > db.Y_Range_fv[1]) ||
          (rdr->vtx_Z < db.Z_Range_fv[0]) || (rdr->vtx_Z > db.Z_Range_fv[1])) {
        continue;
      }
      stop = d_it;
      break;
    }

    if (stop == -1) {
      continue;
    }

    DetBox &stopBox = Detectors[stop];
    OutputEDep.stop = stop;
    OutputEDep.eventCode = rdr->eventCode;
    OutputEDep.vtx[0] = rdr->vtx_X;
    OutputEDep.vtx[1] = rdr->vtx_Y;
    OutputEDep.vtx[2] = rdr->vtx_Z;

    OutputEDep.Enu = rdr->Enu;
    
    OutputEDep.yTrue = rdr->yTrue;
    OutputEDep.Q2True = rdr->Q2True;
    OutputEDep.qTrue = (TLorentzVector*)rdr->qTrue->Clone();
    OutputEDep.W_rest = rdr->W_rest;

    OutputEDep.NuPDG = rdr->nuPDG;
    OutputEDep.LepPDG = rdr->lepPDG;

    OutputEDep.nPi0 = rdr->nPi0;
    OutputEDep.nPiC = rdr->nPiC;
    OutputEDep.nProton = rdr->nProton;
    OutputEDep.nNeutron = rdr->nNeutron;
    OutputEDep.nGamma = rdr->nGamma;
    OutputEDep.eLepTrue = rdr->eLepTrue;
    OutputEDep.ePi0True = rdr->ePi0True;
    OutputEDep.ePiCTrue = rdr->ePiCTrue;
    OutputEDep.eProtonTrue = rdr->eProtonTrue;
    OutputEDep.eNeutronTrue = rdr->eNeutronTrue;
    OutputEDep.eGammaTrue = rdr->eGammaTrue;
    OutputEDep.TotalNonlep_eTrue = rdr->ePi0True + rdr->ePiCTrue +
                                   rdr->eProtonTrue + rdr->eNeutronTrue +
                                   rdr->eGammaTrue;
    OutputEDep.eElectronShowerDepInside = 0.;
    OutputEDep.eElectronShowerDepOutside = 0.;
    double DetXLow  = stopBox.XOffset - stopBox.XWidth_det / 2.0;
    double DetXHigh = stopBox.XOffset + stopBox.XWidth_det / 2.0;
    double DetYLow  = 0. - stopBox.YWidth_det / 2.0;
    double DetYHigh = stopBox.YWidth_det / 2.0;
    double DetZLow  = 0. - stopBox.ZWidth_det / 2.0;
    double DetZHigh = stopBox.ZWidth_det / 2.0;

    // Checking if lepton exits -- muon only
    OutputEDep.flagLepExitBack = false;
    OutputEDep.flagLepExitFront = false;
    OutputEDep.flagLepExitYLow = false;
    OutputEDep.flagLepExitYHigh = false;
    OutputEDep.flagLepExitXLow = false;
    OutputEDep.flagLepExitXHigh = false; 
    OutputEDep.flagLepExit = false;
    if(abs(rdr->lepPDG) == 13){
      
      for(int i = 0; i < 1000; ++i){
         
        if(OutputEDep.flagLepExit){break;}

        double posX = rdr->lepTrackX[i];
        double posY = rdr->lepTrackY[i];
        double posZ = rdr->lepTrackZ[i];
        if(posX == 0 && posY == 0 && posZ == 0){break;} 
        double dX = 0., dY = 0., dZ = 0.;
        bool exitX = false, exitY = false, exitZ = false;

        if(posX <= DetXLow || posX >= DetXHigh){
          dX = fabs(DetXLow - posX);
          exitX = true;
        }
        else{
          exitX = false;       
        }

        if(posY <= DetYLow || posY >= DetYHigh){
          dY = fabs(DetYLow - posY);
          exitY = true;
        }
        else{
          exitY = false;       
        }

        if(posZ <= DetZLow || posZ >= DetZHigh){
          dZ = fabs(DetZLow - posZ);
          exitZ = true;
        }
        else{
          exitZ = false;       
        }

        if(exitX || exitY || exitZ){
        //std::cout << "Det: " << stop << " " << DetXHigh << " " << DetXLow << " " << DetYHigh << " " << DetYLow << " " << DetZHigh << " " << DetZLow << std::endl;
       // std::cout  << i << " "<< posX << " "<< posY << " " << posZ << std::endl;
//        std::cout  << "\t"<< i-1 << " "<< rdr->lepTrackX[i-1] << " "<< rdr->lepTrackY[i-1] << " " << rdr->lepTrackZ[i-1] << std::endl;

          if(dX > dY && dX > dZ){
            OutputEDep.flagLepExitBack = false;
            OutputEDep.flagLepExitFront = false;
            OutputEDep.flagLepExitYLow = false;
            OutputEDep.flagLepExitYHigh = false;

            if(posX <= DetXLow){
              OutputEDep.flagLepExitXLow = true;
              OutputEDep.flagLepExitXHigh = false; 
            }
            else if(posX >= DetXHigh){
              OutputEDep.flagLepExitXLow = false;
              OutputEDep.flagLepExitXHigh = true; 
            }
            else{
              std::cout << "ERROR X" << std::endl;           
            }
          }
          else if(dY > dX && dY > dZ){
            OutputEDep.flagLepExitBack = false;
            OutputEDep.flagLepExitFront = false;
            OutputEDep.flagLepExitXLow = false;
            OutputEDep.flagLepExitXHigh = false;

            if(posY <= DetYLow){
              OutputEDep.flagLepExitYLow = true;
              OutputEDep.flagLepExitYHigh = false; 
            }
            else if(posY >= DetYHigh){
              OutputEDep.flagLepExitYLow = false;
              OutputEDep.flagLepExitYHigh = true; 
            }
            else{
              std::cout << "ERROR Y" << std::endl;           
            }
          }
          else if(dZ > dX && dZ > dY){
            OutputEDep.flagLepExitXLow = false;
            OutputEDep.flagLepExitXHigh = false;
            OutputEDep.flagLepExitYLow = false;
            OutputEDep.flagLepExitYHigh = false;

            if(posZ <= DetZLow){
              OutputEDep.flagLepExitBack = false;
              OutputEDep.flagLepExitFront = true;
            }
            else if(posZ >= DetZHigh){
              OutputEDep.flagLepExitFront = false;
              OutputEDep.flagLepExitBack = true; 
            }
            else{
              std::cout << "ERROR Z" << std::endl;           
            }
          }
          OutputEDep.flagLepExit = true;
          OutputEDep.lepExitingMomX = rdr->lepTrackMomX[i];
          OutputEDep.lepExitingMomY = rdr->lepTrackMomY[i];
          OutputEDep.lepExitingMomZ = rdr->lepTrackMomZ[i];
          OutputEDep.lepExitingPosX = rdr->lepTrackX[i];
          OutputEDep.lepExitingPosY = rdr->lepTrackY[i];
          OutputEDep.lepExitingPosZ = rdr->lepTrackZ[i];
          break;
        }
      }
    }
    ////////End exiting lepton section
    //
    //
    //
    //Splitting up electron shower deposits
    else if(abs(rdr->lepPDG) == 11){
      for(int i = 0; i < 400; ++i){
         if(rdr->eElectronShowerDepInside[i] > 0.){
           double posLow = -3800. + i*10.; 
           double posHigh = -3800. + (i+1)*10.;

           //if (posLow < DetXLow || posHigh > DetXHigh){
           if ( (posLow < stopBox.X_Range_fv[0] && posLow >= DetXLow )     //Between edge of detector and 
             || (posHigh > stopBox.X_Range_fv[1] && posHigh <= DetXHigh ) ){ //fiducial region -- IN VETO
             OutputEDep.eElectronShowerDepOutside += rdr->eElectronShowerDepInside[i];
           /*  std::cout << "OUTSIDE" << std::endl;
             std::cout << "Low: " << DetXLow << " " << posLow << " " << stopBox.X_Range_fv[0] << std::endl;
             std::cout << "High: " << stopBox.X_Range_fv[1] << " " << posHigh << " " << DetXHigh << std::endl;*/
           }
           else{
             OutputEDep.eElectronShowerDepInside += rdr->eElectronShowerDepInside[i];
            /* std::cout << "INSIDE" << std::endl;
             std::cout << "Low: " << DetXLow << " " << posLow << " " << stopBox.X_Range_fv[0] << std::endl;
             std::cout << "High: " << stopBox.X_Range_fv[1] << " " << posHigh << " " << DetXHigh << std::endl;*/
           }
         }
      }
      OutputEDep.eElectronShowerDepOutside += rdr->eElectronShowerDepOutside;
    }
    ///////End electron shower deposits
    
        
    OutputEDep.LepDep_det =
        Jaccumulate(rdr->eLepPrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_right[1], 0) +
        Jaccumulate(rdr->eLepSecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_right[1], 0);

    OutputEDep.ProtonDep_FV = Jaccumulate(rdr->eProtonPrimaryDep,
                                          stopBox.X_fv[0], stopBox.X_fv[1], 0) +
                              Jaccumulate(rdr->eProtonSecondaryDep,
                                          stopBox.X_fv[0], stopBox.X_fv[1], 0);
    OutputEDep.ProtonDep_veto =
        Jaccumulate(rdr->eProtonPrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eProtonSecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eProtonPrimaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eProtonSecondaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eProtonPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true) +
        Jaccumulate(rdr->eProtonSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, false, true);

    OutputEDep.NeutronDep_FV =
        Jaccumulate(rdr->eNeutronPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0) +
        Jaccumulate(rdr->eNeutronSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0);
    OutputEDep.NeutronDep_veto =
        Jaccumulate(rdr->eNeutronPrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eNeutronSecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eNeutronPrimaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eNeutronSecondaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eNeutronPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, false, true) +
        Jaccumulate(rdr->eNeutronSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, false, true);

    OutputEDep.PiCDep_FV =
        Jaccumulate(rdr->ePiCPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0) +
        Jaccumulate(rdr->ePiCSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0);
    OutputEDep.PiCDep_veto =
        Jaccumulate(rdr->ePiCPrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->ePiCSecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->ePiCPrimaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->ePiCSecondaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->ePiCPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true) +
        Jaccumulate(rdr->ePiCSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true);

    OutputEDep.Pi0Dep_FV =
        Jaccumulate(rdr->ePi0PrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0) +
        Jaccumulate(rdr->ePi0SecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0);
    OutputEDep.Pi0Dep_veto =
        Jaccumulate(rdr->ePi0PrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->ePi0SecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->ePi0PrimaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->ePi0SecondaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->ePi0PrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true) +
        Jaccumulate(rdr->ePi0SecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true);

    OutputEDep.OtherDep_FV = Jaccumulate(rdr->eOtherPrimaryDep, stopBox.X_fv[0],
                                         stopBox.X_fv[1], 0) +
                             Jaccumulate(rdr->eOtherSecondaryDep,
                                         stopBox.X_fv[0], stopBox.X_fv[1], 0);
    OutputEDep.OtherDep_veto =
        Jaccumulate(rdr->eOtherPrimaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eOtherSecondaryDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, true, true) +
        Jaccumulate(rdr->eOtherPrimaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eOtherSecondaryDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, true, true) +
        Jaccumulate(rdr->eOtherPrimaryDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    false, true) +
        Jaccumulate(rdr->eOtherSecondaryDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, false, true);

    OutputEDep.TotalNonlep_Dep_FV =
        OutputEDep.ProtonDep_FV + OutputEDep.NeutronDep_FV +
        OutputEDep.PiCDep_FV + OutputEDep.Pi0Dep_FV + OutputEDep.OtherDep_FV;

    OutputEDep.TotalNonlep_Dep_veto =
        OutputEDep.ProtonDep_veto + OutputEDep.NeutronDep_veto +
        OutputEDep.PiCDep_veto + OutputEDep.Pi0Dep_veto +
        OutputEDep.OtherDep_veto;

    if (((OutputEDep.TotalNonlep_Dep_veto != 0) &&
         (!std::isnormal(OutputEDep.TotalNonlep_Dep_veto))) ||
        ((OutputEDep.TotalNonlep_Dep_FV != 0) &&
         (!std::isnormal(OutputEDep.TotalNonlep_Dep_FV)))) {
      std::cout << "\n[INFO] XBins: " << stopBox.X_fv[0] << " -- "
                << stopBox.X_fv[1] << ", Veto left: " << stopBox.X_veto_left[0]
                << " -- " << stopBox.X_veto_left[1]
                << ", Veto right: " << stopBox.X_veto_right[0] << " -- "
                << stopBox.X_veto_right[1] << std::endl;

      std::cout << "FV -- Total: " << OutputEDep.TotalNonlep_Dep_FV
                << ", Proton: " << OutputEDep.ProtonDep_FV
                << ", Neutron: " << OutputEDep.NeutronDep_FV
                << ", PiC: " << OutputEDep.PiCDep_FV
                << ", Pi0: " << OutputEDep.Pi0Dep_FV
                << ", Other: " << OutputEDep.OtherDep_FV << std::endl;

      std::cout << "Veto -- Total: " << OutputEDep.TotalNonlep_Dep_veto
                << ", Proton: " << OutputEDep.ProtonDep_veto
                << ", Neutron: " << OutputEDep.NeutronDep_veto
                << ", PiC: " << OutputEDep.PiCDep_veto
                << ", Pi0: " << OutputEDep.Pi0Dep_veto
                << ", Other: " << OutputEDep.OtherDep_veto << std::endl
                << std::endl;
    }

    OutputTree->Fill();
  }

  std::cout << "\r[INFO]: Read " << rdr->GetEntries() << " entries."
            << std::endl;

  outfile->Write();
  outfile->Close();
}
