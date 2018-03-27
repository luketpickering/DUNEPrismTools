#include "G4ArReader.h"
#include "Utils.hxx"

#include "TApplication.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TProfile2D.h"

std::vector<double> detmin;
std::vector<double> detmax;
std::vector<double> fvgap;
std::vector<int> ndetsteps = {400, 0, 0};
std::string inputG4ArFileName;
int ent = 0;

void SayUsage(char *argv[]) {
  std::cout << "[INFO]: Use like: " << argv[0]
            << " -i <inputg4arbofile> -nx <nxsteps> -dmn <detxmin,ymin,zmin> "
               "-dmx <detxmax,ymax,zmax> -fv <fidgapx,y,z> [-e <entry>]"
            << std::endl;
}

void handleOpts(int argc, char *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-n") {
      ndetsteps = ParseToVect<int>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-dmn") {
      detmin = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-dmx") {
      detmax = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-fv") {
      fvgap = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-i") {
      inputG4ArFileName = argv[++opt];
    } else if (std::string(argv[opt]) == "-e") {
      ent = str2T<int>(argv[++opt]);
    } else if ((std::string(argv[opt]) == "-?") ||
               std::string(argv[opt]) == "--help") {
      SayUsage(argv);
      exit(0);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char *argv[]) {
  handleOpts(argc, argv);

  DetectorAndFVDimensions detdims;
  detdims.NXSteps = ndetsteps[0];
  detdims.NYSteps = ndetsteps[1];
  detdims.NZSteps = ndetsteps[2];
  for (size_t dim_it = 0; dim_it < 3; ++dim_it) {
    detdims.DetMin[dim_it] = detmin[dim_it];
    detdims.DetMax[dim_it] = detmax[dim_it];
    detdims.FVGap[dim_it] = fvgap[dim_it];
  }

  G4ArReader g4ar(inputG4ArFileName, detdims);

  g4ar.GetEvent(ent);

  TH3D *MuonEDep = detdims.BuildDetectorMap();
  TH3D *MuonEDep_sec = detdims.BuildDetectorMap();
  TH3D *OtherEDep = detdims.BuildDetectorMap();
  TH3D *OtherEDep_sec = detdims.BuildDetectorMap();

  Event ev = g4ar.BuildEvent();

  for (DepoTracked &td : ev.TrackedDeposits) {
    if (abs(td.PDG) == 13) {
      MuonEDep->Add(td.Deposits);
      MuonEDep_sec->Add(td.DaughterDeposits);
    }
#ifdef DEBUG

    for (Int_t x_it = 0; x_it < detdims.NXSteps; ++x_it) {
      for (Int_t y_it = 0; y_it < detdims.NYSteps; ++y_it) {
        for (Int_t z_it = 0; z_it < detdims.NZSteps; ++z_it) {
          Int_t gbin = td.Deposits->GetBin(x_it + 1, y_it + 1, z_it + 1);

          if (td.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", td.Deposits content = "
                      << td.Deposits->GetBinContent(gbin) << std::endl;
          }

          if (td.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", td.DaughterDeposits content = "
                      << td.DaughterDeposits->GetBinContent(gbin) << std::endl;
          }
        }
      }
    }
#endif
  }

  for (DepoParticle &td : ev.TotalDeposits) {
    if (abs(td.PDG) == 13) {
      MuonEDep->Add(td.Deposits);
      MuonEDep_sec->Add(td.DaughterDeposits);
    } else {
      OtherEDep->Add(td.Deposits);
      OtherEDep_sec->Add(td.DaughterDeposits);
    }
#ifdef DEBUG

    for (Int_t x_it = 0; x_it < detdims.NXSteps; ++x_it) {
      for (Int_t y_it = 0; y_it < detdims.NYSteps; ++y_it) {
        for (Int_t z_it = 0; z_it < detdims.NZSteps; ++z_it) {
          Int_t gbin = td.Deposits->GetBin(x_it + 1, y_it + 1, z_it + 1);

          if (td.Deposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", td.Deposits content = "
                      << td.Deposits->GetBinContent(gbin) << std::endl;
          }

          if (td.DaughterDeposits->GetBinContent(gbin)) {
            std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                      << ", td.DaughterDeposits content = "
                      << td.DaughterDeposits->GetBinContent(gbin) << std::endl;
          }
        }
      }
    }
#endif
  }

#ifdef DEBUG
  for (Int_t x_it = 0; x_it < detdims.NXSteps; ++x_it) {
    for (Int_t y_it = 0; y_it < detdims.NYSteps; ++y_it) {
      for (Int_t z_it = 0; z_it < detdims.NZSteps; ++z_it) {
        Int_t gbin = MuonEDep->GetBin(x_it + 1, y_it + 1, z_it + 1);

        if (MuonEDep->GetBinContent(gbin)) {
          std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                    << ", MuonEDep content = " << MuonEDep->GetBinContent(gbin)
                    << std::endl;
        }

        if (MuonEDep_sec->GetBinContent(gbin)) {
          std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                    << ", MuonEDep_sec content = "
                    << MuonEDep_sec->GetBinContent(gbin) << std::endl;
        }

        if (OtherEDep->GetBinContent(gbin)) {
          std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                    << ", OtherEDep content = "
                    << OtherEDep->GetBinContent(gbin) << std::endl;
        }

        if (OtherEDep_sec->GetBinContent(gbin)) {
          std::cout << "[INFO]: Bin " << x_it << ", " << y_it << ", " << z_it
                    << ", OtherEDep_sec content = "
                    << OtherEDep_sec->GetBinContent(gbin) << std::endl;
        }
      }
    }
  }
#endif

  TApplication evDispLoop("App", &argc, argv);

  std::cout << "[INFO]: Drawing interaction at {" << ev.VertexPosition.X()
            << ", " << ev.VertexPosition.Y() << ", " << ev.VertexPosition.Z()
            << "}. " << std::endl;

  TCanvas *cmu_yx = new TCanvas("cmu_yx");
  cmu_yx->cd();

  MuonEDep->SetName("mu_e_dep");
  MuonEDep->SetTitle("Muon energy deposit");
  MuonEDep->Project3D("yx")->Draw("COLZ");
  TMarker mmu_yx(ev.VertexPosition.X(), ev.VertexPosition.Y(), 20);
  mmu_yx.Draw();
  cmu_yx->Update();

  TCanvas *cmu_yz = new TCanvas("cmu_yz");
  cmu_yz->cd();

  MuonEDep->Project3D("yz")->Draw("COLZ");
  TMarker mmu_yz(ev.VertexPosition.Z(), ev.VertexPosition.Y(), 20);
  mmu_yz.Draw();
  cmu_yz->Update();

  TCanvas *cother_yx = new TCanvas("cother_yx");
  cother_yx->cd();

  OtherEDep->SetName("other_e_dep");
  OtherEDep->SetTitle("Other energy deposits");
  OtherEDep->Project3D("yx")->Draw("COLZ");
  TMarker mother_yx(ev.VertexPosition.X(), ev.VertexPosition.Y(), 20);
  mother_yx.Draw();
  cother_yx->Update();

  TCanvas *cother_yz = new TCanvas("cother_yz");
  cother_yz->cd();

  OtherEDep->Project3D("yz")->Draw("COLZ");
  TMarker mother_yz(ev.VertexPosition.Z(), ev.VertexPosition.Y(), 20);
  mother_yz.Draw();
  cother_yz->Update();

  evDispLoop.Run();
  return 0;
}
