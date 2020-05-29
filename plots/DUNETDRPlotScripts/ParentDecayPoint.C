#include "TChain.h"
#include "TFile.h"

#include "TGraph.h"
#include "TH2.h"

#include <iostream>

struct NuRay {
  float phase_space_weight;
  float imp_weight;
  short pdg;
  float four_mom[4];
  float ray_end_point[3];

  short parent_pdg;
  float parent_position[3];
};

void ParentDecayPoint(char const *fname) {

  TChain nur("nurays");
  if (!nur.Add(fname)) {
    std::cout << "[ERROR]: Failed to add any files to the chain." << std::endl;
    abort();
  }

  NuRay nr;
  nur.SetBranchAddress("phase_space_weight", &nr.phase_space_weight);
  nur.SetBranchAddress("imp_weight", &nr.imp_weight);
  nur.SetBranchAddress("pdg", &nr.pdg);
  nur.SetBranchAddress("four_mom", &nr.four_mom);
  nur.SetBranchAddress("ray_end_point", &nr.ray_end_point);
  nur.SetBranchAddress("parent_pdg", &nr.parent_pdg);
  nur.SetBranchAddress("parent_position", &nr.parent_position);

  TFile fout("ParentDecayPointPlots.root", "RECREATE");

  TH2D *EOffAxis_unweighted =
      new TH2D("FluxShape2D_unweighted",
               ";E_{#nu} (GeV);Off axis position (m);#Phi_{#nu_{#mu}}", 200, 0,
               8, 105, -2, 33);

  TH2D *EOffAxis_phase_space_weighted =
      new TH2D("FluxShape2D_phase_space_weighted",
               ";E_{#nu} (GeV);Off axis position (m);#Phi_{#nu_{#mu}}", 200, 0,
               8, 105, -2, 33);

  TH2D *EOffAxis = new TH2D(
      "FluxShape2D", ";E_{#nu} (GeV);Off axis position (m);#Phi_{#nu_{#mu}}",
      200, 0, 8, 105, -2, 33);

  TH2D *PDPZ = new TH2D(
      "PDPZ",
      ";Mean Parent Decay Point (m);Off axis position (m);#Phi_{#nu_{#mu}}",
      100, -10, 240, 14, -2, 33);

  for (long i = 0; i < nur.GetEntries(); ++i) {
    nur.GetEntry(i);

    if (!(i % 500000)) {
      std::cout << "[INFO]: Reading entry " << i << " ("
                << (i * 100 / nur.GetEntries()) << " %)" << std::endl;
    }

    EOffAxis_unweighted->Fill(nr.four_mom[3], nr.ray_end_point[0] * 1E-2);
    EOffAxis_phase_space_weighted->Fill(
        nr.four_mom[3], nr.ray_end_point[0] * 1E-2, nr.phase_space_weight);
    EOffAxis->Fill(nr.four_mom[3], nr.ray_end_point[0] * 1E-2,
                   nr.phase_space_weight * nr.imp_weight);
    PDPZ->Fill(nr.parent_position[2] * 1E-2, nr.ray_end_point[0] * 1E-2,
               nr.phase_space_weight * nr.imp_weight);
  }

  TGraph meanpdp(14);
  TGraph zspreadpdp(14);

  for (int i = 0; i < 14; ++i) {
  	auto px= PDPZ->ProjectionX("_px", i + 1, i + 1);
    double m =px->GetMean();
    double w =px->GetStdDev();

    meanpdp.SetPoint(i, -2.0 + (35.0 / 14.0) * double(i), m);
    zspreadpdp.SetPoint(i, -2.0 + (35.0 / 14.0) * double(i), w);

    std::cout << i << " " << -2.0 + (35.0 / 14.0) * double(i) << " " << m << std::endl;
  }

  meanpdp.SetTitle(
      ";Off axis position (m); Mean Parent Decay point Z (m, beam coordiates)");

  fout.WriteObject(&meanpdp, "meanpdp");

  zspreadpdp.SetTitle(
      ";Off axis position (m); StdDev Parent Decay point Z (m, beam coordiates)");

  fout.WriteObject(&zspreadpdp, "zspreadpdp");


  fout.Write();
  fout.Close();
}