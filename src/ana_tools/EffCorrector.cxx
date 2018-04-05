#include "EffCorrector.hxx"

EffCorrector::EffCorrector(EffCorrector::ModeEnum mode,
  std::string const &EffCorrectorFile,
  StopConfig &sc){

  this->mode = mode;

  MuonKinematics_seleff =
    GetHistogram<TH2>(EffCorrectorFile, "MuonKinematics_eff");
  ShowerKinematics_seleff = nullptr;

  switch (mode) {
    case kEHadrAbsPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_abspos_seleff");
      break;
    }
    case kEHadrDetPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_detpos_seleff");
      break;
    }
    case kEHadrVisDetPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "EVisHadr_FV_detpos_seleff");
      break;
    }
  }

  StopActiveRegions = sc.GetStopBoundingBoxes(false);
}

double EffCorrector::GetMuonKinematicsEffWeight(DepositsSummary const & edr){

  TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
  TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                 edr.PrimaryLep_4mom[2]);
  TrDir = TrDir.Unit();

  double MuToWall = CalculateToWall(StopActiveRegions[edr.stop],
    TrStr, TrDir) * 1E-2;

  Int_t mu_eff_x_bin =
      MuonKinematics_seleff->GetXaxis()->FindFixBin(
        edr.GetProjection(DepositsSummary::kEFSLep_True));
  Int_t mu_eff_y_bin =
      MuonKinematics_seleff->GetYaxis()->FindFixBin(MuToWall);

  if (MuonKinematics_seleff->GetBinContent(mu_eff_x_bin, mu_eff_y_bin)) {
    return 1.0 /
        MuonKinematics_seleff->GetBinContent(mu_eff_x_bin, mu_eff_y_bin);
  } else {
    return 1;
  }

}

double EffCorrector::GetHadronKinematicsEffWeight(DepositsSummary const & edr){
  Int_t hadr_eff_x_bin = 0;
  Int_t hadr_eff_y_bin = 0;
  switch (mode) {
    case kEHadrAbsPos:{
      hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
          edr.GetProjection(DepositsSummary::kEHadr_True));
      hadr_eff_y_bin =
          ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtx[0]);
      break;
    }
    case kEHadrDetPos:{
      hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
          edr.GetProjection(DepositsSummary::kEHadr_True));
      hadr_eff_y_bin =
          ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);
      break;
      }
    case kEHadrVisDetPos:{
      hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
          edr.GetProjection(DepositsSummary::kEHadr_vis));
      hadr_eff_y_bin =
          ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);
      break;
    }
  }

  if (ShowerKinematics_seleff->GetBinContent(hadr_eff_x_bin,
    hadr_eff_y_bin)) {
    return 1.0 /
      ShowerKinematics_seleff->GetBinContent(hadr_eff_x_bin, hadr_eff_y_bin);
  } else {
    return 1;
  }
}
