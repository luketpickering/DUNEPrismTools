#include "OscillationHelper.hxx"

#include "ROOTUtility.hxx"

#include "OscillationParametersTreeReader.hxx"

#include <algorithm>
#include <iostream>
#include <string>

static const double deg2rad = asin(1) / 90.0;

OscillationHelper::NuTypes OscillationHelper::GetNuType(int pdg) {
  switch (pdg) {
  case 16:
    return kNutauType;
  case 14:
    return kNumuType;
  case 12:
    return kNueType;
  case -16:
    return kNutaubarType;
  case -14:
    return kNumubarType;
  case -12:
    return kNuebarType;
  default: {
    std::cout << "[ERROR]: Attempting to convert \"neutrino pdg\": " << pdg
              << std::endl;
    exit(1);
  }
  }
}

void OscillationHelper::Setup(std::string const &FileWithConfTree) {
  SetOscillationChannel(14, 14);

  OscillationParameters op(FileWithConfTree);
  DipAngle_degrees = op.DipAngle_degrees;
  std::copy_n(op.OscParams, 6, OscParams);

  std::cout << "[INFO]: Using oscillation parameters: " << std::endl
            << "\tSin^2(Theta_12) = " << this->OscParams[0] << std::endl;
  std::cout << "\tSin^2(Theta_13) = " << this->OscParams[1] << std::endl;
  std::cout << "\tSin^2(Theta_23) = " << this->OscParams[2] << std::endl;

  std::cout << "\tDm^2_21 = " << this->OscParams[3] << " (GeV^2) " << std::endl;
  std::cout << "\t|Dm^2_Atm| = " << this->OscParams[4] << " (GeV^2) "
            << std::endl;

  std::cout << "\tdcp = " << this->OscParams[5] << std::endl;
  std::cout << "\tBeam dip angle = " << DipAngle_degrees << " degrees"
            << std::endl;

  LengthParam = cos((90.0 + DipAngle_degrees) * deg2rad);
  IsSetUp = true;
}
void OscillationHelper::Setup(double OscParams[6], double DipAngle_degrees) {
  std::copy_n(OscParams, 6, this->OscParams);
  this->DipAngle_degrees = DipAngle_degrees;

  std::cout << "[INFO]: Using oscillation parameters: " << std::endl
            << "\tSin^2(Theta_12) = " << this->OscParams[0] << std::endl;
  std::cout << "\tSin^2(Theta_13) = " << this->OscParams[1] << std::endl;
  std::cout << "\tSin^2(Theta_23) = " << this->OscParams[2] << std::endl;

  std::cout << "\tDm^2_21 = " << this->OscParams[3] << " (GeV^2) " << std::endl;
  std::cout << "\t|Dm^2_Atm| = " << this->OscParams[4] << " (GeV^2) "
            << std::endl;

  std::cout << "\tdcp = " << this->OscParams[5] << std::endl;
  std::cout << "\tBeam dip angle = " << DipAngle_degrees << " degrees"
            << std::endl;

  LengthParam = cos((90.0 + DipAngle_degrees) * deg2rad);
  IsSetUp = true;
}

#ifdef USE_FHICL
void OscillationHelper::Setup(fhicl::ParameterSet const &ps) {
  fhicl::ParameterSet const &op_ps =
      ps.get<fhicl::ParameterSet>("OscillationParameters");
  fhicl::ParameterSet const &bl_ps = ps.get<fhicl::ParameterSet>("Baseline");

  OscParams[0] = op_ps.get<double>("S2Th12");
  OscParams[1] = op_ps.get<double>("S2Th13");
  OscParams[2] = op_ps.get<double>("S2Th23");
  OscParams[3] = op_ps.get<double>("Dm2_21");
  OscParams[4] = op_ps.get<double>("Dm2_Atm");
  OscParams[5] = op_ps.get<double>("dcp");

  DipAngle_degrees = bl_ps.get<double>("dip_angle");

  std::cout << "[INFO]: Using oscillation parameters: " << std::endl
            << "\tSin^2(Theta_12) = " << this->OscParams[0] << std::endl;
  std::cout << "\tSin^2(Theta_13) = " << this->OscParams[1] << std::endl;
  std::cout << "\tSin^2(Theta_23) = " << this->OscParams[2] << std::endl;

  std::cout << "\tDm^2_21 = " << this->OscParams[3] << " (GeV^2) " << std::endl;
  std::cout << "\t|Dm^2_Atm| = " << this->OscParams[4] << " (GeV^2) "
            << std::endl;

  std::cout << "\tdcp = " << this->OscParams[5] << std::endl;
  std::cout << "\tBeam dip angle = " << DipAngle_degrees << " degrees"
            << std::endl;

  LengthParam = cos((90.0 + DipAngle_degrees) * deg2rad);
  IsSetUp = true;
}
#endif

void OscillationHelper::SetOscillationChannel(int PDGFrom, int PDGTo) {
  FromPDG = PDGFrom;
  ToPDG = PDGTo;
  FromType = GetNuType(PDGFrom);
  ToType = GetNuType(PDGTo);
}

double OscillationHelper::GetWeight(double ENu_GeV) {
  if (!IsSetUp) {
    std::cout << "[ERROR]: Attempted to calculate oscillation weights without "
                 "setting up the OscillationHelper instance."
              << std::endl;
    throw;
  }
  bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
            OscParams[4], OscParams[5], ENu_GeV, true, FromType);
  bp.DefinePath(LengthParam, 0);
  bp.propagate(ToType);

  return bp.GetProb(FromType, ToType);
}
void OscillationHelper::OscillateHistogram(std::unique_ptr<TH1D> &h) {
  for (Int_t bi_it = 1; bi_it < h->GetXaxis()->GetNbins() + 1; ++bi_it) {
    double weight = GetWeight(h->GetXaxis()->GetBinCenter(bi_it));
    h->SetBinContent(bi_it, weight * h->GetBinContent(bi_it));
    h->SetBinError(bi_it, weight * h->GetBinError(bi_it));
  }
}

void OscillationHelper::WriteConfigTree(TDirectory *dir) {
  dir->cd();
  OscillationParameters *tw = OscillationParameters::MakeTreeWriter();

  tw->DipAngle_degrees = DipAngle_degrees;
  std::copy_n(OscParams, 6, tw->OscParams);
  tw->FromNuPDG = FromPDG;
  tw->ToNuPDG = ToPDG;

  tw->Fill();
}
