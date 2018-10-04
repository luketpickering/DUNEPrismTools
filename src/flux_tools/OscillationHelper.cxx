#include "OscillationHelper.hxx"

#include "ROOTUtility.hxx"

#include "OscillationParametersTreeReader.hxx"

#include <iostream>
#include <string>
#include <algorithm>

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
  std::copy_n(op.OscParams,6,OscParams);


  std::cout << "[INFO]: Using oscillation parameters: " << std::endl
            << "\tSin^2(Theta_12) = " << OscParams[0] << std::endl;
  std::cout << "\tSin^2(Theta_13) = " << OscParams[1] << std::endl;
  std::cout << "\tSin^2(Theta_23) = " << OscParams[2] << std::endl;

  std::cout << "\tDm^2_21 = " << OscParams[3] << " (GeV^2) " << std::endl;
  std::cout << "\t|Dm^2_Atm| = " << OscParams[4] << " (GeV^2) " << std::endl;

  std::cout << "\tdcp = " << OscParams[5] << std::endl;
  std::cout << "\tBeam dip angle = " << DipAngle_degrees << " degrees" << std::endl;

  LengthParam = cos((90.0 + DipAngle_degrees) * deg2rad);
  IsSetUp = true;
}
void OscillationHelper::Setup(double OscParams[6], double DipAngle_degrees){
  std::copy_n(OscParams,6,this->OscParams);
  this->DipAngle_degrees = DipAngle_degrees;

  std::cout << "[INFO]: Using oscillation parameters: " << std::endl
            << "\tSin^2(Theta_12) = " << this->OscParams[0] << std::endl;
  std::cout << "\tSin^2(Theta_13) = " << this->OscParams[1] << std::endl;
  std::cout << "\tSin^2(Theta_23) = " << this->OscParams[2] << std::endl;

  std::cout << "\tDm^2_21 = " << this->OscParams[3] << " (GeV^2) " << std::endl;
  std::cout << "\t|Dm^2_Atm| = " << this->OscParams[4] << " (GeV^2) " << std::endl;

  std::cout << "\tdcp = " << this->OscParams[5] << std::endl;
  std::cout << "\tBeam dip angle = " << DipAngle_degrees << " degrees" << std::endl;

  LengthParam = cos((90.0 + DipAngle_degrees) * deg2rad);
  IsSetUp = true;
}

void OscillationHelper::SetOscillationChannel(int PDGFrom, int PDGTo) {
  FromPDG = PDGFrom;
  ToPDG = PDGTo;
  FromType = GetNuType(PDGFrom);
  ToType = GetNuType(PDGTo);
}

double OscillationHelper::GetWeight(double ENu_GeV) {
  if(!IsSetUp){
    std::cout << "[ERROR]: Attempted to calculate oscillation weights without "
    "setting up the OscillationHelper instance." << std::endl;
    throw;
  }
  bp.SetMNS(OscParams[0], OscParams[1], OscParams[2], OscParams[3],
            OscParams[4], OscParams[5], ENu_GeV, true, FromType);
  bp.DefinePath(LengthParam, 0);
  bp.propagate(ToType);

  return bp.GetProb(FromType, ToType);
}

void OscillationHelper::WriteConfigTree(TDirectory *dir){
  dir->cd();
  OscillationParameters * tw = OscillationParameters::MakeTreeWriter();

  tw->DipAngle_degrees = DipAngle_degrees;
  std::copy_n(OscParams,6,tw->OscParams);
  tw->FromNuPDG = FromPDG;
  tw->ToNuPDG = ToPDG;

  tw->Fill();

}
