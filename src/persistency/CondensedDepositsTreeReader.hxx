#ifndef CONDENSEDTREEREADER_HXX_SEEN
#define CONDENSEDTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"
#include "TObjString.h"

struct CondensedDeposits : public ITreeReader {

  CondensedDeposits();

  CondensedDeposits(std::string const &inputFiles,
    Int_t NXBins = 402, Int_t NMaxTrackSteps = 1000,
                    Double_t timesep_us = 0xdeadbeef);

  size_t GetNPassthroughParts();

  std::pair<Int_t, Double_t *> GetPassthroughPart(size_t i);

  bool AddPassthroughPart(Int_t PDG, Double_t *fourmom);

  Int_t NXBins;
  Int_t NMaxTrackSteps;
  Double_t timesep_us;

  static const Int_t kNYBins = 3;
  static const Int_t kNZBins = 3;
  static const Int_t kNMaxPassthroughParts = 100;

  Int_t EventNum;

  Double_t nu_4mom[4];
  Int_t nu_PDG;
  Double_t VertexPosition[3];
  TObjString *EventCode;

  Int_t PrimaryLepPDG;
  Double_t PrimaryLep_4mom[4];

  Double_t Q2_True;
  Double_t FourMomTransfer_True[4];
  Double_t y_True;
  Double_t W_Rest;

  Int_t NFSParts;
  Int_t FSPart_PDG[kNMaxPassthroughParts];
  Int_t NFSPart4MomEntries;
  Double_t FSPart_4Mom[kNMaxPassthroughParts * 4];

  Int_t NLep;
  Int_t NPi0;
  Int_t NPiC;
  Int_t NProton;
  Int_t NNeutron;
  Int_t NGamma;
  Int_t NOther;
  Int_t NBaryonicRes;
  Int_t NAntiNucleons;

  Double_t EKinPi0_True;
  Double_t EMassPi0_True;
  Double_t EKinPiC_True;
  Double_t EMassPiC_True;
  Double_t EKinProton_True;
  Double_t EMassProton_True;
  Double_t EKinNeutron_True;
  Double_t EMassNeutron_True;
  Double_t EGamma_True;
  Double_t EOther_True;
  Double_t ENonPrimaryLep_True;
  Double_t KENuclearRemnant_True;

  Double_t TotalFS_3mom[3];

  Double_t *LepDep_1D;
  Double_t *HadDep_1D;
  Double_t *ProtonDep_1D;
  Double_t *NeutronDep_1D;
  Double_t *NeutronDep_ChrgWSumTime_1D;
  Double_t *PiCDep_1D;
  Double_t *Pi0Dep_1D;
  Double_t *OtherDep_1D;
  Double_t *NuclRemDep_1D;

  Double_t *LepDep_timesep_1D;
  Double_t *HadDep_timesep_1D;
  Double_t *ProtonDep_timesep_1D;
  Double_t *NeutronDep_timesep_1D;
  Double_t *PiCDep_timesep_1D;
  Double_t *Pi0Dep_timesep_1D;
  Double_t *OtherDep_timesep_1D;

  Double_t *LepDaughterDep_1D;
  Double_t *HadDaughterDep_1D;
  Double_t *ProtonDaughterDep_1D;
  Double_t *NeutronDaughterDep_1D;
  Double_t *NeutronDaughterDep_ChrgWSumTime_1D;
  Double_t *PiCDaughterDep_1D;
  Double_t *Pi0DaughterDep_1D;
  Double_t *OtherDaughterDep_1D;

  Double_t *LepDaughterDep_timesep_1D;
  Double_t *HadDaughterDep_timesep_1D;
  Double_t *ProtonDaughterDep_timesep_1D;
  Double_t *NeutronDaughterDep_timesep_1D;
  Double_t *PiCDaughterDep_timesep_1D;
  Double_t *Pi0DaughterDep_timesep_1D;
  Double_t *OtherDaughterDep_timesep_1D;

  Double_t **LepDep_2D;
  Double_t **HadDep_2D;
  Double_t **ProtonDep_2D;
  Double_t **NeutronDep_2D;
  Double_t **NeutronDep_ChrgWSumTime_2D;
  Double_t **PiCDep_2D;
  Double_t **Pi0Dep_2D;
  Double_t **OtherDep_2D;
  Double_t **NuclRemDep_2D;

  Double_t **LepDep_timesep_2D;
  Double_t **HadDep_timesep_2D;
  Double_t **ProtonDep_timesep_2D;
  Double_t **NeutronDep_timesep_2D;
  Double_t **PiCDep_timesep_2D;
  Double_t **Pi0Dep_timesep_2D;
  Double_t **OtherDep_timesep_2D;

  Double_t **LepDaughterDep_2D;
  Double_t **HadDaughterDep_2D;
  Double_t **ProtonDaughterDep_2D;
  Double_t **NeutronDaughterDep_2D;
  Double_t **NeutronDaughterDep_ChrgWSumTime_2D;
  Double_t **PiCDaughterDep_2D;
  Double_t **Pi0DaughterDep_2D;
  Double_t **OtherDaughterDep_2D;

  Double_t **LepDaughterDep_timesep_2D;
  Double_t **HadDaughterDep_timesep_2D;
  Double_t **ProtonDaughterDep_timesep_2D;
  Double_t **NeutronDaughterDep_timesep_2D;
  Double_t **PiCDaughterDep_timesep_2D;
  Double_t **Pi0DaughterDep_timesep_2D;
  Double_t **OtherDaughterDep_timesep_2D;

  Double_t ***LepDep;
  Double_t ***HadDep;
  Double_t ***ProtonDep;
  Double_t ***NeutronDep;
  Double_t ***NeutronDep_ChrgWSumTime;
  Double_t ***PiCDep;
  Double_t ***Pi0Dep;
  Double_t ***OtherDep;
  Double_t ***NuclRemDep;

  Double_t ***LepDep_timesep;
  Double_t ***HadDep_timesep;
  Double_t ***ProtonDep_timesep;
  Double_t ***NeutronDep_timesep;
  Double_t ***PiCDep_timesep;
  Double_t ***Pi0Dep_timesep;
  Double_t ***OtherDep_timesep;

  Double_t ***LepDaughterDep;
  Double_t ***HadDaughterDep;
  Double_t ***ProtonDaughterDep;
  Double_t ***NeutronDaughterDep;
  Double_t ***NeutronDaughterDep_ChrgWSumTime;
  Double_t ***PiCDaughterDep;
  Double_t ***Pi0DaughterDep;
  Double_t ***OtherDaughterDep;

  Double_t ***LepDaughterDep_timesep;
  Double_t ***HadDaughterDep_timesep;
  Double_t ***ProtonDaughterDep_timesep;
  Double_t ***NeutronDaughterDep_timesep;
  Double_t ***PiCDaughterDep_timesep;
  Double_t ***Pi0DaughterDep_timesep;
  Double_t ***OtherDaughterDep_timesep;

  Double_t *MuonTrackPos_1D;
  Double_t *MuonTrackMom_1D;

  Int_t NMuonTrackSteps;
  Double_t **MuonTrackPos;
  Double_t **MuonTrackMom;

  std::string TreeName();
  void Reset();
  void SetBranchAddresses();

  void AllocateArrays();
  void DeAllocateArrays();
  static CondensedDeposits *MakeTreeWriter(Int_t NXBins,
                                           Int_t NMaxTrackSteps,
                                           Double_t timesep_us = 0xdeadbeef);
  ~CondensedDeposits();
};

#endif
