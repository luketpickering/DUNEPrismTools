#include "CondensedDepositsTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <iostream>
#include <algorithm>

  CondensedDeposits::CondensedDeposits()
      : NXBins(402),
        NMaxTrackSteps(1000),
        timesep_us(0xdeadbeef),
        EventCode(nullptr) {}

  CondensedDeposits::CondensedDeposits(std::string const &inputFiles,
                    Int_t NXBins, Int_t NMaxTrackSteps,
                    Double_t timesep_us) : CondensedDeposits() {

    this->NXBins = NXBins;
    this->NMaxTrackSteps = NMaxTrackSteps;
    this->timesep_us = timesep_us;
    EventCode = nullptr;

    LoadTree(inputFiles);
  }

  size_t CondensedDeposits::GetNPassthroughParts() { return NFSParts; }

  std::pair<Int_t, Double_t *> CondensedDeposits::GetPassthroughPart(size_t i) {
    if (i >= GetNPassthroughParts()) {
      std::cout << "[ERROR]: Asked for passthrough part " << i
                << ", but we only have " << GetNPassthroughParts()
                << " in this event." << std::endl;
      throw;
    }

    return std::pair<Int_t, Double_t *>(FSPart_PDG[i], &FSPart_4Mom[i * 4]);
  }

  bool CondensedDeposits::AddPassthroughPart(Int_t PDG, Double_t *fourmom) {
    if (NFSParts >= kNMaxPassthroughParts) {
      std::cout << "[WARN]: Ignoring passthrough part as it would exceed the "
                   "maximum of "
                << kNMaxPassthroughParts << " particles in this event."
                << std::endl;
      return false;
    }

    FSPart_PDG[NFSParts] = PDG;
    std::copy_n(fourmom, 4, &FSPart_4Mom[NFSParts * 4]);
    NFSParts++;
    NFSPart4MomEntries += 4;
    return true;
  }

  std::string CondensedDeposits::TreeName(){
    return "CondensedDepositsTree";
  }

  void CondensedDeposits::Reset() {
    EventNum = -1;
    std::fill_n(nu_4mom, 4, -1);
    nu_PDG = 0;
    std::fill_n(VertexPosition, 3, 0);
    PrimaryLepPDG = 0;
    std::fill_n(PrimaryLep_4mom, 4, 0);
    Q2_True = -1;
    std::fill_n(FourMomTransfer_True, 4, 0);
    y_True = 0;
    W_Rest = 0;

    NFSParts = 0;
    std::fill_n(FSPart_PDG, kNMaxPassthroughParts, 0);
    NFSPart4MomEntries = 0;
    std::fill_n(FSPart_4Mom, kNMaxPassthroughParts * 4, 0);

    NLep = 0;
    NPi0 = 0;
    NPiC = 0;
    NProton = 0;
    NNeutron = 0;
    NGamma = 0;
    NOther = 0;
    NBaryonicRes = 0;
    NAntiNucleons = 0;
    EKinPi0_True = 0;
    EMassPi0_True = 0;
    EKinPiC_True = 0;
    EMassPiC_True = 0;
    EKinProton_True = 0;
    EMassProton_True = 0;
    EKinNeutron_True = 0;
    EMassNeutron_True = 0;
    EGamma_True = 0;
    EOther_True = 0;
    ENonPrimaryLep_True = 0;

    KENuclearRemnant_True = 0;
    std::fill_n(TotalFS_3mom, 3, 0);

    std::fill_n(LepDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(HadDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(ProtonDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(NeutronDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(NeutronDep_ChrgWSumTime_1D, NXBins * 3 * 3, 0);
    std::fill_n(PiCDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(Pi0Dep_1D, NXBins * 3 * 3, 0);
    std::fill_n(OtherDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(NuclRemDep_1D, NXBins * 3 * 3, 0);

    std::fill_n(LepDaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(HadDaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(ProtonDaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(NeutronDaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(NeutronDaughterDep_ChrgWSumTime_1D, NXBins * 3 * 3, 0);
    std::fill_n(PiCDaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(Pi0DaughterDep_1D, NXBins * 3 * 3, 0);
    std::fill_n(OtherDaughterDep_1D, NXBins * 3 * 3, 0);

    if (timesep_us != 0xdeadbeef) {
      std::fill_n(LepDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(HadDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(ProtonDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(NeutronDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(PiCDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(Pi0Dep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(OtherDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(LepDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(HadDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(ProtonDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(NeutronDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(PiCDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(Pi0DaughterDep_timesep_1D, NXBins * 3 * 3, 0);
      std::fill_n(OtherDaughterDep_timesep_1D, NXBins * 3 * 3, 0);
    }

    NMuonTrackSteps = 0;
    std::fill_n(MuonTrackPos_1D, NMaxTrackSteps * 3, 0xdeadbeef);
    std::fill_n(MuonTrackMom_1D, NMaxTrackSteps * 3, 0xdeadbeef);
  }

  void CondensedDeposits::SetBranchAddresses() {
    AllocateArrays();

    tree->SetBranchAddress("EventNum", &EventNum);
    tree->SetBranchAddress("nu_4mom", &nu_4mom);
    tree->SetBranchAddress("nu_PDG", &nu_PDG);
    tree->SetBranchAddress("VertexPosition", &VertexPosition);
    tree->SetBranchAddress("EventCode", &EventCode);
    tree->SetBranchAddress("PrimaryLepPDG", &PrimaryLepPDG);
    tree->SetBranchAddress("PrimaryLep_4mom", &PrimaryLep_4mom);
    tree->SetBranchAddress("Q2_True", &Q2_True);
    tree->SetBranchAddress("FourMomTransfer_True", &FourMomTransfer_True);
    tree->SetBranchAddress("y_True", &y_True);
    tree->SetBranchAddress("W_Rest", &W_Rest);
    tree->SetBranchAddress("NFSParts", &NFSParts);
    tree->SetBranchAddress("FSPart_PDG", &FSPart_PDG);
    tree->SetBranchAddress("NFSPart4MomEntries", &NFSPart4MomEntries);
    tree->SetBranchAddress("FSPart_4Mom", &FSPart_4Mom);
    tree->SetBranchAddress("NLep", &NLep);
    tree->SetBranchAddress("NPi0", &NPi0);
    tree->SetBranchAddress("NPiC", &NPiC);
    tree->SetBranchAddress("NProton", &NProton);
    tree->SetBranchAddress("NNeutron", &NNeutron);
    tree->SetBranchAddress("NGamma", &NGamma);
    tree->SetBranchAddress("NOther", &NOther);
    tree->SetBranchAddress("NBaryonicRes", &NBaryonicRes);
    tree->SetBranchAddress("NAntiNucleons", &NAntiNucleons);

    tree->SetBranchAddress("EKinPi0_True", &EKinPi0_True);
    tree->SetBranchAddress("EMassPi0_True", &EMassPi0_True);
    tree->SetBranchAddress("EKinPiC_True", &EKinPiC_True);
    tree->SetBranchAddress("EMassPiC_True", &EMassPiC_True);
    tree->SetBranchAddress("EKinProton_True", &EKinProton_True);
    tree->SetBranchAddress("EMassProton_True", &EMassProton_True);
    tree->SetBranchAddress("EKinNeutron_True", &EKinNeutron_True);
    tree->SetBranchAddress("EMassNeutron_True", &EMassNeutron_True);
    tree->SetBranchAddress("EGamma_True", &EGamma_True);
    tree->SetBranchAddress("EOther_True", &EOther_True);
    tree->SetBranchAddress("ENonPrimaryLep_True", &ENonPrimaryLep_True);
    tree->SetBranchAddress("KENuclearRemnant_True", &KENuclearRemnant_True);
    tree->SetBranchAddress("TotalFS_3mom", &TotalFS_3mom);
    tree->SetBranchAddress("LepDep", LepDep_1D);
    tree->SetBranchAddress("HadDep", HadDep_1D);
    tree->SetBranchAddress("ProtonDep", ProtonDep_1D);
    tree->SetBranchAddress("NeutronDep", NeutronDep_1D);
    tree->SetBranchAddress("NeutronDep_ChrgWSumTime",
                           NeutronDep_ChrgWSumTime_1D);
    tree->SetBranchAddress("PiCDep", PiCDep_1D);
    tree->SetBranchAddress("Pi0Dep", Pi0Dep_1D);
    tree->SetBranchAddress("OtherDep", OtherDep_1D);
    tree->SetBranchAddress("NuclRemDep", NuclRemDep_1D);
    tree->SetBranchAddress("LepDaughterDep", LepDaughterDep_1D);
    tree->SetBranchAddress("HadDaughterDep", HadDaughterDep_1D);
    tree->SetBranchAddress("ProtonDaughterDep", ProtonDaughterDep_1D);
    tree->SetBranchAddress("NeutronDaughterDep", NeutronDaughterDep_1D);
    tree->SetBranchAddress("NeutronDaughterDep_ChrgWSumTime",
                           NeutronDaughterDep_ChrgWSumTime_1D);
    tree->SetBranchAddress("PiCDaughterDep", PiCDaughterDep_1D);
    tree->SetBranchAddress("Pi0DaughterDep", Pi0DaughterDep_1D);
    tree->SetBranchAddress("OtherDaughterDep", OtherDaughterDep_1D);

    if (timesep_us != 0xdeadbeef) {
      tree->SetBranchAddress("LepDep_timesep", LepDep_timesep_1D);
      tree->SetBranchAddress("HadDep_timesep", HadDep_timesep_1D);
      tree->SetBranchAddress("ProtonDep_timesep", ProtonDep_timesep_1D);
      tree->SetBranchAddress("NeutronDep_timesep", NeutronDep_timesep_1D);
      tree->SetBranchAddress("PiCDep_timesep", PiCDep_timesep_1D);
      tree->SetBranchAddress("Pi0Dep_timesep", Pi0Dep_timesep_1D);
      tree->SetBranchAddress("OtherDep_timesep", OtherDep_timesep_1D);
      tree->SetBranchAddress("LepDaughterDep_timesep",
                             LepDaughterDep_timesep_1D);
      tree->SetBranchAddress("HadDaughterDep_timesep",
                             HadDaughterDep_timesep_1D);
      tree->SetBranchAddress("ProtonDaughterDep_timesep",
                             ProtonDaughterDep_timesep_1D);
      tree->SetBranchAddress("NeutronDaughterDep_timesep",
                             NeutronDaughterDep_timesep_1D);
      tree->SetBranchAddress("PiCDaughterDep_timesep",
                             PiCDaughterDep_timesep_1D);
      tree->SetBranchAddress("Pi0DaughterDep_timesep",
                             Pi0DaughterDep_timesep_1D);
      tree->SetBranchAddress("OtherDaughterDep_timesep",
                             OtherDaughterDep_timesep_1D);
    }

    tree->SetBranchAddress("NMuonTrackSteps", &NMuonTrackSteps);
    tree->SetBranchAddress("MuonTrackPos", MuonTrackPos_1D);
    tree->SetBranchAddress("MuonTrackMom", MuonTrackMom_1D);
  }

  void CondensedDeposits::AllocateArrays() {
    LepDep_1D = new Double_t[NXBins * 3 * 3];
    HadDep_1D = new Double_t[NXBins * 3 * 3];
    ProtonDep_1D = new Double_t[NXBins * 3 * 3];
    NeutronDep_1D = new Double_t[NXBins * 3 * 3];
    NeutronDep_ChrgWSumTime_1D = new Double_t[NXBins * 3 * 3];
    PiCDep_1D = new Double_t[NXBins * 3 * 3];
    Pi0Dep_1D = new Double_t[NXBins * 3 * 3];
    OtherDep_1D = new Double_t[NXBins * 3 * 3];
    NuclRemDep_1D = new Double_t[NXBins * 3 * 3];

    LepDep_2D = new Double_t *[NXBins * 3];
    HadDep_2D = new Double_t *[NXBins * 3];
    ProtonDep_2D = new Double_t *[NXBins * 3];
    NeutronDep_2D = new Double_t *[NXBins * 3];
    NeutronDep_ChrgWSumTime_2D = new Double_t *[NXBins * 3];
    PiCDep_2D = new Double_t *[NXBins * 3];
    Pi0Dep_2D = new Double_t *[NXBins * 3];
    OtherDep_2D = new Double_t *[NXBins * 3];
    NuclRemDep_2D = new Double_t *[NXBins * 3];

    LepDep = new Double_t **[NXBins];
    HadDep = new Double_t **[NXBins];
    ProtonDep = new Double_t **[NXBins];
    NeutronDep = new Double_t **[NXBins];
    NeutronDep_ChrgWSumTime = new Double_t **[NXBins];
    PiCDep = new Double_t **[NXBins];
    Pi0Dep = new Double_t **[NXBins];
    OtherDep = new Double_t **[NXBins];
    NuclRemDep = new Double_t **[NXBins];

    LepDaughterDep_1D = new Double_t[NXBins * 3 * 3];
    HadDaughterDep_1D = new Double_t[NXBins * 3 * 3];
    ProtonDaughterDep_1D = new Double_t[NXBins * 3 * 3];
    NeutronDaughterDep_1D = new Double_t[NXBins * 3 * 3];
    NeutronDaughterDep_ChrgWSumTime_1D = new Double_t[NXBins * 3 * 3];
    PiCDaughterDep_1D = new Double_t[NXBins * 3 * 3];
    Pi0DaughterDep_1D = new Double_t[NXBins * 3 * 3];
    OtherDaughterDep_1D = new Double_t[NXBins * 3 * 3];

    LepDaughterDep_2D = new Double_t *[NXBins * 3];
    HadDaughterDep_2D = new Double_t *[NXBins * 3];
    ProtonDaughterDep_2D = new Double_t *[NXBins * 3];
    NeutronDaughterDep_2D = new Double_t *[NXBins * 3];
    NeutronDaughterDep_ChrgWSumTime_2D = new Double_t *[NXBins * 3];
    PiCDaughterDep_2D = new Double_t *[NXBins * 3];
    Pi0DaughterDep_2D = new Double_t *[NXBins * 3];
    OtherDaughterDep_2D = new Double_t *[NXBins * 3];

    LepDaughterDep = new Double_t **[NXBins];
    HadDaughterDep = new Double_t **[NXBins];
    ProtonDaughterDep = new Double_t **[NXBins];
    NeutronDaughterDep = new Double_t **[NXBins];
    NeutronDaughterDep_ChrgWSumTime = new Double_t **[NXBins];
    PiCDaughterDep = new Double_t **[NXBins];
    Pi0DaughterDep = new Double_t **[NXBins];
    OtherDaughterDep = new Double_t **[NXBins];

    for (Int_t i = 0; i < NXBins; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        LepDep_2D[i * 3 + j] = &LepDep_1D[i * 3 * 3 + j * 3];
        HadDep_2D[i * 3 + j] = &HadDep_1D[i * 3 * 3 + j * 3];
        ProtonDep_2D[i * 3 + j] = &ProtonDep_1D[i * 3 * 3 + j * 3];
        NeutronDep_2D[i * 3 + j] = &NeutronDep_1D[i * 3 * 3 + j * 3];
        NeutronDep_ChrgWSumTime_2D[i * 3 + j] =
            &NeutronDep_ChrgWSumTime_1D[i * 3 * 3 + j * 3];
        PiCDep_2D[i * 3 + j] = &PiCDep_1D[i * 3 * 3 + j * 3];
        Pi0Dep_2D[i * 3 + j] = &Pi0Dep_1D[i * 3 * 3 + j * 3];
        OtherDep_2D[i * 3 + j] = &OtherDep_1D[i * 3 * 3 + j * 3];
        NuclRemDep_2D[i * 3 + j] = &NuclRemDep_1D[i * 3 * 3 + j * 3];

        LepDaughterDep_2D[i * 3 + j] = &LepDaughterDep_1D[i * 3 * 3 + j * 3];
        HadDaughterDep_2D[i * 3 + j] = &HadDaughterDep_1D[i * 3 * 3 + j * 3];
        ProtonDaughterDep_2D[i * 3 + j] =
            &ProtonDaughterDep_1D[i * 3 * 3 + j * 3];
        NeutronDaughterDep_2D[i * 3 + j] =
            &NeutronDaughterDep_1D[i * 3 * 3 + j * 3];
        NeutronDaughterDep_ChrgWSumTime_2D[i * 3 + j] =
            &NeutronDaughterDep_ChrgWSumTime_1D[i * 3 * 3 + j * 3];
        PiCDaughterDep_2D[i * 3 + j] = &PiCDaughterDep_1D[i * 3 * 3 + j * 3];
        Pi0DaughterDep_2D[i * 3 + j] = &Pi0DaughterDep_1D[i * 3 * 3 + j * 3];
        OtherDaughterDep_2D[i * 3 + j] =
            &OtherDaughterDep_1D[i * 3 * 3 + j * 3];
      }

      // LepDep[NXBins-1][2][2] = (&LepDep_2D[(NXBins-1)*3])[2][2] = ()
      LepDep[i] = &LepDep_2D[i * 3];
      HadDep[i] = &HadDep_2D[i * 3];
      ProtonDep[i] = &ProtonDep_2D[i * 3];
      NeutronDep[i] = &NeutronDep_2D[i * 3];
      NeutronDep_ChrgWSumTime[i] = &NeutronDep_ChrgWSumTime_2D[i * 3];
      PiCDep[i] = &PiCDep_2D[i * 3];
      Pi0Dep[i] = &Pi0Dep_2D[i * 3];
      OtherDep[i] = &OtherDep_2D[i * 3];
      NuclRemDep[i] = &NuclRemDep_2D[i * 3];

      LepDaughterDep[i] = &LepDaughterDep_2D[i * 3];
      HadDaughterDep[i] = &HadDaughterDep_2D[i * 3];
      ProtonDaughterDep[i] = &ProtonDaughterDep_2D[i * 3];
      NeutronDaughterDep[i] = &NeutronDaughterDep_2D[i * 3];
      NeutronDaughterDep_ChrgWSumTime[i] =
          &NeutronDaughterDep_ChrgWSumTime_2D[i * 3];
      PiCDaughterDep[i] = &PiCDaughterDep_2D[i * 3];
      Pi0DaughterDep[i] = &Pi0DaughterDep_2D[i * 3];
      OtherDaughterDep[i] = &OtherDaughterDep_2D[i * 3];
    }

    if (timesep_us != 0xdeadbeef) {
      LepDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      HadDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      ProtonDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      NeutronDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      PiCDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      Pi0Dep_timesep_1D = new Double_t[NXBins * 3 * 3];
      OtherDep_timesep_1D = new Double_t[NXBins * 3 * 3];

      LepDep_timesep_2D = new Double_t *[NXBins * 3];
      HadDep_timesep_2D = new Double_t *[NXBins * 3];
      ProtonDep_timesep_2D = new Double_t *[NXBins * 3];
      NeutronDep_timesep_2D = new Double_t *[NXBins * 3];
      PiCDep_timesep_2D = new Double_t *[NXBins * 3];
      Pi0Dep_timesep_2D = new Double_t *[NXBins * 3];
      OtherDep_timesep_2D = new Double_t *[NXBins * 3];

      LepDep_timesep = new Double_t **[NXBins];
      HadDep_timesep = new Double_t **[NXBins];
      ProtonDep_timesep = new Double_t **[NXBins];
      NeutronDep_timesep = new Double_t **[NXBins];
      PiCDep_timesep = new Double_t **[NXBins];
      Pi0Dep_timesep = new Double_t **[NXBins];
      OtherDep_timesep = new Double_t **[NXBins];

      LepDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      HadDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      ProtonDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      NeutronDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      PiCDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      Pi0DaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];
      OtherDaughterDep_timesep_1D = new Double_t[NXBins * 3 * 3];

      LepDaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      HadDaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      ProtonDaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      NeutronDaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      PiCDaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      Pi0DaughterDep_timesep_2D = new Double_t *[NXBins * 3];
      OtherDaughterDep_timesep_2D = new Double_t *[NXBins * 3];

      LepDaughterDep_timesep = new Double_t **[NXBins];
      HadDaughterDep_timesep = new Double_t **[NXBins];
      ProtonDaughterDep_timesep = new Double_t **[NXBins];
      NeutronDaughterDep_timesep = new Double_t **[NXBins];
      PiCDaughterDep_timesep = new Double_t **[NXBins];
      Pi0DaughterDep_timesep = new Double_t **[NXBins];
      OtherDaughterDep_timesep = new Double_t **[NXBins];

      for (Int_t i = 0; i < NXBins; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          LepDep_timesep_2D[i * 3 + j] = &LepDep_timesep_1D[i * 3 * 3 + j * 3];
          HadDep_timesep_2D[i * 3 + j] = &HadDep_timesep_1D[i * 3 * 3 + j * 3];
          ProtonDep_timesep_2D[i * 3 + j] =
              &ProtonDep_timesep_1D[i * 3 * 3 + j * 3];
          NeutronDep_timesep_2D[i * 3 + j] =
              &NeutronDep_timesep_1D[i * 3 * 3 + j * 3];
          PiCDep_timesep_2D[i * 3 + j] = &PiCDep_timesep_1D[i * 3 * 3 + j * 3];
          Pi0Dep_timesep_2D[i * 3 + j] = &Pi0Dep_timesep_1D[i * 3 * 3 + j * 3];
          OtherDep_timesep_2D[i * 3 + j] =
              &OtherDep_timesep_1D[i * 3 * 3 + j * 3];

          LepDaughterDep_timesep_2D[i * 3 + j] =
              &LepDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          HadDaughterDep_timesep_2D[i * 3 + j] =
              &HadDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          ProtonDaughterDep_timesep_2D[i * 3 + j] =
              &ProtonDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          NeutronDaughterDep_timesep_2D[i * 3 + j] =
              &NeutronDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          PiCDaughterDep_timesep_2D[i * 3 + j] =
              &PiCDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          Pi0DaughterDep_timesep_2D[i * 3 + j] =
              &Pi0DaughterDep_timesep_1D[i * 3 * 3 + j * 3];
          OtherDaughterDep_timesep_2D[i * 3 + j] =
              &OtherDaughterDep_timesep_1D[i * 3 * 3 + j * 3];
        }

        // LepDep_timesep[NXBins-1][2][2] =
        // (&LepDep_timesep_2D[(NXBins-1)*3])[2][2] = ()
        LepDep_timesep[i] = &LepDep_timesep_2D[i * 3];
        HadDep_timesep[i] = &HadDep_timesep_2D[i * 3];
        ProtonDep_timesep[i] = &ProtonDep_timesep_2D[i * 3];
        NeutronDep_timesep[i] = &NeutronDep_timesep_2D[i * 3];
        PiCDep_timesep[i] = &PiCDep_timesep_2D[i * 3];
        Pi0Dep_timesep[i] = &Pi0Dep_timesep_2D[i * 3];
        OtherDep_timesep[i] = &OtherDep_timesep_2D[i * 3];

        LepDaughterDep_timesep[i] = &LepDaughterDep_timesep_2D[i * 3];
        HadDaughterDep_timesep[i] = &HadDaughterDep_timesep_2D[i * 3];
        ProtonDaughterDep_timesep[i] = &ProtonDaughterDep_timesep_2D[i * 3];
        NeutronDaughterDep_timesep[i] = &NeutronDaughterDep_timesep_2D[i * 3];
        PiCDaughterDep_timesep[i] = &PiCDaughterDep_timesep_2D[i * 3];
        Pi0DaughterDep_timesep[i] = &Pi0DaughterDep_timesep_2D[i * 3];
        OtherDaughterDep_timesep[i] = &OtherDaughterDep_timesep_2D[i * 3];
      }
    }

    MuonTrackPos_1D = new Double_t[NMaxTrackSteps * 3];
    MuonTrackMom_1D = new Double_t[NMaxTrackSteps * 3];

    MuonTrackPos = new Double_t *[NMaxTrackSteps];
    MuonTrackMom = new Double_t *[NMaxTrackSteps];

    for (Int_t st_it = 0; st_it < NMaxTrackSteps; ++st_it) {
      MuonTrackPos[st_it] = &MuonTrackPos_1D[st_it * 3];
      MuonTrackMom[st_it] = &MuonTrackMom_1D[st_it * 3];
    }

    Reset();
  }

  void CondensedDeposits::DeAllocateArrays() {
    delete[] LepDep_2D;
    delete[] LepDep_1D;
    delete[] HadDep_2D;
    delete[] HadDep_1D;
    delete[] ProtonDep_2D;
    delete[] ProtonDep_1D;
    delete[] NeutronDep_2D;
    delete[] NeutronDep_ChrgWSumTime_2D;
    delete[] NeutronDep_1D;
    delete[] NeutronDep_ChrgWSumTime_1D;
    delete[] PiCDep_2D;
    delete[] PiCDep_1D;
    delete[] Pi0Dep_2D;
    delete[] Pi0Dep_1D;
    delete[] OtherDep_2D;
    delete[] OtherDep_1D;
    delete[] NuclRemDep_2D;
    delete[] NuclRemDep_1D;
    delete[] LepDaughterDep_2D;
    delete[] LepDaughterDep_1D;
    delete[] HadDaughterDep_2D;
    delete[] HadDaughterDep_1D;
    delete[] ProtonDaughterDep_2D;
    delete[] ProtonDaughterDep_1D;
    delete[] NeutronDaughterDep_2D;
    delete[] NeutronDaughterDep_ChrgWSumTime_2D;
    delete[] NeutronDaughterDep_1D;
    delete[] NeutronDaughterDep_ChrgWSumTime_1D;
    delete[] PiCDaughterDep_2D;
    delete[] PiCDaughterDep_1D;
    delete[] Pi0DaughterDep_2D;
    delete[] Pi0DaughterDep_1D;
    delete[] OtherDaughterDep_2D;
    delete[] OtherDaughterDep_1D;
    delete[] LepDep;
    delete[] HadDep;
    delete[] ProtonDep;
    delete[] NeutronDep;
    delete[] NeutronDep_ChrgWSumTime;
    delete[] PiCDep;
    delete[] Pi0Dep;
    delete[] OtherDep;
    delete[] NuclRemDep;
    delete[] LepDaughterDep;
    delete[] HadDaughterDep;
    delete[] ProtonDaughterDep;
    delete[] NeutronDaughterDep;
    delete[] NeutronDaughterDep_ChrgWSumTime;
    delete[] PiCDaughterDep;
    delete[] Pi0DaughterDep;
    delete[] OtherDaughterDep;

    if (timesep_us != 0xdeadbeef) {
      delete[] LepDep_timesep_2D;
      delete[] LepDep_timesep_1D;
      delete[] HadDep_timesep_2D;
      delete[] HadDep_timesep_1D;
      delete[] ProtonDep_timesep_2D;
      delete[] ProtonDep_timesep_1D;
      delete[] NeutronDep_timesep_2D;
      delete[] NeutronDep_timesep_1D;
      delete[] PiCDep_timesep_2D;
      delete[] PiCDep_timesep_1D;
      delete[] Pi0Dep_timesep_2D;
      delete[] Pi0Dep_timesep_1D;
      delete[] OtherDep_timesep_2D;
      delete[] OtherDep_timesep_1D;
      delete[] LepDaughterDep_timesep_2D;
      delete[] LepDaughterDep_timesep_1D;
      delete[] HadDaughterDep_timesep_2D;
      delete[] HadDaughterDep_timesep_1D;
      delete[] ProtonDaughterDep_timesep_2D;
      delete[] ProtonDaughterDep_timesep_1D;
      delete[] NeutronDaughterDep_timesep_2D;
      delete[] NeutronDaughterDep_timesep_1D;
      delete[] PiCDaughterDep_timesep_2D;
      delete[] PiCDaughterDep_timesep_1D;
      delete[] Pi0DaughterDep_timesep_2D;
      delete[] Pi0DaughterDep_timesep_1D;
      delete[] OtherDaughterDep_timesep_2D;
      delete[] OtherDaughterDep_timesep_1D;
      delete[] LepDep_timesep;
      delete[] HadDep_timesep;
      delete[] ProtonDep_timesep;
      delete[] NeutronDep_timesep;
      delete[] PiCDep_timesep;
      delete[] Pi0Dep_timesep;
      delete[] OtherDep_timesep;
      delete[] LepDaughterDep_timesep;
      delete[] HadDaughterDep_timesep;
      delete[] ProtonDaughterDep_timesep;
      delete[] NeutronDaughterDep_timesep;
      delete[] PiCDaughterDep_timesep;
      delete[] Pi0DaughterDep_timesep;
      delete[] OtherDaughterDep_timesep;
    }

    delete[] MuonTrackPos_1D;
    delete[] MuonTrackMom_1D;
    delete[] MuonTrackPos;
    delete[] MuonTrackMom;
  }

  CondensedDeposits *CondensedDeposits::MakeTreeWriter(Int_t NXBins,
                                           Int_t NMaxTrackSteps,
                                           Double_t timesep_us) {
    CondensedDeposits *fdr = new CondensedDeposits();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->timesep_us = timesep_us;

    fdr->NXBins = NXBins;
    fdr->NMaxTrackSteps = NMaxTrackSteps;

    fdr->AllocateArrays();

    fdr->tree->Branch("EventNum", &fdr->EventNum, "EventNum/I");
    fdr->tree->Branch("nu_4mom", &fdr->nu_4mom, "nu_4mom[4]/D");
    fdr->tree->Branch("nu_PDG", &fdr->nu_PDG, "nu_PDG/I");
    fdr->tree->Branch("VertexPosition", &fdr->VertexPosition, "VertexPosition[3]/D");
    fdr->tree->Branch("EventCode", &fdr->EventCode);
    fdr->tree->Branch("PrimaryLepPDG", &fdr->PrimaryLepPDG, "PrimaryLepPDG/I");
    fdr->tree->Branch("PrimaryLep_4mom", &fdr->PrimaryLep_4mom,
                 "PrimaryLep_4mom[4]/D");
    fdr->tree->Branch("Q2_True", &fdr->Q2_True, "Q2_True/D");
    fdr->tree->Branch("FourMomTransfer_True", &fdr->FourMomTransfer_True,
                 "FourMomTransfer_True[4]/D");
    fdr->tree->Branch("y_True", &fdr->y_True, "y_True/D");
    fdr->tree->Branch("W_Rest", &fdr->W_Rest, "W_Rest/D");
    fdr->tree->Branch("NFSParts", &fdr->NFSParts, "NFSParts/I");
    fdr->tree->Branch("FSPart_PDG", &fdr->FSPart_PDG, "FSPart_PDG[NFSParts]/I");
    fdr->tree->Branch("NFSPart4MomEntries", &fdr->NFSPart4MomEntries,
                 "NFSPart4MomEntries/I");
    fdr->tree->Branch("FSPart_4Mom", &fdr->FSPart_4Mom,
                 "FSPart_4Mom[NFSPart4MomEntries]/D");
    fdr->tree->Branch("NLep", &fdr->NLep, "NLep/I");
    fdr->tree->Branch("NPi0", &fdr->NPi0, "NPi0/I");
    fdr->tree->Branch("NPiC", &fdr->NPiC, "NPiC/I");
    fdr->tree->Branch("NProton", &fdr->NProton, "NProton/I");
    fdr->tree->Branch("NNeutron", &fdr->NNeutron, "NNeutron/I");
    fdr->tree->Branch("NGamma", &fdr->NGamma, "NGamma/I");
    fdr->tree->Branch("NOther", &fdr->NOther, "NOther/I");
    fdr->tree->Branch("NBaryonicRes", &fdr->NBaryonicRes, "NBaryonicRes/I");
    fdr->tree->Branch("NAntiNucleons", &fdr->NAntiNucleons, "NAntiNucleons/I");
    fdr->tree->Branch("EKinPi0_True", &fdr->EKinPi0_True, "EKinPi0_True/D");
    fdr->tree->Branch("EMassPi0_True", &fdr->EMassPi0_True, "EMassPi0_True/D");
    fdr->tree->Branch("EKinPiC_True", &fdr->EKinPiC_True, "EKinPiC_True/D");
    fdr->tree->Branch("EMassPiC_True", &fdr->EMassPiC_True, "EMassPiC_True/D");
    fdr->tree->Branch("EKinProton_True", &fdr->EKinProton_True, "EKinProton_True/D");
    fdr->tree->Branch("EMassProton_True", &fdr->EMassProton_True,
                 "EMassProton_True/D");
    fdr->tree->Branch("EKinNeutron_True", &fdr->EKinNeutron_True,
                 "EKinNeutron_True/D");
    fdr->tree->Branch("EMassNeutron_True", &fdr->EMassNeutron_True,
                 "EMassNeutron_True/D");
    fdr->tree->Branch("EGamma_True", &fdr->EGamma_True, "EGamma_True/D");
    fdr->tree->Branch("EOther_True", &fdr->EOther_True, "EOther_True/D");
    fdr->tree->Branch("ENonPrimaryLep_True", &fdr->ENonPrimaryLep_True,
                 "ENonPrimaryLep_True/D");
    fdr->tree->Branch("KENuclearRemnant_True", &fdr->KENuclearRemnant_True,
                 "KENuclearRemnant_True/D");

    fdr->tree->Branch("TotalFS_3mom", &fdr->TotalFS_3mom, "TotalFS_3mom[3]/D");
    fdr->tree->Branch(
        "LepDep", fdr->LepDep_1D,
        (std::string("LepDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "HadDep", fdr->HadDep_1D,
        (std::string("HadDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "ProtonDep", fdr->ProtonDep_1D,
        (std::string("ProtonDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "NeutronDep", fdr->NeutronDep_1D,
        (std::string("NeutronDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "NeutronDep_ChrgWSumTime", fdr->NeutronDep_ChrgWSumTime_1D,
        (std::string("NeutronDep_ChrgWSumTime[") + to_str(NXBins) + "][3][3]/D")
            .c_str());
    fdr->tree->Branch(
        "PiCDep", fdr->PiCDep_1D,
        (std::string("PiCDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "Pi0Dep", fdr->Pi0Dep_1D,
        (std::string("Pi0Dep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch(
        "OtherDep", fdr->OtherDep_1D,
        (std::string("OtherDep[") + to_str(NXBins) + "][3][3]/D").c_str());
    fdr->tree->Branch("LepDaughterDep", fdr->LepDaughterDep_1D,
                 (std::string("LepDaughterDep[") + to_str(NXBins) + "][3][3]/D")
                     .c_str());
    fdr->tree->Branch("HadDaughterDep", fdr->HadDaughterDep_1D,
                 (std::string("HadDaughterDep[") + to_str(NXBins) + "][3][3]/D")
                     .c_str());
    fdr->tree->Branch(
        "ProtonDaughterDep", fdr->ProtonDaughterDep_1D,
        (std::string("ProtonDaughterDep[") + to_str(NXBins) + "][3][3]/D")
            .c_str());
    fdr->tree->Branch(
        "NeutronDaughterDep", fdr->NeutronDaughterDep_1D,
        (std::string("NeutronDaughterDep[") + to_str(NXBins) + "][3][3]/D")
            .c_str());
    fdr->tree->Branch("NeutronDaughterDep_ChrgWSumTime",
                 fdr->NeutronDaughterDep_ChrgWSumTime_1D,
                 (std::string("NeutronDaughterDep_ChrgWSumTime[") +
                  to_str(NXBins) + "][3][3]/D")
                     .c_str());
    fdr->tree->Branch("PiCDaughterDep", fdr->PiCDaughterDep_1D,
                 (std::string("PiCDaughterDep[") + to_str(NXBins) + "][3][3]/D")
                     .c_str());
    fdr->tree->Branch("Pi0DaughterDep", fdr->Pi0DaughterDep_1D,
                 (std::string("Pi0DaughterDep[") + to_str(NXBins) + "][3][3]/D")
                     .c_str());
    fdr->tree->Branch(
        "OtherDaughterDep", fdr->OtherDaughterDep_1D,
        (std::string("OtherDaughterDep[") + to_str(NXBins) + "][3][3]/D")
            .c_str());

    fdr->tree->Branch(
        "NuclRemDep", fdr->NuclRemDep_1D,
        (std::string("NuclRemDep[") + to_str(NXBins) + "][3][3]/D").c_str());

    if (fdr->timesep_us != 0xdeadbeef) {
      fdr->tree->Branch(
          "LepDep_timesep", fdr->LepDep_timesep_1D,
          (std::string("LepDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "HadDep_timesep", fdr->HadDep_timesep_1D,
          (std::string("HadDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "ProtonDep_timesep", fdr->ProtonDep_timesep_1D,
          (std::string("ProtonDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "NeutronDep_timesep", fdr->NeutronDep_timesep_1D,
          (std::string("NeutronDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "PiCDep_timesep", fdr->PiCDep_timesep_1D,
          (std::string("PiCDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "Pi0Dep_timesep", fdr->Pi0Dep_timesep_1D,
          (std::string("Pi0Dep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch(
          "OtherDep_timesep", fdr->OtherDep_timesep_1D,
          (std::string("OtherDep_timesep[") + to_str(NXBins) + "][3][3]/D")
              .c_str());
      fdr->tree->Branch("LepDaughterDep_timesep", fdr->LepDaughterDep_timesep_1D,
                   (std::string("LepDaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("HadDaughterDep_timesep", fdr->HadDaughterDep_timesep_1D,
                   (std::string("HadDaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("ProtonDaughterDep_timesep",
                   fdr->ProtonDaughterDep_timesep_1D,
                   (std::string("ProtonDaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("NeutronDaughterDep_timesep",
                   fdr->NeutronDaughterDep_timesep_1D,
                   (std::string("NeutronDaughterDep_timesep[") +
                    to_str(NXBins) + "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("PiCDaughterDep_timesep", fdr->PiCDaughterDep_timesep_1D,
                   (std::string("PiCDaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("Pi0DaughterDep_timesep", fdr->Pi0DaughterDep_timesep_1D,
                   (std::string("Pi0DaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
      fdr->tree->Branch("OtherDaughterDep_timesep", fdr->OtherDaughterDep_timesep_1D,
                   (std::string("OtherDaughterDep_timesep[") + to_str(NXBins) +
                    "][3][3]/D")
                       .c_str());
    }

    fdr->tree->Branch("NMuonTrackSteps", &fdr->NMuonTrackSteps, "NMuonTrackSteps/I");

    fdr->tree->Branch(
        "MuonTrackPos", fdr->MuonTrackPos_1D,
        (std::string("MuonTrackPos[") + to_str(NMaxTrackSteps) + "][3]/D")
            .c_str());

    fdr->tree->Branch(
        "MuonTrackMom", fdr->MuonTrackMom_1D,
        (std::string("MuonTrackMom[") + to_str(NMaxTrackSteps) + "][3]/D")
            .c_str());
    fdr->Reset();
    return fdr;
  }

  CondensedDeposits::~CondensedDeposits() {
    if(TreeOwned){
      DeAllocateArrays();
    }
  }
