#include "DepositsSummaryTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>
#include <utility>
#include <iostream>

  DepositsSummary::DepositsSummary() : tree(nullptr), timesep_us(0xdeadbeef), EventCode(nullptr) {}
  DepositsSummary::DepositsSummary(std::string const &treeName, std::string const &inputFiles,
       double timesep_us, bool IsLite)
      : DepositsSummary() {

    tree = OpenTChainWithFileList(treeName, inputFiles, NFiles);

    this->timesep_us = timesep_us;
    this->IsLite = IsLite;

    NEntries = tree->GetEntries();
    SetBranchAddresses();
    GetEntry(0);
  }
  void DepositsSummary::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t DepositsSummary::GetEntry() { return CEnt; }
  UInt_t DepositsSummary::GetEntries() { return NEntries; }

  DepositsSummary::~DepositsSummary() { delete tree; };

  size_t DepositsSummary::GetNPassthroughParts() { return NFSParts; }

  std::pair<Int_t, Double_t *> DepositsSummary::GetPassthroughPart(size_t i) {
    if (i >= GetNPassthroughParts()) {
      std::cout << "[ERROR]: Asked for passthrough part " << i
                << ", but we only have " << GetNPassthroughParts()
                << " in this event." << std::endl;
      throw;
    }

    return std::pair<Int_t, Double_t *>(FSPart_PDG[i], &FSPart_4Mom[i * 4]);
  }

  bool DepositsSummary::AddPassthroughPart(Int_t PDG, Double_t *fourmom) {
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

  void DepositsSummary::Reset() {
    stop = 0;
    stop_weight = 1;
    POT_weight = 1;
    std::fill_n(vtx, 3, 0);
    vtxInDetX = 0;
    XOffset = 0;
    GENIEInteractionTopology = 0;
    std::fill_n(nu_4mom, 4, 0);
    y_True = 0;
    Q2_True = 0;
    std::fill_n(FourMomTransfer_True, 4, 0);
    W_Rest = 0;
    nu_PDG = 0;
    PrimaryLepPDG = 0;
    std::fill_n(PrimaryLep_4mom, 4, 0);
    NFSParts = 0;
    std::fill_n(FSPart_PDG, kNMaxPassthroughParts, 0);
    NFSPart4MomEntries = 0;
    std::fill_n(FSPart_4Mom, kNMaxPassthroughParts * 4, 0);
    NLep = 0;
    NPi0 = 0;
    NPiC = 0;
    NProton = 0;
    NNeutron = 0;
    NAntiNucleons = 0;
    NGamma = 0;
    NOther = 0;
    NBaryonicRes = 0;
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
    Total_ENonPrimaryLep_True = 0;
    std::fill_n(TotalFS_3mom, 3, 0);
    ENonPrimaryLep_KinNucleonTotalOther_True = 0;
    ERecProxy_True = 0;
    LepDep_FV = 0;
    LepDep_veto = 0;
    LepDepDescendent_FV = 0;
    LepDepDescendent_veto = 0;
    ProtonDep_FV = 0;
    ProtonDep_veto = 0;
    NeutronDep_FV = 0;
    NeutronDep_ChrgWAvgTime_FV = 0;
    NeutronDep_veto = 0;
    NeutronDep_ChrgWAvgTime_veto = 0;
    PiCDep_FV = 0;
    PiCDep_veto = 0;
    Pi0Dep_FV = 0;
    Pi0Dep_veto = 0;
    OtherDep_FV = 0;
    OtherDep_veto = 0;
    TotalNonlep_Dep_FV = 0;
    TotalNonlep_Dep_veto = 0;
    LepDep_timesep_FV = 0;
    LepDep_timesep_veto = 0;
    LepDepDescendent_timesep_FV = 0;
    LepDepDescendent_timesep_veto = 0;
    ProtonDep_timesep_FV = 0;
    ProtonDep_timesep_veto = 0;
    NeutronDep_timesep_FV = 0;
    NeutronDep_timesep_veto = 0;
    PiCDep_timesep_FV = 0;
    PiCDep_timesep_veto = 0;
    Pi0Dep_timesep_FV = 0;
    Pi0Dep_timesep_veto = 0;
    OtherDep_timesep_FV = 0;
    OtherDep_timesep_veto = 0;
    TotalNonlep_Dep_timesep_FV = 0;
    TotalNonlep_Dep_timesep_veto = 0;
    LepExit = 0;
    LepExit_AboveThresh = 0;
    LepExitBack = 0;
    LepExitFront = 0;
    LepExitYLow = 0;
    LepExitYHigh = 0;
    LepExitXLow = 0;
    LepExitXHigh = 0;
    LepExitTopology = 0;
    std::fill_n(LepExitingPos, 3, 0);
    std::fill_n(LepExitingMom, 3, 0);
    LepExitKE = 0;
    IsNumu = 0;
    IsAntinu = 0;
    IsCC = 0;
    Is0Pi = 0;
    Is1PiC = 0;
    Is1Pi0 = 0;
    Is1Pi = 0;
    IsNPi = 0;
    IsOther = 0;
    Topology = 0;
    HadrShowerContainedInFV = 0;
    PrimaryLeptonContainedInFV = 0;
  }

  void DepositsSummary::Copy(DepositsSummary const &other) {
    stop = other.stop;
    stop_weight = other.stop_weight;
    POT_weight = other.POT_weight;
    std::copy_n(other.vtx, 3, vtx);
    vtxInDetX = other.vtxInDetX;
    XOffset = other.XOffset;
    GENIEInteractionTopology = other.GENIEInteractionTopology;
    std::copy_n(other.nu_4mom, 4, nu_4mom);
    y_True = other.y_True;
    Q2_True = other.Q2_True;
    std::copy_n(other.FourMomTransfer_True, 4, FourMomTransfer_True);
    W_Rest = other.W_Rest;
    nu_PDG = other.nu_PDG;
    PrimaryLepPDG = other.PrimaryLepPDG;
    std::copy_n(other.PrimaryLep_4mom, 4, PrimaryLep_4mom);
    NFSParts = other.NFSParts;
    std::copy_n(other.FSPart_PDG, kNMaxPassthroughParts, FSPart_PDG);
    NFSPart4MomEntries = other.NFSPart4MomEntries;
    std::copy_n(other.FSPart_4Mom, kNMaxPassthroughParts * 4, FSPart_4Mom);
    NLep = other.NLep;
    NPi0 = other.NPi0;
    NPiC = other.NPiC;
    NProton = other.NProton;
    NNeutron = other.NNeutron;
    NAntiNucleons = other.NAntiNucleons;
    NGamma = other.NGamma;
    NOther = other.NOther;
    NBaryonicRes = other.NBaryonicRes;
    EKinPi0_True = other.EKinPi0_True;
    EMassPi0_True = other.EMassPi0_True;
    EKinPiC_True = other.EKinPiC_True;
    EMassPiC_True = other.EMassPiC_True;
    EKinProton_True = other.EKinProton_True;
    EMassProton_True = other.EMassProton_True;
    EKinNeutron_True = other.EKinNeutron_True;
    EMassNeutron_True = other.EMassNeutron_True;
    EGamma_True = other.EGamma_True;
    EOther_True = other.EOther_True;
    Total_ENonPrimaryLep_True = other.Total_ENonPrimaryLep_True;
    std::copy_n(other.TotalFS_3mom, 3, TotalFS_3mom);
    ENonPrimaryLep_KinNucleonTotalOther_True =
        other.ENonPrimaryLep_KinNucleonTotalOther_True;
    ERecProxy_True = other.ERecProxy_True;
    LepDep_FV = other.LepDep_FV;
    LepDep_veto = other.LepDep_veto;
    LepDepDescendent_FV = other.LepDepDescendent_FV;
    LepDepDescendent_veto = other.LepDepDescendent_veto;
    ProtonDep_FV = other.ProtonDep_FV;
    ProtonDep_veto = other.ProtonDep_veto;
    NeutronDep_FV = other.NeutronDep_FV;
    NeutronDep_ChrgWAvgTime_FV = other.NeutronDep_ChrgWAvgTime_FV;
    NeutronDep_veto = other.NeutronDep_veto;
    NeutronDep_ChrgWAvgTime_veto = other.NeutronDep_ChrgWAvgTime_veto;
    PiCDep_FV = other.PiCDep_FV;
    PiCDep_veto = other.PiCDep_veto;
    Pi0Dep_FV = other.Pi0Dep_FV;
    Pi0Dep_veto = other.Pi0Dep_veto;
    OtherDep_FV = other.OtherDep_FV;
    OtherDep_veto = other.OtherDep_veto;
    TotalNonlep_Dep_FV = other.TotalNonlep_Dep_FV;
    TotalNonlep_Dep_veto = other.TotalNonlep_Dep_veto;
    LepDep_timesep_FV = other.LepDep_timesep_FV;
    LepDep_timesep_veto = other.LepDep_timesep_veto;
    LepDepDescendent_timesep_FV = other.LepDepDescendent_timesep_FV;
    LepDepDescendent_timesep_veto = other.LepDepDescendent_timesep_veto;
    ProtonDep_timesep_FV = other.ProtonDep_timesep_FV;
    ProtonDep_timesep_veto = other.ProtonDep_timesep_veto;
    NeutronDep_timesep_FV = other.NeutronDep_timesep_FV;
    NeutronDep_timesep_veto = other.NeutronDep_timesep_veto;
    PiCDep_timesep_FV = other.PiCDep_timesep_FV;
    PiCDep_timesep_veto = other.PiCDep_timesep_veto;
    Pi0Dep_timesep_FV = other.Pi0Dep_timesep_FV;
    Pi0Dep_timesep_veto = other.Pi0Dep_timesep_veto;
    OtherDep_timesep_FV = other.OtherDep_timesep_FV;
    OtherDep_timesep_veto = other.OtherDep_timesep_veto;
    TotalNonlep_Dep_timesep_FV = other.TotalNonlep_Dep_timesep_FV;
    TotalNonlep_Dep_timesep_veto = other.TotalNonlep_Dep_timesep_veto;
    LepExit = other.LepExit;
    LepExit_AboveThresh = other.LepExit_AboveThresh;
    LepExitBack = other.LepExitBack;
    LepExitFront = other.LepExitFront;
    LepExitYLow = other.LepExitYLow;
    LepExitYHigh = other.LepExitYHigh;
    LepExitXLow = other.LepExitXLow;
    LepExitXHigh = other.LepExitXHigh;
    LepExitTopology = other.LepExitTopology;
    std::copy_n(other.LepExitingPos, 3, LepExitingPos);
    std::copy_n(other.LepExitingMom, 3, LepExitingMom);
    LepExitKE = other.LepExitKE;
    IsNumu = other.IsNumu;
    IsAntinu = other.IsAntinu;
    IsCC = other.IsCC;
    Is0Pi = other.Is0Pi;
    Is1PiC = other.Is1PiC;
    Is1Pi0 = other.Is1Pi0;
    Is1Pi = other.Is1Pi;
    IsNPi = other.IsNPi;
    IsOther = other.IsOther;
    Topology = other.Topology;
    HadrShowerContainedInFV = other.HadrShowerContainedInFV;
    PrimaryLeptonContainedInFV = other.PrimaryLeptonContainedInFV;
  }

  DepositsSummary *DepositsSummary::MakeTreeWriter(TTree *OutputTree, double timesep_us,
                              bool IsLite) {
    DepositsSummary *rtn = new DepositsSummary();
    rtn->timesep_us = timesep_us;
    rtn->IsLite = IsLite;

    OutputTree->Branch("stop", &rtn->stop, "stop/I");
    OutputTree->Branch("stop_weight", &rtn->stop_weight, "stop_weight/D");
    OutputTree->Branch("POT_weight", &rtn->POT_weight, "POT_weight/D");

    OutputTree->Branch("vtx", &rtn->vtx, "vtx[3]/D");
    OutputTree->Branch("vtxInDetX", &rtn->vtxInDetX, "vtxInDetX/D");
    OutputTree->Branch("XOffset", &rtn->XOffset, "XOffset/D");
    OutputTree->Branch("EventCode", &rtn->EventCode);
    OutputTree->Branch("GENIEInteractionTopology",
                       &rtn->GENIEInteractionTopology);

    OutputTree->Branch("nu_4mom", &rtn->nu_4mom, "nu_4mom[4]/D");
    if (!IsLite) {
      OutputTree->Branch("y_True", &rtn->y_True, "y_True/D");
      OutputTree->Branch("W_Rest", &rtn->W_Rest, "W_Rest/D");
      OutputTree->Branch("Q2_True", &rtn->Q2_True, "Q2_True/D");
      OutputTree->Branch("FourMomTransfer_True", &rtn->FourMomTransfer_True,
                         "FourMomTransfer_True[4]/D");
    }

    OutputTree->Branch("nu_PDG", &rtn->nu_PDG, "nu_PDG/I");
    OutputTree->Branch("PrimaryLepPDG", &rtn->PrimaryLepPDG, "PrimaryLepPDG/I");
    OutputTree->Branch("PrimaryLep_4mom", &rtn->PrimaryLep_4mom,
                       "PrimaryLep_4mom[4]/D");

    if (!IsLite) {
      OutputTree->Branch("NFSParts", &rtn->NFSParts, "NFSParts/I");
      OutputTree->Branch("FSPart_PDG", &rtn->FSPart_PDG,
                         "FSPart_PDG[NFSParts]/I");
      OutputTree->Branch("NFSPart4MomEntries", &rtn->NFSPart4MomEntries,
                         "NFSPart4MomEntries/I");
      OutputTree->Branch("FSPart_4Mom", &rtn->FSPart_4Mom,
                         "FSPart_4Mom[NFSPart4MomEntries]/D");

      OutputTree->Branch("NLep", &rtn->NLep, "NLep/I");
      OutputTree->Branch("NPi0", &rtn->NPi0, "NPi0/I");
      OutputTree->Branch("NPiC", &rtn->NPiC, "NPiC/I");
      OutputTree->Branch("NProton", &rtn->NProton, "NProton/I");
      OutputTree->Branch("NNeutron", &rtn->NNeutron, "NNeutron/I");
      OutputTree->Branch("NGamma", &rtn->NGamma, "NGamma/I");
      OutputTree->Branch("NOther", &rtn->NOther, "NOther/I");
      OutputTree->Branch("NBaryonicRes", &rtn->NBaryonicRes, "NBaryonicRes/I");
      OutputTree->Branch("NAntiNucleons", &rtn->NAntiNucleons,
                         "NAntiNucleons/I");

      OutputTree->Branch("EKinPi0_True", &rtn->EKinPi0_True, "EKinPi0_True/D");
      OutputTree->Branch("EMassPi0_True", &rtn->EMassPi0_True,
                         "EMassPi0_True/D");
      OutputTree->Branch("EKinPiC_True", &rtn->EKinPiC_True, "EKinPiC_True/D");
      OutputTree->Branch("EMassPiC_True", &rtn->EMassPiC_True,
                         "EMassPiC_True/D");
      OutputTree->Branch("EKinProton_True", &rtn->EKinProton_True,
                         "EKinProton_True/D");
      OutputTree->Branch("EMassProton_True", &rtn->EMassProton_True,
                         "EMassProton_True/D");
      OutputTree->Branch("EKinNeutron_True", &rtn->EKinNeutron_True,
                         "EKinNeutron_True/D");
      OutputTree->Branch("EMassNeutron_True", &rtn->EMassNeutron_True,
                         "EMassNeutron_True/D");
      OutputTree->Branch("EGamma_True", &rtn->EGamma_True, "EGamma_True/D");
      OutputTree->Branch("EOther_True", &rtn->EOther_True, "EOther_True/D");
    }

    OutputTree->Branch("Total_ENonPrimaryLep_True",
                       &rtn->Total_ENonPrimaryLep_True,
                       "Total_ENonPrimaryLep_True/D");
    OutputTree->Branch("TotalFS_3mom", &rtn->TotalFS_3mom, "TotalFS_3mom[3]/D");
    OutputTree->Branch("ENonPrimaryLep_KinNucleonTotalOther_True",
                       &rtn->ENonPrimaryLep_KinNucleonTotalOther_True,
                       "ENonPrimaryLep_KinNucleonTotalOther_True/D");
    OutputTree->Branch("ERecProxy_True", &rtn->ERecProxy_True,
                       "ERecProxy_True/D");

    OutputTree->Branch("LepDep_FV", &rtn->LepDep_FV, "LepDep_FV/D");
    OutputTree->Branch("LepDep_veto", &rtn->LepDep_veto, "LepDep_veto/D");
    OutputTree->Branch("LepDepDescendent_FV", &rtn->LepDepDescendent_FV,
                       "LepDepDescendent_FV/D");
    OutputTree->Branch("LepDepDescendent_veto", &rtn->LepDepDescendent_veto,
                       "LepDepDescendent_veto/D");
    OutputTree->Branch("ProtonDep_FV", &rtn->ProtonDep_FV, "ProtonDep_FV/D");
    OutputTree->Branch("ProtonDep_veto", &rtn->ProtonDep_veto,
                       "ProtonDep_veto/D");
    OutputTree->Branch("NeutronDep_FV", &rtn->NeutronDep_FV, "NeutronDep_FV/D");
    OutputTree->Branch("NeutronDep_ChrgWAvgTime_FV",
                       &rtn->NeutronDep_ChrgWAvgTime_FV,
                       "NeutronDep_ChrgWAvgTime_FV/D");
    OutputTree->Branch("NeutronDep_veto", &rtn->NeutronDep_veto,
                       "NeutronDep_veto/D");
    OutputTree->Branch("NeutronDep_ChrgWAvgTime_veto",
                       &rtn->NeutronDep_ChrgWAvgTime_veto,
                       "NeutronDep_ChrgWAvgTime_veto/D");
    OutputTree->Branch("PiCDep_FV", &rtn->PiCDep_FV, "PiCDep_FV/D");
    OutputTree->Branch("PiCDep_veto", &rtn->PiCDep_veto, "PiCDep_veto/D");
    OutputTree->Branch("Pi0Dep_FV", &rtn->Pi0Dep_FV, "Pi0Dep_FV/D");
    OutputTree->Branch("Pi0Dep_veto", &rtn->Pi0Dep_veto, "Pi0Dep_veto/D");
    OutputTree->Branch("OtherDep_FV", &rtn->OtherDep_FV, "OtherDep_FV/D");
    OutputTree->Branch("OtherDep_veto", &rtn->OtherDep_veto, "OtherDep_veto/D");

    OutputTree->Branch("TotalNonlep_Dep_FV", &rtn->TotalNonlep_Dep_FV,
                       "TotalNonlep_Dep_FV/D");
    OutputTree->Branch("TotalNonlep_Dep_veto", &rtn->TotalNonlep_Dep_veto,
                       "TotalNonlep_Dep_veto/D");

    if (rtn->timesep_us != 0xdeadbeef) {
      OutputTree->Branch("LepDep_timesep_FV", &rtn->LepDep_timesep_FV,
                         "LepDep_timesep_FV/D");
      OutputTree->Branch("LepDep_timesep_veto", &rtn->LepDep_timesep_veto,
                         "LepDep_timesep_veto/D");
      OutputTree->Branch("LepDepDescendent_timesep_FV",
                         &rtn->LepDepDescendent_timesep_FV,
                         "LepDepDescendent_timesep_FV/D");
      OutputTree->Branch("LepDepDescendent_timesep_veto",
                         &rtn->LepDepDescendent_timesep_veto,
                         "LepDepDescendent_timesep_veto/D");
      OutputTree->Branch("ProtonDep_timesep_FV", &rtn->ProtonDep_timesep_FV,
                         "ProtonDep_timesep_FV/D");
      OutputTree->Branch("ProtonDep_timesep_veto", &rtn->ProtonDep_timesep_veto,
                         "ProtonDep_timesep_veto/D");
      OutputTree->Branch("NeutronDep_timesep_FV", &rtn->NeutronDep_timesep_FV,
                         "NeutronDep_timesep_FV/D");
      OutputTree->Branch("NeutronDep_timesep_veto",
                         &rtn->NeutronDep_timesep_veto,
                         "NeutronDep_timesep_veto/D");
      OutputTree->Branch("PiCDep_timesep_FV", &rtn->PiCDep_timesep_FV,
                         "PiCDep_timesep_FV/D");
      OutputTree->Branch("PiCDep_timesep_veto", &rtn->PiCDep_timesep_veto,
                         "PiCDep_timesep_veto/D");
      OutputTree->Branch("Pi0Dep_timesep_FV", &rtn->Pi0Dep_timesep_FV,
                         "Pi0Dep_timesep_FV/D");
      OutputTree->Branch("Pi0Dep_timesep_veto", &rtn->Pi0Dep_timesep_veto,
                         "Pi0Dep_timesep_veto/D");
      OutputTree->Branch("OtherDep_timesep_FV", &rtn->OtherDep_timesep_FV,
                         "OtherDep_timesep_FV/D");
      OutputTree->Branch("OtherDep_timesep_veto", &rtn->OtherDep_timesep_veto,
                         "OtherDep_timesep_veto/D");

      OutputTree->Branch("TotalNonlep_Dep_timesep_FV",
                         &rtn->TotalNonlep_Dep_timesep_FV,
                         "TotalNonlep_Dep_timesep_FV/D");
      OutputTree->Branch("TotalNonlep_Dep_timesep_veto",
                         &rtn->TotalNonlep_Dep_timesep_veto,
                         "TotalNonlep_Dep_timesep_veto/D");
    }

    OutputTree->Branch("LepExit", &rtn->LepExit, "LepExit/O");
    OutputTree->Branch("LepExit_AboveThresh", &rtn->LepExit_AboveThresh,
                       "LepExit_AboveThresh/O");
    if (!IsLite) {
      OutputTree->Branch("LepExitBack", &rtn->LepExitBack, "LepExitBack/O");
      OutputTree->Branch("LepExitFront", &rtn->LepExitFront, "LepExitFront/O");
      OutputTree->Branch("LepExitYLow", &rtn->LepExitYLow, "LepExitYLow/O");
      OutputTree->Branch("LepExitYHigh", &rtn->LepExitYHigh, "LepExitYHigh/O");
      OutputTree->Branch("LepExitXLow", &rtn->LepExitXLow, "LepExitXLow/O");
      OutputTree->Branch("LepExitXHigh", &rtn->LepExitXHigh, "LepExitXHigh/O");
    }

    OutputTree->Branch("LepExitTopology", &rtn->LepExitTopology,
                       "LepExitTopology/I");

    OutputTree->Branch("LepExitKE", &rtn->LepExitKE, "LepExitKE/D");
    OutputTree->Branch("LepExitingPos", &rtn->LepExitingPos,
                       "LepExitingPos[3]/D");
    OutputTree->Branch("LepExitingMom", &rtn->LepExitingMom,
                       "LepExitingMom[3]/D");

    OutputTree->Branch("IsNumu", &rtn->IsNumu, "IsNumu/O");
    OutputTree->Branch("IsAntinu", &rtn->IsAntinu, "IsAntinu/O");
    OutputTree->Branch("IsCC", &rtn->IsCC, "IsCC/O");
    OutputTree->Branch("Is0Pi", &rtn->Is0Pi, "Is0Pi/O");
    OutputTree->Branch("Is1PiC", &rtn->Is1PiC, "Is1PiC/O");
    OutputTree->Branch("Is1Pi0", &rtn->Is1Pi0, "Is1Pi0/O");
    OutputTree->Branch("Is1Pi", &rtn->Is1Pi, "Is1Pi/O");
    OutputTree->Branch("IsNPi", &rtn->IsNPi, "IsNPi/O");
    OutputTree->Branch("IsOther", &rtn->IsOther, "IsOther/O");
    OutputTree->Branch("Topology", &rtn->Topology, "Topology/I");
    OutputTree->Branch("HadrShowerContainedInFV", &rtn->HadrShowerContainedInFV,
                       "HadrShowerContainedInFV/O");
    OutputTree->Branch("PrimaryLeptonContainedInFV",
                       &rtn->PrimaryLeptonContainedInFV,
                       "PrimaryLeptonContainedInFV/O");
    rtn->Reset();
    return rtn;
  }

  void DepositsSummary::SetBranchAddresses() {
    if (!CheckTTreeHasBranch(tree, "y_True") && !IsLite) {
      std::cout << "[INFO]: Switching EDepReader to LiteMode as missing "
                   "expected branch."
                << std::endl;
      IsLite = true;
    }

    tree->SetBranchAddress("stop", &stop);
    tree->SetBranchAddress("stop_weight", &stop_weight);
    tree->SetBranchAddress("POT_weight", &POT_weight);

    tree->SetBranchAddress("vtx", &vtx);
    tree->SetBranchAddress("vtxInDetX", &vtxInDetX);
    tree->SetBranchAddress("XOffset", &XOffset);
    tree->SetBranchAddress("EventCode", &EventCode);
    tree->SetBranchAddress("GENIEInteractionTopology",
                           &GENIEInteractionTopology);
    tree->SetBranchAddress("nu_4mom", &nu_4mom);
    if (!IsLite) {
      tree->SetBranchAddress("y_True", &y_True);
      tree->SetBranchAddress("Q2_True", &Q2_True);
      tree->SetBranchAddress("FourMomTransfer_True", &FourMomTransfer_True);
      tree->SetBranchAddress("W_Rest", &W_Rest);
    }
    tree->SetBranchAddress("nu_PDG", &nu_PDG);
    tree->SetBranchAddress("PrimaryLepPDG", &PrimaryLepPDG);
    tree->SetBranchAddress("PrimaryLep_4mom", &PrimaryLep_4mom);
    if (!IsLite) {
      tree->SetBranchAddress("NFSParts", &NFSParts);
      tree->SetBranchAddress("FSPart_PDG", FSPart_PDG);
      tree->SetBranchAddress("NFSPart4MomEntries", &NFSPart4MomEntries);
      tree->SetBranchAddress("FSPart_4Mom", FSPart_4Mom);
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
    }
    tree->SetBranchAddress("Total_ENonPrimaryLep_True",
                           &Total_ENonPrimaryLep_True);
    tree->SetBranchAddress("TotalFS_3mom", &TotalFS_3mom);
    tree->SetBranchAddress("ENonPrimaryLep_KinNucleonTotalOther_True",
                           &ENonPrimaryLep_KinNucleonTotalOther_True);
    tree->SetBranchAddress("ERecProxy_True", &ERecProxy_True);
    tree->SetBranchAddress("LepDep_FV", &LepDep_FV);
    tree->SetBranchAddress("LepDep_veto", &LepDep_veto);
    tree->SetBranchAddress("LepDepDescendent_FV", &LepDepDescendent_FV);
    tree->SetBranchAddress("LepDepDescendent_veto", &LepDepDescendent_veto);
    tree->SetBranchAddress("ProtonDep_FV", &ProtonDep_FV);
    tree->SetBranchAddress("ProtonDep_veto", &ProtonDep_veto);
    tree->SetBranchAddress("NeutronDep_FV", &NeutronDep_FV);
    tree->SetBranchAddress("NeutronDep_ChrgWAvgTime_FV",
                           &NeutronDep_ChrgWAvgTime_FV);
    tree->SetBranchAddress("NeutronDep_veto", &NeutronDep_veto);
    tree->SetBranchAddress("NeutronDep_ChrgWAvgTime_veto",
                           &NeutronDep_ChrgWAvgTime_veto);
    tree->SetBranchAddress("PiCDep_FV", &PiCDep_FV);
    tree->SetBranchAddress("PiCDep_veto", &PiCDep_veto);
    tree->SetBranchAddress("Pi0Dep_FV", &Pi0Dep_FV);
    tree->SetBranchAddress("Pi0Dep_veto", &Pi0Dep_veto);
    tree->SetBranchAddress("OtherDep_FV", &OtherDep_FV);
    tree->SetBranchAddress("OtherDep_veto", &OtherDep_veto);
    tree->SetBranchAddress("TotalNonlep_Dep_FV", &TotalNonlep_Dep_FV);
    tree->SetBranchAddress("TotalNonlep_Dep_veto", &TotalNonlep_Dep_veto);
    if (timesep_us != 0xdeadbeef) {
      tree->SetBranchAddress("LepDep_timesep_FV", &LepDep_timesep_FV);
      tree->SetBranchAddress("LepDep_timesep_veto", &LepDep_timesep_veto);
      tree->SetBranchAddress("LepDepDescendent_timesep_FV",
                             &LepDepDescendent_timesep_FV);
      tree->SetBranchAddress("LepDepDescendent_timesep_veto",
                             &LepDepDescendent_timesep_veto);
      tree->SetBranchAddress("ProtonDep_timesep_FV", &ProtonDep_timesep_FV);
      tree->SetBranchAddress("ProtonDep_timesep_veto", &ProtonDep_timesep_veto);
      tree->SetBranchAddress("NeutronDep_timesep_FV", &NeutronDep_timesep_FV);
      tree->SetBranchAddress("NeutronDep_timesep_veto",
                             &NeutronDep_timesep_veto);
      tree->SetBranchAddress("PiCDep_timesep_FV", &PiCDep_timesep_FV);
      tree->SetBranchAddress("PiCDep_timesep_veto", &PiCDep_timesep_veto);
      tree->SetBranchAddress("Pi0Dep_timesep_FV", &Pi0Dep_timesep_FV);
      tree->SetBranchAddress("Pi0Dep_timesep_veto", &Pi0Dep_timesep_veto);
      tree->SetBranchAddress("OtherDep_timesep_FV", &OtherDep_timesep_FV);
      tree->SetBranchAddress("OtherDep_timesep_veto", &OtherDep_timesep_veto);
      tree->SetBranchAddress("TotalNonlep_Dep_timesep_FV",
                             &TotalNonlep_Dep_timesep_FV);
      tree->SetBranchAddress("TotalNonlep_Dep_timesep_veto",
                             &TotalNonlep_Dep_timesep_veto);
    }
    tree->SetBranchAddress("LepExit", &LepExit);
    tree->SetBranchAddress("LepExitKE", &LepExitKE);
    tree->SetBranchAddress("LepExit_AboveThresh", &LepExit_AboveThresh);
    if (!IsLite) {
      tree->SetBranchAddress("LepExitBack", &LepExitBack);
      tree->SetBranchAddress("LepExitFront", &LepExitFront);
      tree->SetBranchAddress("LepExitYLow", &LepExitYLow);
      tree->SetBranchAddress("LepExitYHigh", &LepExitYHigh);
      tree->SetBranchAddress("LepExitXLow", &LepExitXLow);
      tree->SetBranchAddress("LepExitXHigh", &LepExitXHigh);
    }
    tree->SetBranchAddress("LepExitTopology", &LepExitTopology);
    tree->SetBranchAddress("LepExitingPos", &LepExitingPos);
    tree->SetBranchAddress("LepExitingMom", &LepExitingMom);

    tree->SetBranchAddress("IsNumu", &IsNumu);
    tree->SetBranchAddress("IsAntinu", &IsAntinu);
    tree->SetBranchAddress("IsCC", &IsCC);
    tree->SetBranchAddress("Is0Pi", &Is0Pi);
    tree->SetBranchAddress("Is1PiC", &Is1PiC);
    tree->SetBranchAddress("Is1Pi0", &Is1Pi0);
    tree->SetBranchAddress("Is1Pi", &Is1Pi);
    tree->SetBranchAddress("IsNPi", &IsNPi);
    tree->SetBranchAddress("IsOther", &IsOther);
    tree->SetBranchAddress("Topology", &Topology);
    tree->SetBranchAddress("HadrShowerContainedInFV", &HadrShowerContainedInFV);
    tree->SetBranchAddress("PrimaryLeptonContainedInFV",
                           &PrimaryLeptonContainedInFV);
  }
