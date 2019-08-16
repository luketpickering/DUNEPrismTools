#include "dk2nu_TreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>
#include <iostream>
#include <vector>

DK2NuReader::DK2NuReader(std::string treeName, std::string inputFiles,
                         bool DK2NULite)
    : ppfx_tree(nullptr), ppfx_vwgt_tot(nullptr), fUseAllPPFXBranches(false),
      ppfx_vwgt_mipp_pi(nullptr), ppfx_vwgt_mipp_K(nullptr),
      ppfx_vwgt_abs(nullptr), ppfx_vwgt_att(nullptr), ppfx_vwgt_ttpCpi(nullptr),
      ppfx_vwgt_ttpCk(nullptr), ppfx_vwgt_ttnCpi(nullptr),
      ppfx_vwgt_ttpCnu(nullptr), ppfx_vwgt_ttnua(nullptr),
      ppfx_vwgt_ttmesinc(nullptr), ppfx_vwgt_oth(nullptr) {
  tree = OpenTChainWithFileList(treeName, inputFiles, NFiles);
  NEntries = tree->GetEntries();

  SetBranchAddresses(DK2NULite);
  std::cout << "[DK2NuReader]: Loaded TChain: " << NFiles << " files and "
            << NEntries << " entries." << std::endl;
  GetEntry(0);
}

void DK2NuReader::SetBranchAddresses(bool DK2NULite) {
  if (!DK2NULite) {
    tree->SetMakeClass(true);

    tree->SetBranchAddress("job", &job);
    tree->SetBranchAddress("potnum", &potnum);
    tree->SetBranchAddress("decay.norig", &decay_norig);
    tree->SetBranchAddress("decay.ndecay", &decay_ndecay);
    tree->SetBranchAddress("decay.ntype", &decay_ntype);
    tree->SetBranchAddress("decay.vx", &decay_vx);
    tree->SetBranchAddress("decay.vy", &decay_vy);
    tree->SetBranchAddress("decay.vz", &decay_vz);
    tree->SetBranchAddress("decay.pdpx", &decay_pdpx);
    tree->SetBranchAddress("decay.pdpy", &decay_pdpy);
    tree->SetBranchAddress("decay.pdpz", &decay_pdpz);
    tree->SetBranchAddress("decay.ppdxdz", &decay_ppdxdz);
    tree->SetBranchAddress("decay.ppdydz", &decay_ppdydz);
    tree->SetBranchAddress("decay.pppz", &decay_pppz);
    tree->SetBranchAddress("decay.ppenergy", &decay_ppenergy);
    tree->SetBranchAddress("decay.ppmedium", &decay_ppmedium);
    tree->SetBranchAddress("decay.ptype", &decay_ptype);
    tree->SetBranchAddress("decay.muparpx", &decay_muparpx);
    tree->SetBranchAddress("decay.muparpy", &decay_muparpy);
    tree->SetBranchAddress("decay.muparpz", &decay_muparpz);
    tree->SetBranchAddress("decay.mupare", &decay_mupare);
    tree->SetBranchAddress("decay.necm", &decay_necm);
    tree->SetBranchAddress("decay.nimpwt", &decay_nimpwt);
    tree->SetBranchAddress("nuray", &nuray_);
    tree->SetBranchAddress("nuray.px", nuray_px);
    tree->SetBranchAddress("nuray.py", nuray_py);
    tree->SetBranchAddress("nuray.pz", nuray_pz);
    tree->SetBranchAddress("nuray.E", nuray_E);
    tree->SetBranchAddress("nuray.wgt", nuray_wgt);
    tree->SetBranchAddress("ancestor", &ancestor_);
    tree->SetBranchAddress("ancestor.pdg", ancestor_pdg);
    tree->SetBranchAddress("ancestor.startx", ancestor_startx);
    tree->SetBranchAddress("ancestor.starty", ancestor_starty);
    tree->SetBranchAddress("ancestor.startz", ancestor_startz);
    tree->SetBranchAddress("ancestor.startt", ancestor_startt);
    tree->SetBranchAddress("ancestor.startpx", ancestor_startpx);
    tree->SetBranchAddress("ancestor.startpy", ancestor_startpy);
    tree->SetBranchAddress("ancestor.startpz", ancestor_startpz);
    tree->SetBranchAddress("ancestor.stoppx", ancestor_stoppx);
    tree->SetBranchAddress("ancestor.stoppy", ancestor_stoppy);
    tree->SetBranchAddress("ancestor.stoppz", ancestor_stoppz);
    tree->SetBranchAddress("ancestor.polx", ancestor_polx);
    tree->SetBranchAddress("ancestor.poly", ancestor_poly);
    tree->SetBranchAddress("ancestor.polz", ancestor_polz);
    tree->SetBranchAddress("ancestor.pprodpx", ancestor_pprodpx);
    tree->SetBranchAddress("ancestor.pprodpy", ancestor_pprodpy);
    tree->SetBranchAddress("ancestor.pprodpz", ancestor_pprodpz);
    tree->SetBranchAddress("ancestor.nucleus", ancestor_nucleus);
    tree->SetBranchAddress("ancestor.proc", ancestor_proc);
    tree->SetBranchAddress("ancestor.ivol", ancestor_ivol);
    tree->SetBranchAddress("ancestor.imat", ancestor_imat);
    tree->SetBranchAddress("ppvx", &ppvx);
    tree->SetBranchAddress("ppvy", &ppvy);
    tree->SetBranchAddress("ppvz", &ppvz);
    tree->SetBranchAddress("tgtexit.tvx", &tgtexit_tvx);
    tree->SetBranchAddress("tgtexit.tvy", &tgtexit_tvy);
    tree->SetBranchAddress("tgtexit.tvz", &tgtexit_tvz);
    tree->SetBranchAddress("tgtexit.tpx", &tgtexit_tpx);
    tree->SetBranchAddress("tgtexit.tpy", &tgtexit_tpy);
    tree->SetBranchAddress("tgtexit.tpz", &tgtexit_tpz);
    tree->SetBranchAddress("tgtexit.tptype", &tgtexit_tptype);
    tree->SetBranchAddress("tgtexit.tgen", &tgtexit_tgen);
    tree->SetBranchAddress("traj", &traj_);
    tree->SetBranchAddress("traj.trkx", traj_trkx);
    tree->SetBranchAddress("traj.trky", traj_trky);
    tree->SetBranchAddress("traj.trkz", traj_trkz);
    tree->SetBranchAddress("traj.trkpx", traj_trkpx);
    tree->SetBranchAddress("traj.trkpy", traj_trkpy);
    tree->SetBranchAddress("traj.trkpz", traj_trkpz);
    tree->SetBranchAddress("flagbits", &flagbits);
    tree->SetBranchAddress("vint", &vint);
    tree->SetBranchAddress("vdbl", &vdbl);
  } else {
    tree->SetBranchAddress("decay_ntype", &decay_ntype);
    tree->SetBranchAddress("decay_vx", &decay_vx);
    tree->SetBranchAddress("decay_vy", &decay_vy);
    tree->SetBranchAddress("decay_vz", &decay_vz);
    tree->SetBranchAddress("decay_pdpx", &decay_pdpx);
    tree->SetBranchAddress("decay_pdpy", &decay_pdpy);
    tree->SetBranchAddress("decay_pdpz", &decay_pdpz);
    tree->SetBranchAddress("decay_ppdxdz", &decay_ppdxdz);
    tree->SetBranchAddress("decay_ppdydz", &decay_ppdydz);
    tree->SetBranchAddress("decay_pppz", &decay_pppz);
    tree->SetBranchAddress("decay_ppenergy", &decay_ppenergy);
    tree->SetBranchAddress("decay_ptype", &decay_ptype);
    tree->SetBranchAddress("decay_muparpx", &decay_muparpx);
    tree->SetBranchAddress("decay_muparpy", &decay_muparpy);
    tree->SetBranchAddress("decay_muparpz", &decay_muparpz);
    tree->SetBranchAddress("decay_mupare", &decay_mupare);
    tree->SetBranchAddress("decay_necm", &decay_necm);
    tree->SetBranchAddress("decay_nimpwt", &decay_nimpwt);
  }
}

void DK2NuReader::SetPPFXBranchAddresses(size_t NUniverses) {
  ppfx_vwgt_tot = new double[NUniverses];

  tree->SetBranchAddress("ppfx_cvwgt", &ppfx_cvwgt);
  tree->SetBranchAddress("ppfx_vwgt_tot", ppfx_vwgt_tot);

  if (CheckTTreeHasBranch(tree, "ppfx_vwgt_mipp_pi")) {
    fUseAllPPFXBranches = true;

    ppfx_vwgt_mipp_pi = new double[NUniverses];
    ppfx_vwgt_mipp_K = new double[NUniverses];
    ppfx_vwgt_abs = new double[NUniverses];
    ppfx_vwgt_att = new double[NUniverses];
    ppfx_vwgt_ttpCpi = new double[NUniverses];
    ppfx_vwgt_ttpCk = new double[NUniverses];
    ppfx_vwgt_ttnCpi = new double[NUniverses];
    ppfx_vwgt_ttpCnu = new double[NUniverses];
    ppfx_vwgt_ttnua = new double[NUniverses];
    ppfx_vwgt_ttmesinc = new double[NUniverses];
    ppfx_vwgt_oth = new double[NUniverses];

    tree->SetBranchAddress("ppfx_vwgt_mipp_pi", ppfx_vwgt_mipp_pi);
    tree->SetBranchAddress("ppfx_vwgt_mipp_K", ppfx_vwgt_mipp_K);
    tree->SetBranchAddress("ppfx_vwgt_abs", ppfx_vwgt_abs);
    tree->SetBranchAddress("ppfx_vwgt_att", ppfx_vwgt_att);
    tree->SetBranchAddress("ppfx_vwgt_ttpCpi", ppfx_vwgt_ttpCpi);
    tree->SetBranchAddress("ppfx_vwgt_ttpCk", ppfx_vwgt_ttpCk);
    tree->SetBranchAddress("ppfx_vwgt_ttnCpi", ppfx_vwgt_ttnCpi);
    tree->SetBranchAddress("ppfx_vwgt_ttpCnu", ppfx_vwgt_ttpCnu);
    tree->SetBranchAddress("ppfx_vwgt_ttnua", ppfx_vwgt_ttnua);
    tree->SetBranchAddress("ppfx_vwgt_ttmesinc", ppfx_vwgt_ttmesinc);
    tree->SetBranchAddress("ppfx_vwgt_oth", ppfx_vwgt_oth);

    tree->SetBranchAddress("ppfx_cvwgt_mipp_pi", &ppfx_cvwgt_mipp_pi);
    tree->SetBranchAddress("ppfx_cvwgt_mipp_K", &ppfx_cvwgt_mipp_K);
    tree->SetBranchAddress("ppfx_cvwgt_abs", &ppfx_cvwgt_abs);
    tree->SetBranchAddress("ppfx_cvwgt_att", &ppfx_cvwgt_att);
    tree->SetBranchAddress("ppfx_cvwgt_ttpCpi", &ppfx_cvwgt_ttpCpi);
    tree->SetBranchAddress("ppfx_cvwgt_ttpCk", &ppfx_cvwgt_ttpCk);
    tree->SetBranchAddress("ppfx_cvwgt_ttnCpi", &ppfx_cvwgt_ttnCpi);
    tree->SetBranchAddress("ppfx_cvwgt_ttpCnu", &ppfx_cvwgt_ttpCnu);
    tree->SetBranchAddress("ppfx_cvwgt_ttnua", &ppfx_cvwgt_ttnua);
    tree->SetBranchAddress("ppfx_cvwgt_ttmesinc", &ppfx_cvwgt_ttmesinc);
    tree->SetBranchAddress("ppfx_cvwgt_oth", &ppfx_cvwgt_oth);
  }
}

double DK2NuReader::GetParentMass() {
  static double const pimass = 0.13957; // in GeV
  static double const kmass = 0.49368;
  static double const k0mass = 0.49767;
  static double const mumass = 0.105658389;

  if ((decay_ptype == 211) || (decay_ptype == -211)) {
    return pimass;
  } else if ((decay_ptype == 321) || (decay_ptype == -321)) {
    return kmass;
  } else if ((decay_ptype == 311) || (decay_ptype == 310) ||
             (decay_ptype == 130)) {
    return k0mass;
  } else if ((decay_ptype == 13) || (decay_ptype == -13)) {
    return mumass;
  } else {
    std::cout << "[WARN] Unexpected neutrino deday-parent PDG code: "
              << decay_ptype << std::endl;
    return 0xdeadbeef;
  }
}

void DK2NuReader::GetEntry(UInt_t e) {
  CEnt = e;
  tree->GetEntry(CEnt);
}

void DK2NuReader::AddPPFXFriend(std::string treeName, std::string inputFiles,
                                size_t NUniverses) {

  ppfx_tree = OpenTChainWithFileList(treeName, inputFiles, ppfx_NFiles);

  if (NFiles != ppfx_NFiles) {
    std::cout << "[DK2NuReader]: Tried adding PPFX friend tree, but the DK2Nu "
                 "tree used a different number of files, unlikely that weights "
                 "will be in sync. Exiting hard for safety."
              << std::endl;
    exit(1);
  }

  ppfx_NEntries = ppfx_tree->GetEntries();
  ppfx_NUniverses = NUniverses;

  if (ppfx_NEntries != NEntries) {
    std::cout
        << "[DK2NuReader]: When adding PPFX friend tree, the DK2Nu tree "
           "contained a different number of entries than the PPFX friend, "
           "unlikely that weights will be in sync. Exiting hard for safety."
        << std::endl;
    exit(1);
  }

  tree->AddFriend(ppfx_tree);

  ppfx_vwgt_tot = new double[ppfx_NUniverses];

  tree->SetBranchAddress("cvwgt", &ppfx_cvwgt);
  tree->SetBranchAddress("vwgt_tot", ppfx_vwgt_tot);

  if (CheckTTreeHasBranch(ppfx_tree, "vwgt_mipp_pi")) {
    std::cout << "[INFO]: Input PPFX Friend tree has separate hadron "
                 "production weight branches."
              << std::endl;

    fUseAllPPFXBranches = true;

    ppfx_vwgt_mipp_pi = new double[NUniverses];
    ppfx_vwgt_mipp_K = new double[NUniverses];
    ppfx_vwgt_abs = new double[NUniverses];
    ppfx_vwgt_att = new double[NUniverses];
    ppfx_vwgt_ttpCpi = new double[NUniverses];
    ppfx_vwgt_ttpCk = new double[NUniverses];
    ppfx_vwgt_ttnCpi = new double[NUniverses];
    ppfx_vwgt_ttpCnu = new double[NUniverses];
    ppfx_vwgt_ttnua = new double[NUniverses];
    ppfx_vwgt_ttmesinc = new double[NUniverses];
    ppfx_vwgt_oth = new double[NUniverses];

    tree->SetBranchAddress("vwgt_mipp_pi", ppfx_vwgt_mipp_pi);
    tree->SetBranchAddress("vwgt_mipp_K", ppfx_vwgt_mipp_K);
    tree->SetBranchAddress("vwgt_abs", ppfx_vwgt_abs);
    tree->SetBranchAddress("vwgt_att", ppfx_vwgt_att);
    tree->SetBranchAddress("vwgt_ttpCpi", ppfx_vwgt_ttpCpi);
    tree->SetBranchAddress("vwgt_ttpCk", ppfx_vwgt_ttpCk);
    tree->SetBranchAddress("vwgt_ttnCpi", ppfx_vwgt_ttnCpi);
    tree->SetBranchAddress("vwgt_ttpCnu", ppfx_vwgt_ttpCnu);
    tree->SetBranchAddress("vwgt_ttnua", ppfx_vwgt_ttnua);
    tree->SetBranchAddress("vwgt_ttmesinc", ppfx_vwgt_ttmesinc);
    tree->SetBranchAddress("vwgt_oth", ppfx_vwgt_oth);

    tree->SetBranchAddress("cvwgt_mipp_pi", &ppfx_cvwgt_mipp_pi);
    tree->SetBranchAddress("cvwgt_mipp_K", &ppfx_cvwgt_mipp_K);
    tree->SetBranchAddress("cvwgt_abs", &ppfx_cvwgt_abs);
    tree->SetBranchAddress("cvwgt_att", &ppfx_cvwgt_att);
    tree->SetBranchAddress("cvwgt_ttpCpi", &ppfx_cvwgt_ttpCpi);
    tree->SetBranchAddress("cvwgt_ttpCk", &ppfx_cvwgt_ttpCk);
    tree->SetBranchAddress("cvwgt_ttnCpi", &ppfx_cvwgt_ttnCpi);
    tree->SetBranchAddress("cvwgt_ttpCnu", &ppfx_cvwgt_ttpCnu);
    tree->SetBranchAddress("cvwgt_ttnua", &ppfx_cvwgt_ttnua);
    tree->SetBranchAddress("cvwgt_ttmesinc", &ppfx_cvwgt_ttmesinc);
    tree->SetBranchAddress("cvwgt_oth", &ppfx_cvwgt_oth);
  }
}

UInt_t DK2NuReader::GetEntry() { return CEnt; }
UInt_t DK2NuReader::GetEntries() { return NEntries; }

void DK2NuReader::WriteOutLiteTree(TTree *outtree) {
  Int_t _decay_ntype;
  Double_t _decay_vx;
  Double_t _decay_vy;
  Double_t _decay_vz;
  Double_t _decay_pdpx;
  Double_t _decay_pdpy;
  Double_t _decay_pdpz;
  Double_t _decay_ppdxdz;
  Double_t _decay_ppdydz;
  Double_t _decay_pppz;
  Double_t _decay_ppenergy;
  Int_t _decay_ptype;
  Double_t _decay_muparpx;
  Double_t _decay_muparpy;
  Double_t _decay_muparpz;
  Double_t _decay_mupare;
  Double_t _decay_necm;
  Double_t _decay_nimpwt;

  outtree->Branch("decay_ntype", &_decay_ntype, "decay_ntype/I");
  outtree->Branch("decay_vx", &_decay_vx, "decay_vx/D");
  outtree->Branch("decay_vy", &_decay_vy, "decay_vy/D");
  outtree->Branch("decay_vz", &_decay_vz, "decay_vz/D");
  outtree->Branch("decay_pdpx", &_decay_pdpx, "decay_pdpx/D");
  outtree->Branch("decay_pdpy", &_decay_pdpy, "decay_pdpy/D");
  outtree->Branch("decay_pdpz", &_decay_pdpz, "decay_pdpz/D");
  outtree->Branch("decay_ppdxdz", &_decay_ppdxdz, "decay_ppdxdz/D");
  outtree->Branch("decay_ppdydz", &_decay_ppdydz, "decay_ppdydz/D");
  outtree->Branch("decay_pppz", &_decay_pppz, "decay_pppz/D");
  outtree->Branch("decay_ppenergy", &_decay_ppenergy, "decay_ppenergy/D");
  outtree->Branch("decay_ptype", &_decay_ptype, "decay_ptype/I");
  outtree->Branch("decay_muparpx", &_decay_muparpx, "decay_muparpx/D");
  outtree->Branch("decay_muparpy", &_decay_muparpy, "decay_muparpy/D");
  outtree->Branch("decay_muparpz", &_decay_muparpz, "decay_muparpz/D");
  outtree->Branch("decay_mupare", &_decay_mupare, "decay_mupare/D");
  outtree->Branch("decay_necm", &_decay_necm, "decay_necm/D");
  outtree->Branch("decay_nimpwt", &_decay_nimpwt, "decay_nimpwt/D");

  Double_t _ppfx_cvwgt;
  Double_t _ppfx_cvwgt_mipp_pi;
  Double_t _ppfx_cvwgt_mipp_K;
  Double_t _ppfx_cvwgt_abs;
  Double_t _ppfx_cvwgt_att;
  Double_t _ppfx_cvwgt_ttpCpi;
  Double_t _ppfx_cvwgt_ttpCk;
  Double_t _ppfx_cvwgt_ttnCpi;
  Double_t _ppfx_cvwgt_ttpCnu;
  Double_t _ppfx_cvwgt_ttnua;
  Double_t _ppfx_cvwgt_ttmesinc;
  Double_t _ppfx_cvwgt_oth;

  Double_t *_ppfx_vwgt_tot;
  TBranch *ppfx_vwgt_tot_branch;

  Double_t *_ppfx_vwgt_mipp_pi;
  TBranch *ppfx_vwgt_mipp_pi_branch;
  Double_t *_ppfx_vwgt_mipp_K;
  TBranch *ppfx_vwgt_mipp_K_branch;
  Double_t *_ppfx_vwgt_abs;
  TBranch *ppfx_vwgt_abs_branch;
  Double_t *_ppfx_vwgt_att;
  TBranch *ppfx_vwgt_att_branch;
  Double_t *_ppfx_vwgt_ttpCpi;
  TBranch *ppfx_vwgt_ttpCpi_branch;
  Double_t *_ppfx_vwgt_ttpCk;
  TBranch *ppfx_vwgt_ttpCk_branch;
  Double_t *_ppfx_vwgt_ttnCpi;
  TBranch *ppfx_vwgt_ttnCpi_branch;
  Double_t *_ppfx_vwgt_ttpCnu;
  TBranch *ppfx_vwgt_ttpCnu_branch;
  Double_t *_ppfx_vwgt_ttnua;
  TBranch *ppfx_vwgt_ttnua_branch;
  Double_t *_ppfx_vwgt_ttmesinc;
  TBranch *ppfx_vwgt_ttmesinc_branch;
  Double_t *_ppfx_vwgt_oth;
  TBranch *ppfx_vwgt_oth_branch;

  if (ppfx_tree) {
    outtree->Branch("ppfx_cvwgt", &_ppfx_cvwgt, "ppfx_cvwgt/D");

    _ppfx_vwgt_tot = new Double_t[ppfx_NUniverses];
    ppfx_vwgt_tot_branch =
        outtree->Branch("ppfx_vwgt_tot", _ppfx_vwgt_tot,
                        (std::string("ppfx_vwgt_tot[") +
                         std::to_string(ppfx_NUniverses) + "]/D")
                            .c_str());

    if (fUseAllPPFXBranches) {

      outtree->Branch("ppfx_cvwgt_mipp_pi", &_ppfx_cvwgt_mipp_pi,
                      "ppfx_cvwgt_mipp_pi/D");
      outtree->Branch("ppfx_cvwgt_mipp_K", &_ppfx_cvwgt_mipp_K,
                      "ppfx_cvwgt_mipp_K/D");
      outtree->Branch("ppfx_cvwgt_abs", &_ppfx_cvwgt_abs, "ppfx_cvwgt_abs/D");
      outtree->Branch("ppfx_cvwgt_att", &_ppfx_cvwgt_att, "ppfx_cvwgt_att/D");
      outtree->Branch("ppfx_cvwgt_ttpCpi", &_ppfx_cvwgt_ttpCpi,
                      "ppfx_cvwgt_ttpCpi/D");
      outtree->Branch("ppfx_cvwgt_ttpCk", &_ppfx_cvwgt_ttpCk,
                      "ppfx_cvwgt_ttpCk/D");
      outtree->Branch("ppfx_cvwgt_ttnCpi", &_ppfx_cvwgt_ttnCpi,
                      "ppfx_cvwgt_ttnCpi/D");
      outtree->Branch("ppfx_cvwgt_ttpCnu", &_ppfx_cvwgt_ttpCnu,
                      "ppfx_cvwgt_ttpCnu/D");
      outtree->Branch("ppfx_cvwgt_ttnua", &_ppfx_cvwgt_ttnua,
                      "ppfx_cvwgt_ttnua/D");
      outtree->Branch("ppfx_cvwgt_ttmesinc", &_ppfx_cvwgt_ttmesinc,
                      "ppfx_cvwgt_ttmesinc/D");
      outtree->Branch("ppfx_cvwgt_oth", &_ppfx_cvwgt_oth, "ppfx_cvwgt_oth/D");

      _ppfx_vwgt_mipp_pi = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_mipp_pi_branch =
          outtree->Branch("ppfx_vwgt_mipp_pi", _ppfx_vwgt_mipp_pi,
                          (std::string("ppfx_vwgt_mipp_pi[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_mipp_K = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_mipp_K_branch =
          outtree->Branch("ppfx_vwgt_mipp_K", _ppfx_vwgt_mipp_K,
                          (std::string("ppfx_vwgt_mipp_K[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_abs = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_abs_branch =
          outtree->Branch("ppfx_vwgt_abs", _ppfx_vwgt_abs,
                          (std::string("ppfx_vwgt_abs[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_att = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_att_branch =
          outtree->Branch("ppfx_vwgt_att", _ppfx_vwgt_att,
                          (std::string("ppfx_vwgt_att[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttpCpi = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttpCpi_branch =
          outtree->Branch("ppfx_vwgt_ttpCpi", _ppfx_vwgt_ttpCpi,
                          (std::string("ppfx_vwgt_ttpCpi[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttpCk = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttpCk_branch =
          outtree->Branch("ppfx_vwgt_ttpCk", _ppfx_vwgt_ttpCk,
                          (std::string("ppfx_vwgt_ttpCk[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttnCpi = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttnCpi_branch =
          outtree->Branch("ppfx_vwgt_ttnCpi", _ppfx_vwgt_ttnCpi,
                          (std::string("ppfx_vwgt_ttnCpi[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttpCnu = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttpCnu_branch =
          outtree->Branch("ppfx_vwgt_ttpCnu", _ppfx_vwgt_ttpCnu,
                          (std::string("ppfx_vwgt_ttpCnu[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttnua = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttnua_branch =
          outtree->Branch("ppfx_vwgt_ttnua", _ppfx_vwgt_ttnua,
                          (std::string("ppfx_vwgt_ttnua[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_ttmesinc = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_ttmesinc_branch =
          outtree->Branch("ppfx_vwgt_ttmesinc", _ppfx_vwgt_ttmesinc,
                          (std::string("ppfx_vwgt_ttmesinc[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());

      _ppfx_vwgt_oth = new Double_t[ppfx_NUniverses];
      ppfx_vwgt_oth_branch =
          outtree->Branch("ppfx_vwgt_oth", _ppfx_vwgt_oth,
                          (std::string("ppfx_vwgt_oth[") +
                           std::to_string(ppfx_NUniverses) + "]/D")
                              .c_str());
    }
  }

  for (UInt_t ent = 0; ent < NEntries; ++ent) {
    tree->GetEntry(ent);

    _decay_ntype = decay_ntype;
    _decay_vx = decay_vx;
    _decay_vy = decay_vy;
    _decay_vz = decay_vz;
    _decay_pdpx = decay_pdpx;
    _decay_pdpy = decay_pdpy;
    _decay_pdpz = decay_pdpz;
    _decay_ppdxdz = decay_ppdxdz;
    _decay_ppdydz = decay_ppdydz;
    _decay_pppz = decay_pppz;
    _decay_ppenergy = decay_ppenergy;
    _decay_ptype = decay_ptype;
    _decay_muparpx = decay_muparpx;
    _decay_muparpy = decay_muparpy;
    _decay_muparpz = decay_muparpz;
    _decay_mupare = decay_mupare;
    _decay_necm = decay_necm;
    _decay_nimpwt = decay_nimpwt;

    if (ppfx_tree) {
      _ppfx_cvwgt = ppfx_cvwgt;
      std::copy_n(ppfx_vwgt_tot, ppfx_NUniverses, _ppfx_vwgt_tot);
      if (fUseAllPPFXBranches) {
        std::copy_n(ppfx_vwgt_mipp_pi, ppfx_NUniverses, _ppfx_vwgt_mipp_pi);
        std::copy_n(ppfx_vwgt_mipp_K, ppfx_NUniverses, _ppfx_vwgt_mipp_K);
        std::copy_n(ppfx_vwgt_abs, ppfx_NUniverses, _ppfx_vwgt_abs);
        std::copy_n(ppfx_vwgt_att, ppfx_NUniverses, _ppfx_vwgt_att);
        std::copy_n(ppfx_vwgt_ttpCpi, ppfx_NUniverses, _ppfx_vwgt_ttpCpi);
        std::copy_n(ppfx_vwgt_ttpCk, ppfx_NUniverses, _ppfx_vwgt_ttpCk);
        std::copy_n(ppfx_vwgt_ttnCpi, ppfx_NUniverses, _ppfx_vwgt_ttnCpi);
        std::copy_n(ppfx_vwgt_ttpCnu, ppfx_NUniverses, _ppfx_vwgt_ttpCnu);
        std::copy_n(ppfx_vwgt_ttnua, ppfx_NUniverses, _ppfx_vwgt_ttnua);
        std::copy_n(ppfx_vwgt_ttmesinc, ppfx_NUniverses, _ppfx_vwgt_ttmesinc);
        std::copy_n(ppfx_vwgt_oth, ppfx_NUniverses, _ppfx_vwgt_oth);

        _ppfx_cvwgt_mipp_pi = ppfx_cvwgt_mipp_pi;
        _ppfx_cvwgt_mipp_K = ppfx_cvwgt_mipp_K;
        _ppfx_cvwgt_abs = ppfx_cvwgt_abs;
        _ppfx_cvwgt_att = ppfx_cvwgt_att;
        _ppfx_cvwgt_ttpCpi = ppfx_cvwgt_ttpCpi;
        _ppfx_cvwgt_ttpCk = ppfx_cvwgt_ttpCk;
        _ppfx_cvwgt_ttnCpi = ppfx_cvwgt_ttnCpi;
        _ppfx_cvwgt_ttpCnu = ppfx_cvwgt_ttpCnu;
        _ppfx_cvwgt_ttnua = ppfx_cvwgt_ttnua;
        _ppfx_cvwgt_ttmesinc = ppfx_cvwgt_ttmesinc;
        _ppfx_cvwgt_oth = ppfx_cvwgt_oth;
      }
    }

    outtree->Fill();
  }

  if (ppfx_tree) {
    outtree->ResetBranchAddress(ppfx_vwgt_tot_branch);
    delete[] _ppfx_vwgt_tot;
    if (fUseAllPPFXBranches) {
      outtree->ResetBranchAddress(ppfx_vwgt_mipp_pi_branch);
      delete[] _ppfx_vwgt_mipp_pi;
      outtree->ResetBranchAddress(ppfx_vwgt_mipp_K_branch);
      delete[] _ppfx_vwgt_mipp_K;
      outtree->ResetBranchAddress(ppfx_vwgt_abs_branch);
      delete[] _ppfx_vwgt_abs;
      outtree->ResetBranchAddress(ppfx_vwgt_att_branch);
      delete[] _ppfx_vwgt_att;
      outtree->ResetBranchAddress(ppfx_vwgt_ttpCpi_branch);
      delete[] _ppfx_vwgt_ttpCpi;
      outtree->ResetBranchAddress(ppfx_vwgt_ttpCk_branch);
      delete[] _ppfx_vwgt_ttpCk;
      outtree->ResetBranchAddress(ppfx_vwgt_ttnCpi_branch);
      delete[] _ppfx_vwgt_ttnCpi;
      outtree->ResetBranchAddress(ppfx_vwgt_ttpCnu_branch);
      delete[] _ppfx_vwgt_ttpCnu;
      outtree->ResetBranchAddress(ppfx_vwgt_ttnua_branch);
      delete[] _ppfx_vwgt_ttnua;
      outtree->ResetBranchAddress(ppfx_vwgt_ttmesinc_branch);
      delete[] _ppfx_vwgt_ttmesinc;
      outtree->ResetBranchAddress(ppfx_vwgt_oth_branch);
      delete[] _ppfx_vwgt_oth;
    }
  }

  tree->GetEntry(CEnt);
}

DK2NuReader::~DK2NuReader() {
  delete tree;
  if (ppfx_tree) {
    delete ppfx_tree;
    delete[] ppfx_vwgt_tot;
  }
}

void DK2NuReader::DumpMaxWeights() {
  std::cout << "\tppfx_maxweight: " << ppfx_maxweight
            << "\n\tppfx_abs_maxweight: " << ppfx_abs_maxweight
            << "\n\tppfx_att_maxweight: " << ppfx_att_maxweight
            << "\n\tppfx_ttpCpi_maxweight: " << ppfx_ttpCpi_maxweight
            << "\n\tppfx_ttpCk_maxweight: " << ppfx_ttpCk_maxweight
            << "\n\tppfx_ttnCpi_maxweight: " << ppfx_ttnCpi_maxweight
            << "\n\tppfx_ttpCnu_maxweight: " << ppfx_ttpCnu_maxweight
            << "\n\tppfx_ttnua_maxweight: " << ppfx_ttnua_maxweight
            << "\n\tppfx_ttmesinc_maxweight: " << ppfx_ttmesinc_maxweight
            << "\n\tppfx_oth_maxweight: " << ppfx_oth_maxweight << std::endl;
}

DKMetaReader::DKMetaReader(std::string treeName, std::string inputFiles,
                           bool DK2NULite)
/*: beamsim(0),
  physics(0),
  physcuts(0),
  tgtcfg(0),
  horncfg(0),
  dkvolcfg(0) */
{
  tree = OpenTChainWithFileList(treeName, inputFiles, NFiles);
  NEntries = tree->GetEntries();

  SetBranchAddresses(DK2NULite);
  std::cout << "[DKMetaReader]: Loaded TChain: " << NFiles << " files and "
            << NEntries << " entries." << std::endl;
  GetEntry(0);
}

void DKMetaReader::SetBranchAddresses(bool DK2NULite) {
  if (!DK2NULite) {
    tree->SetMakeClass(true);

    tree->SetBranchAddress("job", &job);
    tree->SetBranchAddress("pots", &pots);
    // tree->SetBranchAddress("beamsim", &beamsim);
    // tree->SetBranchAddress("physics", &physics);
    // tree->SetBranchAddress("physcuts", &physcuts);
    // tree->SetBranchAddress("tgtcfg", &tgtcfg);
    // tree->SetBranchAddress("horncfg", &horncfg);
    // tree->SetBranchAddress("dkvolcfg", &dkvolcfg);
    tree->SetBranchAddress("beam0x", &beam0x);
    tree->SetBranchAddress("beam0y", &beam0y);
    tree->SetBranchAddress("beam0z", &beam0z);
    tree->SetBranchAddress("beamhwidth", &beamhwidth);
    tree->SetBranchAddress("beamvwidth", &beamvwidth);
    tree->SetBranchAddress("beamdxdz", &beamdxdz);
    tree->SetBranchAddress("beamdydz", &beamdydz);
    tree->SetBranchAddress("location", &location_);
    tree->SetBranchAddress("location.x", location_x);
    tree->SetBranchAddress("location.y", location_y);
    tree->SetBranchAddress("location.z", location_z);
    // tree->SetBranchAddress("location.name", location_name);
    // tree->SetBranchAddress("vintnames", &vintnames);
    // tree->SetBranchAddress("vdblnames", &vdblnames);
  } else {
    tree->SetBranchAddress("pots", &pots);
    // tree->SetBranchAddress("beamsim", &beamsim);
    // tree->SetBranchAddress("physics", &physics);
    // tree->SetBranchAddress("physcuts", &physcuts);
    // tree->SetBranchAddress("tgtcfg", &tgtcfg);
    // tree->SetBranchAddress("horncfg", &horncfg);
    // tree->SetBranchAddress("dkvolcfg", &dkvolcfg);
    tree->SetBranchAddress("beam0x", &beam0x);
    tree->SetBranchAddress("beam0y", &beam0y);
    tree->SetBranchAddress("beam0z", &beam0z);
    tree->SetBranchAddress("beamhwidth", &beamhwidth);
    tree->SetBranchAddress("beamvwidth", &beamvwidth);
  }
}

void DKMetaReader::GetEntry(UInt_t e) {
  CEnt = e;
  tree->GetEntry(CEnt);
}
UInt_t DKMetaReader::GetEntry() { return CEnt; }
UInt_t DKMetaReader::GetEntries() { return NEntries; }

void DKMetaReader::WriteOutLiteTree(TTree *outtree) {
  double _pots;
  // std::string *_beamsim;
  // std::string *_physics;
  // std::string *_physcuts;
  // std::string *_tgtcfg;
  // std::string *_horncfg;
  // std::string *_dkvolcfg;
  double _beam0x;
  double _beam0y;
  double _beam0z;
  double _beamhwidth;
  double _beamvwidth;

  outtree->Branch("pots", &_pots);
  // outtree->Branch("beamsim", &_beamsim);
  // outtree->Branch("physics", &_physics);
  // outtree->Branch("physcuts", &_physcuts);
  // outtree->Branch("tgtcfg", &_tgtcfg);
  // outtree->Branch("horncfg", &_horncfg);
  // outtree->Branch("dkvolcfg", &_dkvolcfg);
  outtree->Branch("beam0x", &_beam0x);
  outtree->Branch("beam0y", &_beam0y);
  outtree->Branch("beam0z", &_beam0z);
  outtree->Branch("beamhwidth", &_beamhwidth);
  outtree->Branch("beamvwidth", &_beamvwidth);

  for (UInt_t ent = 0; ent < NEntries; ++ent) {
    tree->GetEntry(ent);

    _pots = pots;
    // (*_beamsim) = (*beamsim);
    // (*_physics) = (*physics);
    // (*_physcuts) = (*physcuts);
    // (*_tgtcfg) = (*tgtcfg);
    // (*_horncfg) = (*horncfg);
    // (*_dkvolcfg) = (*dkvolcfg);
    _beam0x = beam0x;
    _beam0y = beam0y;
    _beam0z = beam0z;
    _beamhwidth = beamhwidth;
    _beamvwidth = beamvwidth;

    outtree->Fill();
  }

  tree->GetEntry(CEnt);
}

DKMetaReader::~DKMetaReader() { delete tree; }

std::string GetPPFXHistName(Int_t PPFXUniv, Int_t NPPFXUniv) {
  if (PPFXUniv == 0) {
    return "_Nom";

  } else if (PPFXUniv == 1) {
    return "_CV";
  }

  PPFXUniv -= 2;

  int univ_cycle = PPFXUniv / NPPFXUniv;
  int univ = PPFXUniv % NPPFXUniv;

  if (univ_cycle > DK2NuReader::kNPPFXAllWeights) {
    std::cout << "[ERROR]: Tried to build histogram name for NPPFXU = "
              << PPFXUniv << " but only have " << NPPFXUniv << " = cycle "
              << univ_cycle << std::endl;
    throw;
  }

  std::string const cycle_names[] = {
      // "",       "_mipp_pi", "_mipp_K", "_abs",   "_att",      "_ttpCpi",
      "",        "_abs",    "_att",   "_ttpCpi",   "_ttpCk",
      "_ttnCpi", "_ttpCnu", "_ttnua", "_ttmesinc", "_oth"};

  std::stringstream ss("");
  ss << "_univ_" << univ << cycle_names[univ_cycle];

  return ss.str();
}

double GetPPFXWeight(Int_t PPFXUniv, Int_t NPPFXUniv, DK2NuReader &dk2nuRdr) {
  if (PPFXUniv == 0) {
    return 1;
  } else if (PPFXUniv == 1) {
    dk2nuRdr.ppfx_maxweight =
        std::max(dk2nuRdr.ppfx_maxweight, dk2nuRdr.ppfx_cvwgt);
    return dk2nuRdr.ppfx_cvwgt;
  }

  PPFXUniv -= 2;

  int univ_cycle = PPFXUniv / NPPFXUniv;
  int univ = PPFXUniv % NPPFXUniv;

  double W = 0;

  switch (univ_cycle) {
  case 0: {
    W = dk2nuRdr.ppfx_vwgt_tot[univ];
    dk2nuRdr.ppfx_maxweight = std::max(dk2nuRdr.ppfx_maxweight, W);
    break;
  }
  case 1: {
    W = dk2nuRdr.ppfx_vwgt_abs[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_abs);
    dk2nuRdr.ppfx_abs_maxweight = std::max(dk2nuRdr.ppfx_abs_maxweight, W);
    break;
  }
  case 2: {
    W = dk2nuRdr.ppfx_vwgt_att[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_att);
    dk2nuRdr.ppfx_att_maxweight = std::max(dk2nuRdr.ppfx_att_maxweight, W);
    break;
  }
  case 3: {
    W = dk2nuRdr.ppfx_vwgt_ttpCpi[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttpCpi);
    dk2nuRdr.ppfx_ttpCpi_maxweight =
        std::max(dk2nuRdr.ppfx_ttpCpi_maxweight, W);
    break;
  }
  case 4: {
    W = dk2nuRdr.ppfx_vwgt_ttpCk[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttpCk);
    dk2nuRdr.ppfx_ttpCk_maxweight = std::max(dk2nuRdr.ppfx_ttpCk_maxweight, W);
    break;
  }
  case 5: {
    W = dk2nuRdr.ppfx_vwgt_ttnCpi[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttnCpi);
    dk2nuRdr.ppfx_ttnCpi_maxweight =
        std::max(dk2nuRdr.ppfx_ttnCpi_maxweight, W);
    break;
  }
  case 6: {
    W = dk2nuRdr.ppfx_vwgt_ttpCnu[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttpCnu);
    dk2nuRdr.ppfx_ttpCnu_maxweight =
        std::max(dk2nuRdr.ppfx_ttpCnu_maxweight, W);
    break;
  }
  case 7: {
    W = dk2nuRdr.ppfx_vwgt_ttnua[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttnua);
    dk2nuRdr.ppfx_ttnua_maxweight = std::max(dk2nuRdr.ppfx_ttnua_maxweight, W);
    break;
  }
  case 8: {
    W = dk2nuRdr.ppfx_vwgt_ttmesinc[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_ttmesinc);
    dk2nuRdr.ppfx_ttmesinc_maxweight =
        std::max(dk2nuRdr.ppfx_ttmesinc_maxweight, W);
    break;
  }
  case 9: {
    W = dk2nuRdr.ppfx_vwgt_oth[univ] *
        (dk2nuRdr.ppfx_cvwgt / dk2nuRdr.ppfx_cvwgt_oth);
    dk2nuRdr.ppfx_oth_maxweight = std::max(dk2nuRdr.ppfx_oth_maxweight, W);
    break;
  }
  default: { throw; }
  }

  return std::isnormal(W) ? W : 0;
}

void DumpPPFXWeights(Int_t PPFXUniv, DK2NuReader &dk2nuRdr) {
  std::cout << "[INFO]: cvwgt[" << PPFXUniv << "] = " << dk2nuRdr.ppfx_cvwgt
            << std::endl;
  std::cout << "\tppfx_vwgt_tot[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_tot[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_mipp_pi[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_mipp_pi[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_mipp_K[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_mipp_K[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_abs[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_abs[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_att[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_att[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttpCpi[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttpCpi[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttpCk[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttpCk[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttnCpi[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttnCpi[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttpCnu[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttpCnu[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttnua[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttnua[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_ttmesinc[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_ttmesinc[PPFXUniv] << std::endl;
  std::cout << "\tppfx_vwgt_oth[" << PPFXUniv
            << "] = " << dk2nuRdr.ppfx_vwgt_oth[PPFXUniv] << std::endl;
}
