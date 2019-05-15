#ifndef DK2NU_TREEREADER_HXX_SEEN
#define DK2NU_TREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>
#include <vector>

struct DK2NuReader {
  DK2NuReader(std::string treeName, std::string inputFiles, bool DK2NULite);

  static const Int_t kMaxnuray = 3;
  static const Int_t kMaxancestor = 13;
  static const Int_t kMaxtraj = 10;
  static const Int_t kNPPFXAllWeights = 10;

  Int_t job;
  Int_t potnum;
  Int_t decay_norig;
  Int_t decay_ndecay;
  Int_t decay_ntype;
  Double_t decay_vx;
  Double_t decay_vy;
  Double_t decay_vz;
  Double_t decay_pdpx;
  Double_t decay_pdpy;
  Double_t decay_pdpz;
  Double_t decay_ppdxdz;
  Double_t decay_ppdydz;
  Double_t decay_pppz;
  Double_t decay_ppenergy;
  Int_t decay_ppmedium;
  Int_t decay_ptype;
  Double_t decay_muparpx;
  Double_t decay_muparpy;
  Double_t decay_muparpz;
  Double_t decay_mupare;
  Double_t decay_necm;
  Double_t decay_nimpwt;
  Int_t nuray_;
  Double_t nuray_px[kMaxnuray];  //[nuray_]
  Double_t nuray_py[kMaxnuray];  //[nuray_]
  Double_t nuray_pz[kMaxnuray];  //[nuray_]
  Double_t nuray_E[kMaxnuray];   //[nuray_]
  Double_t nuray_wgt[kMaxnuray]; //[nuray_]
  Int_t ancestor_;
  Int_t ancestor_pdg[kMaxancestor];        //[ancestor_]
  Double_t ancestor_startx[kMaxancestor];  //[ancestor_]
  Double_t ancestor_starty[kMaxancestor];  //[ancestor_]
  Double_t ancestor_startz[kMaxancestor];  //[ancestor_]
  Double_t ancestor_startt[kMaxancestor];  //[ancestor_]
  Double_t ancestor_startpx[kMaxancestor]; //[ancestor_]
  Double_t ancestor_startpy[kMaxancestor]; //[ancestor_]
  Double_t ancestor_startpz[kMaxancestor]; //[ancestor_]
  Double_t ancestor_stoppx[kMaxancestor];  //[ancestor_]
  Double_t ancestor_stoppy[kMaxancestor];  //[ancestor_]
  Double_t ancestor_stoppz[kMaxancestor];  //[ancestor_]
  Double_t ancestor_polx[kMaxancestor];    //[ancestor_]
  Double_t ancestor_poly[kMaxancestor];    //[ancestor_]
  Double_t ancestor_polz[kMaxancestor];    //[ancestor_]
  Double_t ancestor_pprodpx[kMaxancestor]; //[ancestor_]
  Double_t ancestor_pprodpy[kMaxancestor]; //[ancestor_]
  Double_t ancestor_pprodpz[kMaxancestor]; //[ancestor_]
  Int_t ancestor_nucleus[kMaxancestor];    //[ancestor_]
  std::string ancestor_proc[kMaxancestor];
  std::string ancestor_ivol[kMaxancestor];
  std::string ancestor_imat[kMaxancestor];
  Double_t ppvx;
  Double_t ppvy;
  Double_t ppvz;
  Double_t tgtexit_tvx;
  Double_t tgtexit_tvy;
  Double_t tgtexit_tvz;
  Double_t tgtexit_tpx;
  Double_t tgtexit_tpy;
  Double_t tgtexit_tpz;
  Int_t tgtexit_tptype;
  Int_t tgtexit_tgen;
  Int_t traj_;
  Double_t traj_trkx[kMaxtraj];  //[traj_]
  Double_t traj_trky[kMaxtraj];  //[traj_]
  Double_t traj_trkz[kMaxtraj];  //[traj_]
  Double_t traj_trkpx[kMaxtraj]; //[traj_]
  Double_t traj_trkpy[kMaxtraj]; //[traj_]
  Double_t traj_trkpz[kMaxtraj]; //[traj_]
  Int_t flagbits;
  std::vector<Int_t> vint;
  std::vector<Double_t> vdbl;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void SetBranchAddresses(bool DK2NULite = false);
  void SetPPFXBranchAddresses(size_t NUniverses = 100);

  TChain *ppfx_tree;
  UInt_t ppfx_NFiles;
  UInt_t ppfx_NEntries;
  double ppfx_cvwgt;
  double *ppfx_vwgt_tot;

  bool fUseAllPPFXBranches;
  double *ppfx_vwgt_mipp_pi;
  double *ppfx_vwgt_mipp_K;
  double *ppfx_vwgt_abs;
  double *ppfx_vwgt_att;
  double *ppfx_vwgt_ttpCpi;
  double *ppfx_vwgt_ttpCk;
  double *ppfx_vwgt_ttnCpi;
  double *ppfx_vwgt_ttpCnu;
  double *ppfx_vwgt_ttnua;
  double *ppfx_vwgt_ttmesinc;
  double *ppfx_vwgt_oth;

  double ppfx_cvwgt_mipp_pi;
  double ppfx_cvwgt_mipp_K;
  double ppfx_cvwgt_abs;
  double ppfx_cvwgt_att;
  double ppfx_cvwgt_ttpCpi;
  double ppfx_cvwgt_ttpCk;
  double ppfx_cvwgt_ttnCpi;
  double ppfx_cvwgt_ttpCnu;
  double ppfx_cvwgt_ttnua;
  double ppfx_cvwgt_ttmesinc;
  double ppfx_cvwgt_oth;

  UInt_t ppfx_NUniverses;

  void AddPPFXFriend(std::string treeName, std::string inputFiles,
                     size_t NUniverses);

  double GetParentMass();

  void GetEntry(UInt_t e);
  UInt_t GetEntry();
  UInt_t GetEntries();

  void WriteOutLiteTree(TTree *outtree);

  ~DK2NuReader();
};

struct DKMetaReader {
  DKMetaReader(std::string treeName, std::string inputFiles, bool DK2NULite);

  static const Int_t kMaxlocation = 3;

  Int_t job;
  Double_t pots;
  // std::string *beamsim;
  // std::string *physics;
  // std::string *physcuts;
  // std::string *tgtcfg;
  // std::string *horncfg;
  // std::string *dkvolcfg;
  Double_t beam0x;
  Double_t beam0y;
  Double_t beam0z;
  Double_t beamhwidth;
  Double_t beamvwidth;
  Double_t beamdxdz;
  Double_t beamdydz;
  Int_t location_;
  Double_t location_x[kMaxlocation]; //[location_]
  Double_t location_y[kMaxlocation]; //[location_]
  Double_t location_z[kMaxlocation]; //[location_]
  // std::string location_name[kMaxlocation];
  // std::vector<std::string> vintnames;
  // std::vector<std::string> vdblnames;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void SetBranchAddresses(bool DK2NULite);

  void GetEntry(UInt_t e);
  UInt_t GetEntry();
  UInt_t GetEntries();

  void WriteOutLiteTree(TTree *outtree);

  ~DKMetaReader();
};

std::string GetPPFXHistName(Int_t PPFXUniv, Int_t NPPFXUniv);

double GetPPFXWeight(Int_t PPFXUniv, Int_t NPPFXUniv, DK2NuReader &dk2nuRdr);
void DumpPPFXWeights(Int_t PPFXUniv, DK2NuReader &dk2nuRdr);

#endif
