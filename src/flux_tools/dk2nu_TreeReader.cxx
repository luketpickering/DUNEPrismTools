#include "dk2nu_TreeReader.hxx"

#include <vector>
#include <iostream>

  DK2NuReader::DK2NuReader(std::string treeName, std::string inputFiles, bool DK2NULite) {
    tree = new TChain(treeName.c_str());
    NFiles = tree->Add(inputFiles.c_str());
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

  double DK2NuReader::GetParentMass() {
    static double const pimass = 0.13957;  // in GeV
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
      std::cout << "[WARN] Unexpected neutrino deday-parent PDG code: " << decay_ptype << std::endl;
      return 0xdeadbeef;
    }
  }

  void DK2NuReader::GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
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

      outtree->Fill();
    }

    tree->GetEntry(CEnt);
  }

  DK2NuReader::~DK2NuReader() { delete tree; }


  DKMetaReader::DKMetaReader(std::string treeName, std::string inputFiles, bool DK2NULite)
      /*: beamsim(0),
        physics(0),
        physcuts(0),
        tgtcfg(0),
        horncfg(0),
        dkvolcfg(0) */{
    tree = new TChain(treeName.c_str());
    NFiles = tree->Add(inputFiles.c_str());
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
