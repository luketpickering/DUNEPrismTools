#include "TChain.h"
#include "TObjString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>

struct FullDetTreeReader {
  FullDetTreeReader(std::string treeName, std::string inputFiles) {
    tree = new TChain(treeName.c_str());

    NFiles = tree->Add(inputFiles.c_str());
    NEntries = tree->GetEntries();

    SetBranchAddresses();
    std::cout << "[FullDetTreeReader]: Loaded TChain: " << NFiles
              << " files and " << NEntries << " entries." << std::endl;
    GetEntry(0);
  }

  static const Int_t kNXBins = 400;
  static const Int_t kNYBins = 3;
  static const Int_t kNZBins = 3;

  Double_t Enu;
  Int_t nuPDG;
  Double_t vtx_X;
  Double_t vtx_Y;
  Double_t vtx_Z;
  Int_t eventNum;
  TObjString * eventCode;
  Int_t lepPDG;
  Int_t nLep;
  Double_t Q2True;
  TLorentzVector * qTrue;
  Double_t yTrue;
  Double_t W_rest;
  Int_t nPi0;
  Int_t nPiC;
  Int_t nProton;
  Int_t nNeutron;
  Int_t nGamma;

  Double_t eLepTrue;
  Double_t ePi0True;
  Double_t ePiCTrue;
  Double_t eProtonTrue;
  Double_t eNeutronTrue;
  Double_t eGammaTrue;
  Double_t eHadTrueTotal;

  Double_t eLepPrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eHadPrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eProtonPrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eNeutronPrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t ePiCPrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t ePi0PrimaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eOtherPrimaryDep[kNXBins][kNYBins][kNZBins];

  Double_t eLepSecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eHadSecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eProtonSecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eNeutronSecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t ePiCSecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t ePi0SecondaryDep[kNXBins][kNYBins][kNZBins];
  Double_t eOtherSecondaryDep[kNXBins][kNYBins][kNZBins];

  Double_t eElectronShowerDepInside[kNXBins];
  Double_t eElectronShowerDepOutside;

  Double_t lepTrackX[1000];
  Double_t lepTrackY[1000];
  Double_t lepTrackZ[1000];
  Double_t lepTrackMomX[1000];
  Double_t lepTrackMomY[1000];
  Double_t lepTrackMomZ[1000];

  bool flagLepExit;
  bool flagLepExitBack;
  bool flagLepExitFront;
  bool flagLepExitY;

  double lepExitingPosX;
  double lepExitingPosY;
  double lepExitingPosZ;
  double lepExitingMomX;
  double lepExitingMomY;
  double lepExitingMomZ;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void SetBranchAddresses() {
    tree->SetBranchAddress("Enu", &Enu);
    tree->SetBranchAddress("nuPDG", &nuPDG);
    tree->SetBranchAddress("vtx_X", &vtx_X);
    tree->SetBranchAddress("vtx_Y", &vtx_Y);
    tree->SetBranchAddress("vtx_Z", &vtx_Z);
    tree->SetBranchAddress("eventNum", &eventNum);
    tree->SetBranchAddress("eventCode", &eventCode);
    tree->SetBranchAddress("lepPDG", &lepPDG);
    tree->SetBranchAddress("nLep", &nLep);
    tree->SetBranchAddress("Q2True", &Q2True);
    tree->SetBranchAddress("qTrue", &qTrue);
    tree->SetBranchAddress("yTrue", &yTrue);
    tree->SetBranchAddress("W_rest", &W_rest);
    tree->SetBranchAddress("nPi0", &nPi0);
    tree->SetBranchAddress("nPiC", &nPiC);
    tree->SetBranchAddress("nProton", &nProton);
    tree->SetBranchAddress("nNeutron", &nNeutron);
    tree->SetBranchAddress("nGamma", &nGamma);

    tree->SetBranchAddress("eLepTrue", & eLepTrue);
    tree->SetBranchAddress("eHadTrueTotal", & eHadTrueTotal);
    tree->SetBranchAddress("ePi0True", & ePi0True);
    tree->SetBranchAddress("ePiCTrue", & ePiCTrue);
    tree->SetBranchAddress("eProtonTrue", & eProtonTrue);
    tree->SetBranchAddress("eNeutronTrue", & eNeutronTrue);
    tree->SetBranchAddress("eGammaTrue", & eGammaTrue);

    tree->SetBranchAddress("eLepPrimaryDepFull", &eLepPrimaryDep);
    tree->SetBranchAddress("eHadPrimaryDep", &eHadPrimaryDep);
    tree->SetBranchAddress("eProtonPrimaryDep", &eProtonPrimaryDep);
    tree->SetBranchAddress("eNeutronPrimaryDep", &eNeutronPrimaryDep);
    tree->SetBranchAddress("ePiCPrimaryDep", &ePiCPrimaryDep);
    tree->SetBranchAddress("ePi0PrimaryDep", &ePi0PrimaryDep);
    tree->SetBranchAddress("eOtherPrimaryDep", &eOtherPrimaryDep);

    tree->SetBranchAddress("eLepSecondaryDepFull", &eLepSecondaryDep);
    tree->SetBranchAddress("eHadSecondaryDep", &eHadSecondaryDep);
    tree->SetBranchAddress("eProtonSecondaryDep", &eProtonSecondaryDep);
    tree->SetBranchAddress("eNeutronSecondaryDep", &eNeutronSecondaryDep);
    tree->SetBranchAddress("ePiCSecondaryDep", &ePiCSecondaryDep);
    tree->SetBranchAddress("ePi0SecondaryDep", &ePi0SecondaryDep);
    tree->SetBranchAddress("eOtherSecondaryDep", &eOtherSecondaryDep);

    tree->SetBranchAddress("flagLepExitBack", &flagLepExitBack);
    tree->SetBranchAddress("flagLepExitFront", &flagLepExitFront);
    tree->SetBranchAddress("flagLepExitY", &flagLepExitY);
    tree->SetBranchAddress("flagLepExit", &flagLepExit);
    tree->SetBranchAddress("lepExitingPosX",&lepExitingPosX);
    tree->SetBranchAddress("lepExitingPosY",&lepExitingPosY);
    tree->SetBranchAddress("lepExitingPosZ",&lepExitingPosZ); 
    tree->SetBranchAddress("lepExitingMomX",&lepExitingMomX);
    tree->SetBranchAddress("lepExitingMomY",&lepExitingMomY);
    tree->SetBranchAddress("lepExitingMomZ",&lepExitingMomZ);
    tree->SetBranchAddress("lepTrackX", &lepTrackX);
    tree->SetBranchAddress("lepTrackY", &lepTrackY);
    tree->SetBranchAddress("lepTrackZ", &lepTrackZ);
    tree->SetBranchAddress("lepTrackMomX", &lepTrackMomX);
    tree->SetBranchAddress("lepTrackMomY", &lepTrackMomY);
    tree->SetBranchAddress("lepTrackMomZ", &lepTrackMomZ);
    tree->SetBranchAddress("eElectronShowerDepInside", &eElectronShowerDepInside);
    tree->SetBranchAddress("eElectronShowerDepOutside", &eElectronShowerDepOutside);
  }

  void GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t GetEntry() { return CEnt; }
  UInt_t GetEntries() { return NEntries; }

  ~FullDetTreeReader() { delete tree; };
};
