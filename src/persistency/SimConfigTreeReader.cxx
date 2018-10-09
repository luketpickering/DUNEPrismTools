#include "SimConfigTreeReader.hxx"

#include "ROOTUtility.hxx"

#include <algorithm>

SimConfig::SimConfig(std::string const &inputFile) {

    LoadTree(inputFile);
  }

  std::string SimConfig::TreeName(){
    return "SimConfigTree";
  }

  void SimConfig::Reset() {
    NXSteps = 0;
    std::fill_n(DetMin,3,0);
    std::fill_n(DetMax,3,0);
    std::fill_n(VetoGap,3,0);
    NMaxTrackSteps = 0;
    POTPerFile = 0xdeadbeef;
    timesep_us = 0xdeadbeef;
  }

  void SimConfig::Copy(SimConfig const &other) {
    NXSteps = other.NXSteps;
    std::copy_n(other.DetMin,3,DetMin);
    std::copy_n(other.DetMax,3,DetMax);
    std::copy_n(other.VetoGap,3,VetoGap);
    NMaxTrackSteps = other.NMaxTrackSteps;
    POTPerFile = other.POTPerFile;
    timesep_us = other.timesep_us;
  }

  void SimConfig::SetBranchAddresses() {
    tree->SetBranchAddress("NXSteps", &NXSteps);
    tree->SetBranchAddress("DetMin", &DetMin);
    tree->SetBranchAddress("DetMax", &DetMax);
    tree->SetBranchAddress("VetoGap", &VetoGap);
    tree->SetBranchAddress("NMaxTrackSteps", &NMaxTrackSteps);
    tree->SetBranchAddress("POTPerFile", &POTPerFile);
    tree->SetBranchAddress("timesep_us", &timesep_us);
  }

  SimConfig *SimConfig::MakeTreeWriter() {
    SimConfig *fdr = new SimConfig();
    fdr->tree = new TTree(fdr->TreeName().c_str(),"");
    fdr->TreeOwned = false;

    fdr->tree->Branch("NXSteps", &fdr->NXSteps, "NXSteps/I");
    fdr->tree->Branch("DetMin", &fdr->DetMin, "DetMin[3]/D");
    fdr->tree->Branch("DetMax", &fdr->DetMax, "DetMax[3]/D");
    fdr->tree->Branch("VetoGap", &fdr->VetoGap, "VetoGap[3]/D");
    fdr->tree->Branch("NMaxTrackSteps", &fdr->NMaxTrackSteps, "NMaxTrackSteps/I");
    fdr->tree->Branch("POTPerFile", &fdr->POTPerFile, "POTPerFile/D");
    fdr->tree->Branch("timesep_us", &fdr->timesep_us, "timesep_us/D");
    fdr->Reset();
    return fdr;
  }

  SimConfig::~SimConfig() {}
