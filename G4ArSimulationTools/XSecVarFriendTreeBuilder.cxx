#include "EDepTreeReader.h"
#include "VALORModelClassifier.h"

#include "Utils.hxx"

#include <string>
#include <vector>

std::string inpfile;
std::string oupfile;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>   : TChain descriptor for"
               " input tree. \n"
               "\t-o <friendoutput.root>     : TFile for output friend tree. \n"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      oupfile = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  EDep edr("EDeps", inpfile);

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");
  TTree *friendtree = new TTree("XSecWeights", "");
  TTree *configtree = new TTree("ConfigTree", "");

  size_t NWeights = 1;
  configtree->Branch("NWeights", &NWeights, "NWeights/I");
  configtree->Fill();

  double *xsecweights = new double[NWeights];
  friendtree->Branch(
      "XSecWeights", xsecweights,
      (std::string("XSecWeights[") + to_str(NWeights) + "]/D").c_str());

  size_t loud_every = edr.GetEntries() / 10;
  size_t NFills = 0;
  Long64_t NEntries = edr.GetEntries();
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    std::string evc = edr.EventCode->GetString().Data();

    std::vector<VALORModel::TrueClass> appdial = GetApplicableDials(edr);

    std::cout << "[INFO]: event " << e_it << ", code: " << evc << std::endl;
    size_t d_it = 0;
    xsecweights[0] = 1;
    for (auto &d : appdial) {
      std::cout << "\t [" << d_it << "]: " << d
                << ", +1 sigma weight: " << GetVALORWeight(d, 1, edr)
                << std::endl;
      xsecweights[0] *= GetVALORWeight(d, 1, edr);
      d_it++;
    }

    NFills++;
    friendtree->Fill();
  }
  of->Write();
  of->Close();
}
