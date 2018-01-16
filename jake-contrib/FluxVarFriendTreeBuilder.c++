#include "DetectorStop.hxx"
#include "TTree.h"
#include "Utils.hxx"

#include <string>
#include <vector>

std::string inpfile;
std::string nominphistfile, nominphistname;
std::string varinphistfile, varinphistname;
std::string oupfile;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>             : TChain descriptor for"
               " input tree. \n"
               "\t-hn <weighthist.root,weighthistname> : Input file and "
               "histogram name to use for nominal flux in Enu:Lateral offset.\n"
               "\t-hv <weighthist.root,weighthistname> : Input file and "
               "histogram name to use for varied flux in Enu:Lateral offset.\n"
               "\t-o <outputfile.root>                 : Output file to write "
               "friend tree to.\n"
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
    } else if (std::string(argv[opt]) == "-hn") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      if (params.size() != 2) {
        std::cout << "[ERROR]: Expected to find "
                     "\"weighthist.root,weighthistname\" passed to -h, but "
                     "found\""
                  << argv[opt - 1] << "\"" << std::endl;
        SayUsage(argv);
        exit(1);
      }

      nominphistfile = params[0];
      nominphistname = params[1];

    } else if (std::string(argv[opt]) == "-hv") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      if (params.size() != 2) {
        std::cout << "[ERROR]: Expected to find "
                     "\"weighthist.root,weighthistname\" passed to -h, but "
                     "found\""
                  << argv[opt - 1] << "\"" << std::endl;
        SayUsage(argv);
        exit(1);
      }

      varinphistfile = params[0];
      varinphistname = params[1];

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

  if (!nominphistfile.size() || !nominphistname.size()) {
    std::cout << "[ERROR]: Expected to recieve an input nominal flux histogram "
                 "descriptor like: \"weighthist.root,weighthistname\", but "
                 "found: \""
              << nominphistfile << "," << nominphistname << "\"." << std::endl;
    return 1;
  }

  if (!varinphistfile.size() || !varinphistname.size()) {
    std::cout << "[ERROR]: Expected to recieve an input varied flux histogram "
                 "descriptor like: \"weighthist.root,weighthistname\", but "
                 "found: \""
              << varinphistfile << "," << varinphistname << "\"." << std::endl;
    return 1;
  }
  // Load weighting TH2
  TH2 *weightinghist_nom = GetHistogram<TH2D>(nominphistfile, nominphistname);
  TH2 *weightinghist_var = GetHistogram<TH2D>(varinphistfile, varinphistname);

  // Hookup input tree
  TChain *Input = new TChain("EDeps");
  int NAdded = Input->Add(inpfile.c_str());

  if (!NAdded) {
    std::cout << "[ERROR]: Failed to add \"" << inpfile
              << "\" to input TChain -- does it exist and does it include the "
                 "\"EDeps\" TTree."
              << std::endl;
    return 1;
  }

  size_t NInputEntries = Input->GetEntries();
  if (!NAdded) {
    std::cout << "[ERROR]: Found no EDeps entries in:\"" << inpfile << "\"."
              << std::endl;
    return 1;
  }

  double Enu;
  double vtx[3];
  Input->SetBranchAddress("Enu", &Enu);
  Input->SetBranchAddress("vtx", &vtx);

  // Make output tree
  TFile *oupf = CheckOpenFile(oupfile, "RECREATE");
  TTree *oupt = new TTree("FluxVarEDepsFriend", "");
  double weight;
  oupt->Branch("weight", &weight, "weight/D");

  // Gen weights
  for (size_t e_it = 0; e_it < NInputEntries; ++e_it) {
    Input->GetEntry(e_it);

    Int_t ebin = weightinghist_nom->GetXaxis()->FindFixBin(Enu);
    Int_t xbin =
        weightinghist_nom->GetYaxis()->FindFixBin(fabs(vtx[0]) / 100.0);

    weight = weightinghist_var->GetBinContent(ebin, xbin) /
             weightinghist_nom->GetBinContent(ebin, xbin);
    oupt->Fill();
  }

  oupf->Write();
  oupf->Close();
}
