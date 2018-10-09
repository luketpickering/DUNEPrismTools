#include "DetectorStop.hxx"
#include "EDepTreeReader.h"
#include "TTree.h"
#include "Utils.hxx"

#include <string>
#include <vector>

std::string inpfile;
std::string nominphistfile, nominphistname;
std::vector<std::string> varinphistfile, varinphistname;
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

      varinphistfile.push_back(params[0]);
      varinphistname.push_back(params[1]);

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
                 "descriptor like: \"weighthist.root,weighthistname\"."
              << std::endl;
    return 1;
  }
  // Load weighting TH2
  TH2 *weightinghist_nom = GetHistogram<TH2D>(nominphistfile, nominphistname);
  std::vector<TH2 *> weightinghist_var;

  size_t nvar = varinphistfile.size();
  for (size_t i = 0; i < nvar; ++i) {
    weightinghist_var.push_back(
        GetHistogram<TH2D>(varinphistfile[i], varinphistname[i]));
  }

  EDep edr("EDeps", inpfile);

  // Make output tree
  TFile *oupf = CheckOpenFile(oupfile, "RECREATE");
  TTree *oupt = new TTree("FluxVarEDepsFriend", "");
  std::vector<double> weight;
  weight.resize(varinphistfile.size());
  for (size_t i = 0; i < nvar; ++i) {
    oupt->Branch((std::string("weight_var") + to_str(i)).c_str(), &weight[i],
                 (std::string("weight_var") + to_str(i) + "/D").c_str());

    // Gen weights
    size_t NInputEntries = edr.GetEntries();
    for (size_t e_it = 0; e_it < NInputEntries; ++e_it) {
      edr.GetEntry(e_it);

      Int_t ebin = weightinghist_nom->GetXaxis()->FindFixBin(edr.nu_4mom[3]);
      Int_t xbin =
          weightinghist_nom->GetYaxis()->FindFixBin(fabs(edr.vtx[0]) / 100.0);

      for (size_t i = 0; i < nvar; ++i) {
        weight[i] = weightinghist_var[i]->GetBinContent(ebin, xbin) /
                    weightinghist_nom->GetBinContent(ebin, xbin);
      }
      oupt->Fill();
    }

    oupf->Write();
    oupf->Close();
  }
}
