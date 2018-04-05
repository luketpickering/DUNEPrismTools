#ifndef ITREEREADER_HXX_SEEN
#define ITREEREADER_HXX_SEEN

#include "ROOTUtility.hxx"

#include "TTree.h"

#include <string>

struct ITreeReader {

  ITreeReader(ITreeReader const &) = delete;
  ITreeReader &operator=(ITreeReader const &) = delete;
  ITreeReader(ITreeReader &) = delete;
  ITreeReader &operator=(ITreeReader &) = delete;

  ITreeReader() : tree(nullptr), NFiles(0), NEntries(0), CEnt (0),
    TreeOwned(false) {};

  void LoadTree(std::string const &inputFile) {
      tree = OpenTChainWithFileList(TreeName(), inputFile);

      if(!tree){
        std::cout << "[" << TreeName() << "]: Failed to read input tree from "
          "file." << std::endl;
          throw;
      }
      TreeOwned = true;

      NEntries = tree->GetEntries();
      SetBranchAddresses();
      std::cout << "[" << TreeName() << "]: Loaded TChain with " << NEntries
        << " entries." << std::endl;

      if(NEntries > 0){
        GetEntry(0);
      }
  }

  TTree *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;
  bool TreeOwned;

  virtual std::string TreeName() = 0;
  virtual void Reset() = 0;
  virtual void SetBranchAddresses() = 0;

  void GetEntry(UInt_t e) {
    CEnt = e;
    tree->GetEntry(CEnt);
  }

  UInt_t GetEntry() { return CEnt; }
  UInt_t GetEntries() { return NEntries; }

  void ReleaseInputFile() { if (TreeOwned){ delete tree; } tree = nullptr; }
  void Fill() {
    if (!TreeOwned){ tree->Fill(); } else {
      std::cout << "[WARN]: Attempted to call fill on an input tree of type: "
        << TreeName() << "." << std::endl;
    }
  }
  void SetDirectory(TDirectory *tdir){ tree->SetDirectory(tdir); }

  virtual ~ITreeReader() { ReleaseInputFile(); }
};

#endif
