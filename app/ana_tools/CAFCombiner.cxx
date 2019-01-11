#include "FileSystemUtility.hxx"

#include "CAFReader.hxx"

#include "DUNETDRDetHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

// Expected FHiCL like
// Inputs: [ {
//   InputDirectory: "/path/to/files"
//   FilePattern: "file_pattern*.root"
//   NMaxEvents: 1E6
//   POTPerFileOverride: 1E17
// }, ]
// OutputFile: "Output_CAF.root"

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string oupf = ps.get<std::string>("OutputFile");

  CAFReader *wrtr = nullptr;

  for (fhicl::ParameterSet input :
       ps.get<std::vector<fhicl::ParameterSet>>("Inputs")) {
    std::string idir =
        EnsureTrailingSlash(input.get<std::string>("InputDirectory"));
    std::string pat = input.get<std::string>("FilePattern");

    std::vector<std::string> MatchingFileList = GetMatchingFiles(idir, pat);

    size_t NMaxEvents =
        input.get<size_t>("NMaxEvents", std::numeric_limits<size_t>::max());
    double POTPerFileOverride = input.get<double>(
        "POTPerFileOverride", std::numeric_limits<double>::max());

    for (std::string const &inpf : MatchingFileList) {
      CAFReader rdr(idir + inpf);

      size_t fents = std::min(NMaxEvents, rdr.GetEntries());
      double POTScale = double(fents) / double(rdr.GetEntries());

      if (fents) {
        std::cout << "[INFO]: Reading " << fents << "/" << rdr.GetEntries()
                  << " CAF events from file " << inpf
                  << " with POTScale = " << POTScale << std::endl;
      }

      size_t nsel = 0;
      size_t nfail_cath = 0;
      size_t nfail_wall = 0;
      size_t nfail_desert = 0;
      size_t nfail_fv = 0;
      for (size_t ent = 0; ent < fents; ++ent) {
        rdr.GetEntry(ent);

        if (!wrtr) {
          wrtr = CAFReader::MakeWriter(oupf, rdr.isNDFile);
        }

        if (ent == 0) {
          if (POTPerFileOverride != std::numeric_limits<double>::max()) {

            std::cout << "[INFO]: POT Overriden to " << POTPerFileOverride
                      << std::endl;
            rdr.FilePOT = POTPerFileOverride * POTScale;
          } else {
            std::cout << "[INFO]: File reported total POT = " << rdr.FilePOT
                      << std::endl;
            rdr.FilePOT *= POTScale;
          }
        }

        (*wrtr) = rdr;

        if (ent == 0) {
          wrtr->NewFile(rdr.isFD ? FD_XFV_Select : ND_XFV_Select);
        }

        if (FV_Select(rdr)) {
          wrtr->Fill();
          nsel++;
        } else if (!rdr.isFD) {
          if (!rdr.isFD && !PRISMVertex_Select(rdr.vtx_x)) {
            nfail_desert += 1;
          } else if (rdr.isFD ? FD_IsWall_Select(rdr.vtx_x)
                              : ND_IsWall_Select(rdr.vtx_x)) {
            nfail_wall += 1;
          } else if (rdr.isFD ? FD_IsCathode_Select(rdr.vtx_x)
                              : ND_IsCathode_Select(rdr.vtx_x)) {
            nfail_cath += 1;
          } else {
            nfail_fv += 1;
          }
        }
      }

      if (fents) {
        std::cout << "\tSelected " << nsel << "/" << fents << std::endl;
        std::cout << "\tNFail Vtx Desert " << nfail_desert << std::endl;
        std::cout << "\tNFail Module Wall " << nfail_wall << std::endl;
        std::cout << "\tNFail Module Cathode " << nfail_cath << std::endl;
        std::cout << "\tNFail FV " << nfail_fv << std::endl;
      }
    }
  }

  std::cout << "[INFO]: Wrote " << wrtr->GetEntries() << " events from "
            << " " << wrtr->GetNFiles() << " input files" << std::endl;

  wrtr->file->Write();
  wrtr->file->Close();
}
