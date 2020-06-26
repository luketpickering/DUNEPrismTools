#include "CovarianceHelper.hxx"
#include "EffectiveFluxUncertaintyHelper.hxx"
#include "ROOTUtility.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "Eigen/Dense"

#include "TH2Jagged.h"

#include <memory>
#include <random>

// Expect a fhicl like
// Nominal: {
// //   %C is replaced by ND_{nu,nubar}
//   NDFileName: %CFlux.root
// //   %S is replaced by the species name, {numu, numubar, nue, nuebar}
//   NDHistName: LBNF_%S_flux_Nom
//   FDFileName: %CFlux.root
//   FDHistName: LBNF_%S_flux_Nom
// }
//
// Varied: [ {
// //   Whether this should be applied to the Nominal fluxes or used as is.
//   IsRelative: false
//   Name: "var_name"
//   NDFileName: %CFlux.root
//   NDHistName: LBNF_%S_flux_Nom
//   FDFileName: %CFlux.root
//   FDHistName: LBNF_%S_flux_Nom
// }, ]
//
// EffectiveVariations: [
//   {
//     NThrows: 10000
//     Name: "Total"
//     NEffectiveParametersToRead:20
//     ND_detector_tag: "ND"
//     FD_detector_tag: "FD"
//     nu_mode_beam_tag: "nu"
//     nubar_mode_beam_tag: "nubar"
//     numu_species_tag: "numu"
//     nue_species_tag: "nue"
//     numubar_species_tag: "numubar"
//     nuebar_species_tag: "nuebar"
//     InputFile: "FluxCovmat_Total.root"
//   }
// ]
int main(int argc, char const *argv[]) {

  if (argc < 2 || argc > 5 || (std::string(argv[1]) == "-?") ||
      (std::string(argv[1]) == "--help")) {
    std::cout << "[RUNLIKE]: " << argv[0]
              << " config.fcl [input_file.root] [output_file.root]"
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string def_ifile = "";
  if (argc > 2) {
    def_ifile = argv[2];
    std::cout << "Reading nominal from " << def_ifile << std::endl;
  }

  std::string def_efile = "";
  if (argc > 3) {
    def_efile = argv[3];
    std::cout << "Reading errors from " << def_efile << std::endl;
  }

  std::string def_ofile = "";
  if (argc > 4) {
    def_ofile = argv[4];
    std::cout << "Writing to " << def_ofile << std::endl;
  }

  std::string oupf = ps.get<std::string>("OutputFile", def_ofile);
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  fhicl::ParameterSet nom_ps = ps.get<fhicl::ParameterSet>("Nominal");

  for (std::string const &beam_config :
       ps.get<std::vector<std::string>>("BeamModes")) {

    std::cout << "[INFO]: Running configured ratios for beam mode: "
              << beam_config << std::endl;

    TDirectory *beam_dir = oupFile->mkdir(beam_config.c_str());

    std::string ND_NominalFile_Template =
        nom_ps.get<std::string>("NDFileName", def_ifile);
    std::string FD_NominalFile_Template =
        nom_ps.get<std::string>("FDFileName", def_ifile);

    std::string ND_NominalFile = str_replace(ND_NominalFile_Template, "%C",
                                             std::string("ND_") + beam_config);
    std::string FD_NominalFile = str_replace(FD_NominalFile_Template, "%C",
                                             std::string("FD_") + beam_config);

    for (std::string const &species_config :
         ps.get<std::vector<std::string>>("Species")) {

      std::cout << "[INFO]: \tRunning configured ratios for neutrino species: "
                << species_config << std::endl;

      TDirectory *spec_dir = beam_dir->mkdir(species_config.c_str());

      std::string ND_NominalHist_Template =
          nom_ps.get<std::string>("NDHistName");
      std::string FD_NominalHist_Template =
          nom_ps.get<std::string>("FDHistName");

      std::string ND_NominalHist =
          str_replace(ND_NominalHist_Template, "%S", species_config);
      std::string FD_NominalHist =
          str_replace(FD_NominalHist_Template, "%S", species_config);

      ND_NominalHist =
          str_replace(ND_NominalHist, "%C", std::string("ND_") + beam_config);
      FD_NominalHist =
          str_replace(FD_NominalHist, "%C", std::string("FD_") + beam_config);

      std::unique_ptr<TH2JaggedD> NDHist_nom_2D =
          GetHistogram_uptr<TH2JaggedD>(ND_NominalFile, ND_NominalHist);

      // assumes 1m binning from where bin 1  is -3 -- -2
      std::unique_ptr<TH1> NDHist_nom(
          NDHist_nom_2D
              ->UniformRange((ND_NominalHist + "_OnAxisND").c_str(), 2, 6, true)
              ->NonUniformSlice(1));
      NDHist_nom->SetDirectory(spec_dir);

      std::vector<Int_t> OffAxisBins;
      std::vector<double> OffAxisPositions;
      std::vector<std::unique_ptr<TH1>> OffAxis_Noms;

      for (double OAPos : ps.get<std::vector<double>>("OffAxisPositions")) {
        Int_t OABin = NDHist_nom_2D->GetYaxis()->FindFixBin(OAPos);
        if ((OABin == 0) ||
            (OABin == (NDHist_nom_2D->GetYaxis()->GetNbins() + 1))) {
          std::cout << "Couldn't find y bin for oa pos = " << OAPos
                    << std::endl;
          throw;
          continue;
        }
        OffAxisPositions.push_back(OAPos);
        OffAxisBins.push_back(OABin);
      }

      for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
        OffAxis_Noms.emplace_back(
            NDHist_nom_2D
                ->UniformRange((ND_NominalHist + "_ND_OA" +
                                std::to_string(OffAxisPositions[OA_it]))
                                   .c_str(),
                               OffAxisBins[OA_it] - 2, OffAxisBins[OA_it] + 2,
                               true)
                ->NonUniformSlice(1));
        OffAxis_Noms.back()->SetDirectory(spec_dir);
      }

      std::unique_ptr<TH1> FDHist_nom =
          GetHistogram_uptr<TH1>(FD_NominalFile, FD_NominalHist);
      FDHist_nom->SetName((FD_NominalHist + "_FD").c_str());
      FDHist_nom->SetDirectory(spec_dir);

      std::unique_ptr<TH1> RatioHist_nom(dynamic_cast<TH1 *>(FDHist_nom->Clone(
          ("FD_ND_" + species_config + "_Ratio_Nom").c_str())));
      RatioHist_nom->Divide(NDHist_nom.get());
      RatioHist_nom->SetDirectory(spec_dir);

      for (fhicl::ParameterSet eps :
           ps.get<std::vector<fhicl::ParameterSet>>("VariationSets", {})) {

        std::string const &VarName = eps.get<std::string>("Name");
        TDirectory *var_dir = spec_dir->mkdir(VarName.c_str());

        std::cout << "[INFO]: \t\tRunning configured ratios for Decomposed "
                     "error component: "
                  << VarName << std::endl;

        if (!eps.has_key("InputFile") && def_efile.size()) {
          eps.put("InputFile", def_efile);
        }

        // Instantiate Helper object
        EffectiveFluxUncertaintyHelper eph;
        eph.Initialize(eps);

        int nu_pdg;
        if (species_config == "numu") {
          nu_pdg = 14;
        } else if (species_config == "nue") {
          nu_pdg = 12;
        } else if (species_config == "numubar") {
          nu_pdg = -14;
        } else if (species_config == "nuebar") {
          nu_pdg = -12;
        }
        bool IsNuMode = false;
        if (beam_config == "nu") {
          IsNuMode = true;
        }

        if (!eph.GetNParameters()) {
          std::cout << "[ERROR]: Found no parameters." << std::endl;
          abort();
        }

        std::cout << "\t\t--NParameters: " << eph.GetNParameters() << std::endl;

        size_t NBins = NDHist_nom->GetXaxis()->GetNbins();
        std::vector<std::pair<int, int>> EPHBins;

        for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
          double enu = NDHist_nom->GetXaxis()->GetBinCenter(bi_it + 1);
          EPHBins.emplace_back(eph.GetBin(nu_pdg, enu, 0, 0, true, IsNuMode),
                               eph.GetBin(nu_pdg, enu, 0, 0, false, IsNuMode));
        }

        std::vector<size_t> OANBins;
        // Using the 1D projections means we don't have to take care over jagged
        // TH2s
        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          OANBins.push_back(OffAxis_Noms[OA_it]->GetXaxis()->GetNbins());
        }

        // oa bin, enu bin
        std::vector<std::vector<int>> EPH_OABins;
        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          EPH_OABins.emplace_back();
          for (size_t bi_it = 0; bi_it < OANBins[OA_it]; ++bi_it) {
            double enu =
                OffAxis_Noms[OA_it]->GetXaxis()->GetBinCenter(bi_it + 1);
            EPH_OABins.back().emplace_back(eph.GetBin(
                nu_pdg, enu, OffAxisPositions[OA_it], 0, true, IsNuMode));
          }
        }

        int NuConfig_ND = eph.GetNuConfig(nu_pdg, true, IsNuMode);
        int NuConfig_FD = eph.GetNuConfig(nu_pdg, false, IsNuMode);

        std::random_device r;

        std::mt19937 RNEngine(r());
        std::normal_distribution<double> RNJesus(0, 1);

        std::vector<std::vector<double>> RelativeTweaks_ND;
        // pos,throw,enubin
        std::vector<std::vector<std::vector<double>>> RelativeTweaks_ND_OA(
            OffAxisBins.size());
        std::vector<std::vector<double>> RelativeTweaks_FD;
        std::vector<std::vector<double>> RelativeTweaks_Ratio;
        size_t NThrows = eps.get<size_t>("NThrows");
        for (size_t t_it = 0; t_it < NThrows; ++t_it) {
          RelativeTweaks_ND.emplace_back(NBins, 1);
          RelativeTweaks_FD.emplace_back(NBins, 1);
          RelativeTweaks_Ratio.emplace_back(NBins, 1);
          for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
            RelativeTweaks_ND_OA[OA_it].emplace_back(OANBins[OA_it], 1);
          }

          std::vector<double> pvals(eph.GetNParameters());
          for (size_t p_it = 0; p_it < pvals.size(); ++p_it) {
            pvals[p_it] = NThrows == 1 ? 1 : RNJesus(RNEngine);
          }

#pragma omp parallel for
          for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
            for (size_t p_it = 0; p_it < pvals.size(); ++p_it) {

              double NDWeight = eph.GetFluxWeight(
                  p_it, pvals[p_it], EPHBins[bi_it].first, NuConfig_ND);
              double FDWeight = eph.GetFluxWeight(
                  p_it, pvals[p_it], EPHBins[bi_it].second, NuConfig_FD);

              RelativeTweaks_ND.back()[bi_it] *= NDWeight;
              RelativeTweaks_FD.back()[bi_it] *= FDWeight;
              RelativeTweaks_Ratio.back()[bi_it] *= (NDWeight / FDWeight);
            }

            RelativeTweaks_ND.back()[bi_it] -= 1;
            RelativeTweaks_FD.back()[bi_it] -= 1;
            RelativeTweaks_Ratio.back()[bi_it] -= 1;
          }

          for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
#pragma omp parallel for
            for (size_t bi_it = 0; bi_it < OANBins[OA_it]; ++bi_it) {
              for (size_t p_it = 0; p_it < pvals.size(); ++p_it) {
                RelativeTweaks_ND_OA[OA_it].back()[bi_it] *= eph.GetFluxWeight(
                    p_it, pvals[p_it], EPH_OABins[OA_it][bi_it], NuConfig_ND);
              }
              RelativeTweaks_ND_OA[OA_it].back()[bi_it] -= 1;
            }
          }
        }

        CovarianceBuilder cb_ND(RelativeTweaks_ND, true);
        CovarianceBuilder cb_FD(RelativeTweaks_FD, true);
        CovarianceBuilder cb_Ratio(RelativeTweaks_Ratio, true);

        std::string HistNameStub = "LBNF_" + species_config + "_flux";

        std::unique_ptr<TH1> NDHist_rel(dynamic_cast<TH1 *>(
            NDHist_nom->Clone((HistNameStub + "_ND_relative").c_str())));
        FillHistFromEigenVector(NDHist_rel.get(), cb_ND.GetStdDevVector());
        NDHist_rel->SetDirectory(var_dir);

        std::unique_ptr<TH1> FDHist_rel(dynamic_cast<TH1 *>(
            FDHist_nom->Clone((HistNameStub + "_FD_relative").c_str())));
        FillHistFromEigenVector(FDHist_rel.get(), cb_FD.GetStdDevVector());
        FDHist_rel->SetDirectory(var_dir);

        std::unique_ptr<TH1> RatioHist_rel(dynamic_cast<TH1 *>(
            RatioHist_nom->Clone((HistNameStub + "_ND_FD_Ratio").c_str())));
        FillHistFromEigenVector(RatioHist_rel.get(),
                                cb_Ratio.GetStdDevVector());
        RatioHist_rel->SetDirectory(var_dir);

        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          CovarianceBuilder cb_OA(RelativeTweaks_ND_OA[OA_it]);

          std::unique_ptr<TH1> OA_rel(
              dynamic_cast<TH1 *>(OffAxis_Noms[OA_it]->Clone(
                  (HistNameStub + "_ND_OA" +
                   std::to_string(OffAxisPositions[OA_it]) + "_relative")
                      .c_str())));
          FillHistFromEigenVector(OA_rel.get(), cb_OA.GetStdDevVector());
          OA_rel->SetDirectory(var_dir);

          OA_rel.release();
        }

        NDHist_rel.release();
        FDHist_rel.release();
        RatioHist_rel.release();
      }

      NDHist_nom.release();
      FDHist_nom.release();
      RatioHist_nom.release();

      for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
        OffAxis_Noms[OA_it].release();
      }
    } // end species
  }   // end beam modes

  oupFile->Write();
  oupFile->Close();
}
