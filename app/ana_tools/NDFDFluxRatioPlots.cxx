#include "CovarianceHelper.hxx"
#include "EffectiveFluxUncertaintyHelper.hxx"
#include "ROOTUtility.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "Eigen/Dense"

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

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }
  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  fhicl::ParameterSet nom_ps = ps.get<fhicl::ParameterSet>("Nominal");

  for (std::string const &beam_config :
       ps.get<std::vector<std::string>>("BeamModes")) {

    std::cout << "[INFO]: Running configured ratios for beam mode: "
              << beam_config << std::endl;

    TDirectory *beam_dir = oupFile->mkdir(beam_config.c_str());

    std::string ND_NominalFile_Template = nom_ps.get<std::string>("NDFileName");
    std::string FD_NominalFile_Template = nom_ps.get<std::string>("FDFileName");

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

      std::unique_ptr<TH2> NDHist_nom_2D =
          GetHistogram_uptr<TH2>(ND_NominalFile, ND_NominalHist);

      std::unique_ptr<TH1> NDHist_nom(
          dynamic_cast<TH1 *>(NDHist_nom_2D->ProjectionX(
              (ND_NominalHist + "_OnAxisND").c_str(), 1, 1)));
      NDHist_nom->SetDirectory(spec_dir);

      std::vector<Int_t> OffAxisBins;
      std::vector<double> OffAxisPositions;
      std::vector<std::unique_ptr<TH1>> OffAxis_Noms;

      for (double OAPos : ps.get<std::vector<double>>("OffAxisPositions")) {
        Int_t OABin = NDHist_nom_2D->GetYaxis()->FindFixBin(OAPos);
        if ((OABin == 0) ||
            (OABin == (NDHist_nom_2D->GetYaxis()->GetNbins() + 1))) {
          continue;
        }
        OffAxisPositions.push_back(OAPos);
        OffAxisBins.push_back(OABin);
      }

      for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
        OffAxis_Noms.emplace_back(dynamic_cast<TH1 *>(
            NDHist_nom_2D->ProjectionX((ND_NominalHist + "_ND_OA" +
                                        std::to_string(OffAxisPositions[OA_it]))
                                           .c_str(),
                                       OffAxisBins[OA_it],
                                       OffAxisBins[OA_it])));
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

      for (fhicl::ParameterSet const &vps :
           ps.get<std::vector<fhicl::ParameterSet>>("Varied", {})) {

        bool IsRelative = vps.get<bool>("IsRelative", false);

        std::string const &VarName = vps.get<std::string>("Name");
        TDirectory *var_dir = spec_dir->mkdir(VarName.c_str());

        std::cout << "[INFO]: \t\tRunning configured ratios for tweak: "
                  << VarName << std::endl;

        std::string ND_File_Template = vps.get<std::string>("NDFileName");
        std::string FD_File_Template = vps.get<std::string>("FDFileName");

        std::string ND_File = str_replace(ND_File_Template, "%C",
                                          std::string("ND_") + beam_config);
        std::string FD_File = str_replace(FD_File_Template, "%C",
                                          std::string("FD_") + beam_config);

        ND_File = str_replace(ND_File, "%N", VarName);
        FD_File = str_replace(FD_File, "%N", VarName);

        std::string ND_Hist_Template = vps.get<std::string>("NDHistName");
        std::string FD_Hist_Template = vps.get<std::string>("FDHistName");

        std::string ND_Hist = str_replace(ND_Hist_Template, "%C",
                                          std::string("ND_") + beam_config);
        std::string FD_Hist = str_replace(FD_Hist_Template, "%C",
                                          std::string("FD_") + beam_config);
        ND_Hist = str_replace(ND_Hist, "%S", species_config);
        FD_Hist = str_replace(FD_Hist, "%S", species_config);
        ND_Hist = str_replace(ND_Hist, "%N", VarName);
        FD_Hist = str_replace(FD_Hist, "%N", VarName);

        std::unique_ptr<TH2> NDHist_var_2D =
            GetHistogram_uptr<TH2>(ND_File, ND_Hist);

        std::unique_ptr<TH1> NDHist_var(dynamic_cast<TH1 *>(
            NDHist_var_2D->ProjectionX((ND_Hist + "_OnAxisND").c_str(), 1, 1)));
        NDHist_var->SetDirectory(var_dir);

        if (IsRelative) {
          for (Int_t bi_it = 0; bi_it < NDHist_var->GetXaxis()->GetNbins();
               ++bi_it) {
            NDHist_var->SetBinContent(
                bi_it + 1, (1 + NDHist_var->GetBinContent(bi_it + 1)) *
                               NDHist_nom->GetBinContent(bi_it + 1));
          }
        }

        std::unique_ptr<TH1> NDHist_rel(dynamic_cast<TH1 *>(
            NDHist_var->Clone((ND_Hist + "_OnAxisND" + "_relative").c_str())));
        NDHist_rel->SetDirectory(var_dir);
        NDHist_rel->Divide(NDHist_nom.get());
        for (Int_t bi_it = 0; bi_it < NDHist_rel->GetXaxis()->GetNbins();
             ++bi_it) {
          NDHist_rel->SetBinContent(bi_it + 1,
                                    NDHist_rel->GetBinContent(bi_it + 1) - 1);
        }

        std::vector<std::unique_ptr<TH1>> OffAxis_var;
        std::vector<std::unique_ptr<TH1>> OffAxis_rel;

        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          OffAxis_var.emplace_back(
              dynamic_cast<TH1 *>(NDHist_var_2D->ProjectionX(
                  (ND_Hist + "_ND_OA" + std::to_string(OffAxisPositions[OA_it]))
                      .c_str(),
                  OffAxisBins[OA_it], OffAxisBins[OA_it])));
          OffAxis_var.back()->SetDirectory(spec_dir);

          if (IsRelative) {
            for (Int_t bi_it = 0;
                 bi_it < OffAxis_var.back()->GetXaxis()->GetNbins(); ++bi_it) {
              OffAxis_var.back()->SetBinContent(
                  bi_it + 1,
                  (1 + OffAxis_var.back()->GetBinContent(bi_it + 1)) *
                      OffAxis_Noms[OA_it]->GetBinContent(bi_it + 1));
            }
          }

          OffAxis_rel.emplace_back(dynamic_cast<TH1 *>(NDHist_var->Clone(
              (ND_Hist + "_ND_OA" + std::to_string(OffAxisPositions[OA_it]) +
               "_relative")
                  .c_str())));
          OffAxis_rel.back()->SetDirectory(var_dir);
          OffAxis_rel.back()->Divide(OffAxis_Noms[OA_it].get());
          for (Int_t bi_it = 0;
               bi_it < OffAxis_rel.back()->GetXaxis()->GetNbins(); ++bi_it) {
            OffAxis_rel.back()->SetBinContent(
                bi_it + 1, OffAxis_rel.back()->GetBinContent(bi_it + 1) - 1);
          }
        }

        std::unique_ptr<TH1> FDHist_var =
            GetHistogram_uptr<TH1>(FD_File, FD_Hist);
        FDHist_var->SetName((FD_Hist + "_FD").c_str());
        FDHist_var->SetDirectory(var_dir);

        if (IsRelative) {
          for (Int_t bi_it = 0; bi_it < FDHist_var->GetXaxis()->GetNbins();
               ++bi_it) {
            FDHist_var->SetBinContent(
                bi_it + 1, (1 + FDHist_var->GetBinContent(bi_it + 1)) *
                               FDHist_nom->GetBinContent(bi_it + 1));
          }
        }

        std::unique_ptr<TH1> FDHist_rel(dynamic_cast<TH1 *>(
            FDHist_var->Clone((FD_Hist + "_FD_relative").c_str())));
        FDHist_rel->SetDirectory(var_dir);
        FDHist_rel->Divide(FDHist_nom.get());
        for (Int_t bi_it = 0; bi_it < FDHist_rel->GetXaxis()->GetNbins();
             ++bi_it) {
          FDHist_rel->SetBinContent(bi_it + 1,
                                    FDHist_rel->GetBinContent(bi_it + 1) - 1);
        }

        std::unique_ptr<TH1> RatioHist_var(dynamic_cast<TH1 *>(
            FDHist_var->Clone(("FD_ND_" + species_config + "_Ratio").c_str())));
        RatioHist_var->Divide(NDHist_var.get());
        RatioHist_var->SetDirectory(var_dir);

        std::unique_ptr<TH1> RatioHist_rel(
            dynamic_cast<TH1 *>(RatioHist_var->Clone(
                ("FD_ND_" + species_config + "_Ratio" + "_relative").c_str())));
        RatioHist_rel->SetDirectory(var_dir);
        RatioHist_rel->Divide(RatioHist_nom.get());
        for (Int_t bi_it = 0; bi_it < RatioHist_rel->GetXaxis()->GetNbins();
             ++bi_it) {
          RatioHist_rel->SetBinContent(
              bi_it + 1, RatioHist_rel->GetBinContent(bi_it + 1) - 1);
        }

        NDHist_var.release();
        NDHist_rel.release();

        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          OffAxis_var[OA_it].release();
          OffAxis_rel[OA_it].release();
        }

        FDHist_var.release();
        FDHist_rel.release();
        RatioHist_var.release();
        RatioHist_rel.release();
      }

      for (fhicl::ParameterSet const &eps :
           ps.get<std::vector<fhicl::ParameterSet>>("EffectiveVariations",
                                                    {})) {

        std::string const &VarName = eps.get<std::string>("Name");
        TDirectory *var_dir = spec_dir->mkdir(VarName.c_str());

        std::cout << "[INFO]: \t\tRunning configured ratios for Decomposed "
                     "error component: "
                  << VarName << std::endl;

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

        std::cout << "\t\t--NParameters: " << eph.GetNParameters() << std::endl;

        size_t NBins = NDHist_nom->GetXaxis()->GetNbins();
        std::vector<std::pair<int, int>> EPHBins;

        for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
          double enu = NDHist_nom->GetXaxis()->GetBinCenter(bi_it + 1);
          EPHBins.emplace_back(eph.GetBin(nu_pdg, enu, 0, 0, true, IsNuMode),
                               eph.GetBin(nu_pdg, enu, 0, 0, false, IsNuMode));
        }

        // oa bin, enu bin
        std::vector<std::vector<int>> EPH_OABins;
        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          EPH_OABins.emplace_back();
          for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
            double enu = NDHist_nom->GetXaxis()->GetBinCenter(bi_it + 1);
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
            RelativeTweaks_ND_OA[OA_it].emplace_back(NBins, 1);
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
              // if (!bi_it) {
              //   std::cout << "\t\t--Excursion for ND bin 0, parameter " <<
              //   p_it
              //             << " = " << NDWeight << "\n";
              // }
              double FDWeight = eph.GetFluxWeight(
                  p_it, pvals[p_it], EPHBins[bi_it].second, NuConfig_FD);
              RelativeTweaks_ND.back()[bi_it] *= NDWeight;
              RelativeTweaks_FD.back()[bi_it] *= FDWeight;
              RelativeTweaks_Ratio.back()[bi_it] *= (NDWeight / FDWeight);

              for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
                RelativeTweaks_ND_OA[OA_it].back()[bi_it] *= eph.GetFluxWeight(
                    p_it, pvals[p_it], EPH_OABins[OA_it][bi_it], NuConfig_ND);
              }
            }

            RelativeTweaks_ND.back()[bi_it] -= 1;
            RelativeTweaks_FD.back()[bi_it] -= 1;
            RelativeTweaks_Ratio.back()[bi_it] -= 1;
            for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
              RelativeTweaks_ND_OA[OA_it].back()[bi_it] -= 1;
            }
          }

          // if (!t_it) {
          //   std::cout << "\t\t--Throw 0 relative vector [";
          //   for (size_t bi_it = 0; bi_it < NBins; ++bi_it) {
          //     std::cout << RelativeTweaks_ND.back()[bi_it]
          //               << (((bi_it + 1) == NBins) ? " ]" : ", ");
          //   }
          //   std::cout << std::endl;
          // }
        }

        CovarianceBuilder cb_ND(RelativeTweaks_ND, true);
        CovarianceBuilder cb_FD(RelativeTweaks_FD, true);
        CovarianceBuilder cb_Ratio(RelativeTweaks_Ratio, true);

        std::unique_ptr<TH1> NDHist_rel(dynamic_cast<TH1 *>(
            NDHist_nom->Clone((ND_NominalHist + "_ND_relative").c_str())));
        FillHistFromEigenVector(NDHist_rel.get(), cb_ND.GetStdDevVector());
        NDHist_rel->SetDirectory(var_dir);

        std::unique_ptr<TH1> FDHist_rel(dynamic_cast<TH1 *>(
            FDHist_nom->Clone((FD_NominalHist + "_FD_relative").c_str())));
        FillHistFromEigenVector(FDHist_rel.get(), cb_FD.GetStdDevVector());
        FDHist_rel->SetDirectory(var_dir);

        std::unique_ptr<TH1> RatioHist_rel(dynamic_cast<TH1 *>(
            RatioHist_nom->Clone((ND_NominalHist + "_ND_FD_Ratio").c_str())));
        FillHistFromEigenVector(RatioHist_rel.get(),
                                cb_Ratio.GetStdDevVector());
        RatioHist_rel->SetDirectory(var_dir);

        for (size_t OA_it = 0; OA_it < OffAxisBins.size(); ++OA_it) {
          CovarianceBuilder cb_OA(RelativeTweaks_ND_OA[OA_it]);

          std::unique_ptr<TH1> OA_rel(
              dynamic_cast<TH1 *>(OffAxis_Noms[OA_it]->Clone(
                  (ND_NominalHist + "_ND_OA" +
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
    }
  }

  oupFile->Write();
  oupFile->Close();
}
