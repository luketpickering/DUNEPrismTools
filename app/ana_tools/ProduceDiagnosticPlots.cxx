TH1D *CoeffWeightingHelper = nullptr;
TH2D *OffSetVsEProxy = nullptr;
TH2D *OffSetVsETrue = nullptr;
std::vector<TH1D *> SliceEProxy_True_LinCombWeighted;
std::vector<TH1D *> SliceETrue_True;
std::vector<TH1D *> SliceERec_Selected;
std::vector<TH1D *> SliceEMu_Selected;
std::vector<TH1D *> SliceEMu_True;
std::vector<TH1D *> SliceEHadr_Selected;
std::vector<TH1D *> SliceETHadr_Selected;

if (lincombfile.size()) {
  TChain *CoeffTree = new TChain("CoeffTree");

  CoeffTree->Add(lincombfile.c_str());

  double XRange[2];
  double Coeff;

  CoeffTree->SetBranchAddress("XRange", &XRange);
  CoeffTree->SetBranchAddress("Coeff", &Coeff);

  std::vector<double> XRangeBins;
  std::vector<double> Coeffs;

  std::cout << "[INFO]: XRange bins: " << std::flush;
  CoeffTree->GetEntry(0);
  XRangeBins.push_back(XRange[0]);
  Coeffs.push_back(Coeff);
  XRangeBins.push_back(XRange[1]);
  std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
  for (Long64_t i = 1; i < CoeffTree->GetEntries(); ++i) {
    CoeffTree->GetEntry(i);

    // If non-contiguous, must push an empty bit between.
    if (fabs(XRangeBins.back() - XRange[0]) > 1E-5) {
      Coeffs.push_back(0);
      XRangeBins.push_back(XRange[0]);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }

    Coeffs.push_back(Coeff);
    XRangeBins.push_back(XRange[1]);
    std::cout << ", " << XRangeBins.back() << std::flush;
  }
  std::cout << std::endl;

  CoeffWeightingHelper = new TH1D("CoeffWeightingHelper", "",
                                  (XRangeBins.size() - 1), XRangeBins.data());
  CoeffWeightingHelper->SetDirectory(nullptr);

  for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
    CoeffWeightingHelper->SetBinContent(bin_it, Coeffs[bin_it - 1]);

    SliceEProxy_True_LinCombWeighted.push_back(
        new TH1D((std::string("SliceEProxy_True_LinCombWeighted_Slice") +
                  to_str(bin_it - 1))
                     .c_str(),
                 ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));

    SliceETrue_True.push_back(new TH1D(
        (std::string("SliceETrue_True") + to_str(bin_it - 1)).c_str(),
        ";E_{#nu} (GeV);Count", 100, 0, 10));
    SliceERec_Selected.push_back(new TH1D(
        (std::string("SliceERec_Selected") + to_str(bin_it - 1)).c_str(),
        ";E_{Rec} (GeV);Count", 100, 0, 10));
    SliceEMu_Selected.push_back(new TH1D(
        (std::string("SliceEMu_Selected") + to_str(bin_it - 1)).c_str(),
        ";E_{#mu} (GeV);Count", 100, 0, 10));
    SliceEMu_True.push_back(
        new TH1D((std::string("SliceEMu_True") + to_str(bin_it - 1)).c_str(),
                 ";E_{#mu} (GeV);Count", 100, 0, 10));
    SliceEHadr_Selected.push_back(new TH1D(
        (std::string("SliceEHadr_Selected") + to_str(bin_it - 1)).c_str(),
        ";E_{Hadr} (GeV);Count", 100, 0, 10));
    SliceETHadr_Selected.push_back(new TH1D(
        (std::string("SliceETHadr_Selected") + to_str(bin_it - 1)).c_str(),
        ";E_{Hadr} (GeV);Count", 100, 0, 10));
  }

  OffSetVsEProxy =
      new TH2D("OffSetVsEProxy", ";E_{#nu,Proxy} (GeV);Offset (cm);Count",
               100, 0, 10, XRangeBins.size() - 1, XRangeBins.data());
  OffSetVsETrue =
      new TH2D("OffSetVsETrue", ";E_{#nu} (GeV);Offset (cm);Count", 100, 0,
               10, XRangeBins.size() - 1, XRangeBins.data());

  OffSetVsEProxy->SetDirectory(nullptr);
  OffSetVsETrue->SetDirectory(nullptr);
}




  TH1D *EventRates_Selected =
      new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
  TH1D *EventRates_True =
      new TH1D("EventRates_True", ";Offset (cm);Count", 400, -250, 3750);

  TH1D *ERec_Selected = new TH1D(
      "ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count", 100, 0, 10);
  TH1D *ETrue_Selected =
      new TH1D("ETrue_Selected", ";E_{#nu} (GeV);Count", 100, 0, 10);
  TH1D *EHadr_Selected =
      new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count", 100, 0, 10);
  TH1D *EHadr_True =
      new TH1D("EHadr_True", ";E_{Hadr} (GeV);Count", 100, 0, 10);

  TH2D *EvRateERec_Selected = new TH2D(
      "EvRateERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
      100, 0, 10, 400, -250, 3750);

  TH2D *EvRateEHadr_Selected =
      new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);
  TH2D *EvRateEHadr_True =
      new TH2D("EvRateEHadr_True", ";E_{Hadr} (GeV);Offset (cm)", 100, 0, 10,
               400, -250, 3750);

  TH2D *EvRateEProxy_Selected =
      new TH2D("EvRateEProxy_Selected", ";E_{#nu,proxy} (GeV);Offset (cm)", 100,
               0, 10, 400, -250, 3750);
  TH2D *EvRateEProxy_True =
      new TH2D("EvRateEProxy_True", ";E_{#nu,proxy} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);

  TH2D *EvRateEVeto = new TH2D("EvRateEVeto", ";E_{Hadr} (GeV);Offset (cm)",
                               100, 0, 2, 400, -250, 3750);

  TH2D *EProxyERec_Selected =
      new TH2D("EProxyERec_Selected", ";E_{#nu,proxy} (GeV);E_{Rec} (GeV)", 100,
               1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueEProxy_Selected =
      new TH2D("ETrueEProxy_Selected", ";E_{#nu} (GeV);E_{#nu,proxy} (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERec_Selected =
      new TH2D("ETrueERec_Selected", ";E_{#nu} (GeV);E_{Rec} (GeV)", 100, 1E-5,
               10, 100, 1E-5, 10);

  TH2D *ETrueEmu_Selected =
      new TH2D("ETrueEmu_Selected", ";E_{#nu} (GeV);E_{#mu} (GeV)", 100, 1E-5,
               10, 100, 1E-5, 10);

  TH2D *ETrueERecNeutron_Selected = new TH2D(
      "ETrueERecNeutron_Selected", ";E_{#nu} (GeV);E_{Rec} Neutron (GeV)", 100,
      1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecNeutron_timesep_Selected = new TH2D(
      "ETrueERecNeutron_timesep_Selected",
      ";E_{#nu} (GeV);E_{Rec} Neutron (GeV)", 100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecProton_Selected = new TH2D(
      "ETrueERecProton_Selected", ";E_{#nu} (GeV);E_{Rec} Proton (GeV)", 100,
      1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecProton_timesep_Selected = new TH2D(
      "ETrueERecProton_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Proton (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecPion_Selected =
      new TH2D("ETrueERecPion_Selected", ";E_{#nu} (GeV);E_{Rec} Pion (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecPion_timesep_Selected = new TH2D(
      "ETrueERecPion_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Pion (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecOther_Selected =
      new TH2D("ETrueERecOther_Selected", ";E_{#nu} (GeV);E_{Rec} Other (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecOther_timesep_Selected = new TH2D(
      "ETrueERecOther_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Other (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TDirectory *oupD = of;

  TDirectory *wD = oupD->mkdir("StopEventRates");
  wD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEProxy_Selected;
  std::vector<TH1D *> StopEProxy_True;
  std::vector<TH1D *> StopETrue_Selected;
  std::vector<TH1D *> StopETrue_True;
  std::vector<TH1D *> StopERec_Selected;
  std::vector<TH1D *> StopEventRates_True;
  std::vector<TH1D *> StopEventRates_TrueUnWeighted;

  for (Int_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEProxy_Selected.push_back(new TH1D(
        (std::string("StopEProxy_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEProxy_True.push_back(new TH1D(
        (std::string("StopEProxy_True_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEventRates_True.push_back(new TH1D(
        (std::string("StopEventRates_True_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_TrueUnWeighted.push_back(new TH1D(
        (std::string("StopEventRates_TrueUnWeighted_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));

    StopETrue_Selected.push_back(
        new TH1D((std::string("StopETrue_Selected") + to_str(stop_it)).c_str(),
                 ";E_{#nu} (GeV);Count", 100, 0, 10));
    StopETrue_True.push_back(
        new TH1D((std::string("StopETrue_True") + to_str(stop_it)).c_str(),
                 ";E_{#nu} (GeV);Count", 100, 0, 10));
    StopERec_Selected.push_back(
        new TH1D((std::string("StopERec_Selected") + to_str(stop_it)).c_str(),
                 ";E_{Rec} (GeV);Count", 100, 0, 10));
  }
  oupD->cd();


  // Fill truth histos
  EventRates_True->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EHadr_True->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                   DepSumRdr.stop_weight);
  EvRateEHadr_True->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                         -1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EvRateEProxy_True->Fill(DepSumRdr.ERecProxy_True, -1 * DepSumRdr.vtx[0],
                          DepSumRdr.stop_weight);
  StopEProxy_True[DepSumRdr.stop]->Fill(DepSumRdr.ERecProxy_True, DepSumRdr.stop_weight);
  StopEventRates_True[DepSumRdr.stop]->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  StopEventRates_TrueUnWeighted[DepSumRdr.stop]->Fill(-1 * DepSumRdr.vtx[0]);
  EvRateEVeto->Fill(DepSumRdr.TotalNonlep_Dep_veto, -1 * DepSumRdr.vtx[0]);
  StopETrue_True[DepSumRdr.stop]->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);

  Int_t xb = -1;
  if (CoeffWeightingHelper) {
    OffSetVsEProxy->Fill(DepSumRdr.ERecProxy_True, -1 * DepSumRdr.vtx[0]);
    OffSetVsETrue->Fill(DepSumRdr.nu_4mom[3], -1 * DepSumRdr.vtx[0]);

    xb = CoeffWeightingHelper->GetXaxis()->FindFixBin(-1 * DepSumRdr.vtx[0]);
    double CoeffWeight = CoeffWeightingHelper->GetBinContent(xb);

    if (xb && (xb <= Int_t(SliceEProxy_True_LinCombWeighted.size()))) {
      SliceEProxy_True_LinCombWeighted[xb - 1]->Fill(
          DepSumRdr.ERecProxy_True, DepSumRdr.stop_weight * CoeffWeight);

      SliceETrue_True[xb - 1]->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);

      SliceEMu_True[xb - 1]->Fill(DepSumRdr.PrimaryLep_4mom[3], DepSumRdr.stop_weight);
    }
  }


  // Fill selected histos
  EventRates_Selected->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  ERec_Selected->Fill(DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                          DepSumRdr.TotalNonlep_Dep_veto,
                      DepSumRdr.stop_weight);
  EHadr_Selected->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                       DepSumRdr.stop_weight);
  ETrue_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);

  EvRateERec_Selected->Fill(DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                                DepSumRdr.TotalNonlep_Dep_veto,
                            -1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EvRateEHadr_Selected->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                             -1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EvRateEProxy_Selected->Fill(DepSumRdr.ERecProxy_True, -1 * DepSumRdr.vtx[0],
                              DepSumRdr.stop_weight);
  StopEProxy_Selected[DepSumRdr.stop]->Fill(DepSumRdr.ERecProxy_True, DepSumRdr.stop_weight);
  StopEventRates_Selected[DepSumRdr.stop]->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  StopETrue_Selected[DepSumRdr.stop]->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);
  StopERec_Selected[DepSumRdr.stop]->Fill(DepSumRdr.PrimaryLep_4mom[3] +
                                        DepSumRdr.TotalNonlep_Dep_FV +
                                        DepSumRdr.TotalNonlep_Dep_veto,
                                    DepSumRdr.stop_weight);

  EProxyERec_Selected->Fill(DepSumRdr.ERecProxy_True,
                            DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                                DepSumRdr.TotalNonlep_Dep_veto,
                            DepSumRdr.stop_weight);
  ETrueEProxy_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.ERecProxy_True,
                             DepSumRdr.stop_weight);
  ETrueERec_Selected->Fill(DepSumRdr.nu_4mom[3],
                           DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                               DepSumRdr.TotalNonlep_Dep_veto,
                           DepSumRdr.stop_weight);

  ETrueEmu_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.PrimaryLep_4mom[3],
                          DepSumRdr.stop_weight);
  ETrueERecNeutron_Selected->Fill(DepSumRdr.nu_4mom[3],
                                  DepSumRdr.NeutronDep_FV + DepSumRdr.NeutronDep_veto,
                                  DepSumRdr.stop_weight);
  ETrueERecNeutron_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.NeutronDep_timesep_FV + DepSumRdr.NeutronDep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecProton_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.ProtonDep_FV + DepSumRdr.ProtonDep_veto, DepSumRdr.stop_weight);
  ETrueERecProton_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.ProtonDep_timesep_FV + DepSumRdr.ProtonDep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecPion_Selected->Fill(
      DepSumRdr.nu_4mom[3],
      DepSumRdr.PiCDep_FV + DepSumRdr.PiCDep_veto + DepSumRdr.Pi0Dep_FV + DepSumRdr.Pi0Dep_veto,
      DepSumRdr.stop_weight);
  ETrueERecPion_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.PiCDep_timesep_FV + DepSumRdr.PiCDep_timesep_veto +
                          DepSumRdr.Pi0Dep_timesep_FV + DepSumRdr.Pi0Dep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecOther_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.OtherDep_FV + DepSumRdr.OtherDep_veto, DepSumRdr.stop_weight);
  ETrueERecOther_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.OtherDep_timesep_FV + DepSumRdr.OtherDep_timesep_veto,
      DepSumRdr.stop_weight);

  if (CoeffWeightingHelper) {
    if (xb && (xb <= Int_t(SliceEProxy_True_LinCombWeighted.size()))) {
      SliceERec_Selected[xb - 1]->Fill(DepSumRdr.PrimaryLep_4mom[3] +
                                           DepSumRdr.TotalNonlep_Dep_FV +
                                           DepSumRdr.TotalNonlep_Dep_veto,
                                       DepSumRdr.stop_weight);
      SliceEMu_Selected[xb - 1]->Fill(DepSumRdr.PrimaryLep_4mom[3],
                                      DepSumRdr.stop_weight);
      SliceEHadr_Selected[xb - 1]->Fill(
          DepSumRdr.TotalNonlep_Dep_FV + DepSumRdr.TotalNonlep_Dep_veto, DepSumRdr.stop_weight);
      SliceETHadr_Selected[xb - 1]->Fill(
          DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3], DepSumRdr.stop_weight);
    }
  }


  // Fill selected histos
  EventRates_Selected->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  ERec_Selected->Fill(DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                          DepSumRdr.TotalNonlep_Dep_veto,
                      DepSumRdr.stop_weight);
  EHadr_Selected->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                       DepSumRdr.stop_weight);
  ETrue_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);

  EvRateERec_Selected->Fill(DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                                DepSumRdr.TotalNonlep_Dep_veto,
                            -1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EvRateEHadr_Selected->Fill(DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3],
                             -1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  EvRateEProxy_Selected->Fill(DepSumRdr.ERecProxy_True, -1 * DepSumRdr.vtx[0],
                              DepSumRdr.stop_weight);
  StopEProxy_Selected[DepSumRdr.stop]->Fill(DepSumRdr.ERecProxy_True, DepSumRdr.stop_weight);
  StopEventRates_Selected[DepSumRdr.stop]->Fill(-1 * DepSumRdr.vtx[0], DepSumRdr.stop_weight);
  StopETrue_Selected[DepSumRdr.stop]->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.stop_weight);
  StopERec_Selected[DepSumRdr.stop]->Fill(DepSumRdr.PrimaryLep_4mom[3] +
                                        DepSumRdr.TotalNonlep_Dep_FV +
                                        DepSumRdr.TotalNonlep_Dep_veto,
                                    DepSumRdr.stop_weight);

  EProxyERec_Selected->Fill(DepSumRdr.ERecProxy_True,
                            DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                                DepSumRdr.TotalNonlep_Dep_veto,
                            DepSumRdr.stop_weight);
  ETrueEProxy_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.ERecProxy_True,
                             DepSumRdr.stop_weight);
  ETrueERec_Selected->Fill(DepSumRdr.nu_4mom[3],
                           DepSumRdr.PrimaryLep_4mom[3] + DepSumRdr.TotalNonlep_Dep_FV +
                               DepSumRdr.TotalNonlep_Dep_veto,
                           DepSumRdr.stop_weight);

  ETrueEmu_Selected->Fill(DepSumRdr.nu_4mom[3], DepSumRdr.PrimaryLep_4mom[3],
                          DepSumRdr.stop_weight);
  ETrueERecNeutron_Selected->Fill(DepSumRdr.nu_4mom[3],
                                  DepSumRdr.NeutronDep_FV + DepSumRdr.NeutronDep_veto,
                                  DepSumRdr.stop_weight);
  ETrueERecNeutron_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.NeutronDep_timesep_FV + DepSumRdr.NeutronDep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecProton_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.ProtonDep_FV + DepSumRdr.ProtonDep_veto, DepSumRdr.stop_weight);
  ETrueERecProton_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.ProtonDep_timesep_FV + DepSumRdr.ProtonDep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecPion_Selected->Fill(
      DepSumRdr.nu_4mom[3],
      DepSumRdr.PiCDep_FV + DepSumRdr.PiCDep_veto + DepSumRdr.Pi0Dep_FV + DepSumRdr.Pi0Dep_veto,
      DepSumRdr.stop_weight);
  ETrueERecPion_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.PiCDep_timesep_FV + DepSumRdr.PiCDep_timesep_veto +
                          DepSumRdr.Pi0Dep_timesep_FV + DepSumRdr.Pi0Dep_timesep_veto,
      DepSumRdr.stop_weight);
  ETrueERecOther_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.OtherDep_FV + DepSumRdr.OtherDep_veto, DepSumRdr.stop_weight);
  ETrueERecOther_timesep_Selected->Fill(
      DepSumRdr.nu_4mom[3], DepSumRdr.OtherDep_timesep_FV + DepSumRdr.OtherDep_timesep_veto,
      DepSumRdr.stop_weight);

  if (CoeffWeightingHelper) {
    if (xb && (xb <= Int_t(SliceEProxy_True_LinCombWeighted.size()))) {
      SliceERec_Selected[xb - 1]->Fill(DepSumRdr.PrimaryLep_4mom[3] +
                                           DepSumRdr.TotalNonlep_Dep_FV +
                                           DepSumRdr.TotalNonlep_Dep_veto,
                                       DepSumRdr.stop_weight);
      SliceEMu_Selected[xb - 1]->Fill(DepSumRdr.PrimaryLep_4mom[3],
                                      DepSumRdr.stop_weight);
      SliceEHadr_Selected[xb - 1]->Fill(
          DepSumRdr.TotalNonlep_Dep_FV + DepSumRdr.TotalNonlep_Dep_veto, DepSumRdr.stop_weight);
      SliceETHadr_Selected[xb - 1]->Fill(
          DepSumRdr.ERecProxy_True - DepSumRdr.PrimaryLep_4mom[3], DepSumRdr.stop_weight);
    }
  }


  if (CoeffWeightingHelper) {
    TDirectory *oupD = of;

    TDirectory *wD = oupD->mkdir("TruthLinearCombinations");
    wD->cd();

    TH1D *LinCombEProxy = new TH1D(
        "LinCombEProxy", ";E_{#nu,Proxy} (GeV);Events / GeV", 100, 0, 10);
    TH1D *LinCombETrue =
        new TH1D("LinCombETrue", ";E_{#nu} (GeV);Events / GeV", 100, 0, 10);

    for (Int_t xbin = 1; xbin < OffSetVsETrue->GetXaxis()->GetNbins() + 1;
         ++xbin) {
      double bc_etrue = 0, be_etrue = 0, bc_eprox = 0, be_eprox = 0;

      for (Int_t ybin = 1; ybin < OffSetVsETrue->GetYaxis()->GetNbins() + 1;
           ++ybin) {
        bc_etrue += OffSetVsETrue->GetBinContent(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);

        be_etrue += OffSetVsETrue->GetBinError(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);

        bc_eprox += OffSetVsEProxy->GetBinContent(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);
        be_eprox += OffSetVsEProxy->GetBinError(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);
      }

      LinCombETrue->SetBinContent(xbin, bc_etrue);
      LinCombETrue->SetBinError(xbin, be_etrue);

      LinCombEProxy->SetBinContent(xbin, bc_eprox);
      LinCombEProxy->SetBinError(xbin, be_eprox);
    }

    CoeffWeightingHelper->SetDirectory(wD);
    OffSetVsEProxy->SetDirectory(wD);
    OffSetVsETrue->SetDirectory(wD);

    TDirectory *sliceD = oupD->mkdir("SliceObservables");

    for (size_t i = 0; i < SliceEProxy_True_LinCombWeighted.size(); ++i) {
      SliceEProxy_True_LinCombWeighted[i]->SetDirectory(wD);
      SliceETrue_True[i]->SetDirectory(sliceD);
      SliceERec_Selected[i]->SetDirectory(sliceD);
      SliceEMu_True[i]->SetDirectory(sliceD);
      SliceEMu_Selected[i]->SetDirectory(sliceD);
      SliceEHadr_Selected[i]->SetDirectory(sliceD);
      SliceETHadr_Selected[i]->SetDirectory(sliceD);
    }
  }

  TH2D *EProxyERec_Selected_ETrueNorm = static_cast<TH2D *>(
      EProxyERec_Selected->Clone("EProxyERec_SelectedETrueNorm"));
  TH2D *ETrueEProxy_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueEProxy_Selected->Clone("ETrueEProxy_SelectedETrueNorm"));
  TH2D *ETrueERec_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERec_Selected->Clone("ETrueERec_SelectedETrueNorm"));
  TH2D *ETrueEmu_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueEmu_Selected->Clone("ETrueEmu_SelectedETrueNorm"));
  TH2D *ETrueERecNeutron_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecNeutron_Selected->Clone("ETrueERecNeutron_SelectedETrueNorm"));
  TH2D *ETrueERecNeutron_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecNeutron_timesep_Selected->Clone(
          "ETrueERecNeutron_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecProton_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecProton_Selected->Clone("ETrueERecProton_SelectedETrueNorm"));
  TH2D *ETrueERecProton_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecProton_timesep_Selected->Clone(
          "ETrueERecProton_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecPion_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecPion_Selected->Clone("ETrueERecPion_SelectedETrueNorm"));
  TH2D *ETrueERecPion_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecPion_timesep_Selected->Clone(
          "ETrueERecPion_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecOther_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecOther_Selected->Clone("ETrueERecOther_SelectedETrueNorm"));
  TH2D *ETrueERecOther_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecOther_timesep_Selected->Clone(
          "ETrueERecOther_timesep_SelectedETrueNorm"));

  for (Int_t bin_it = 0; bin_it < ETrue_Selected->GetXaxis()->GetNbins();
       ++bin_it) {
    double NETrueSlice = ETrue_Selected->GetBinContent(bin_it + 1);
    if (!NETrueSlice) {
      continue;
    }
    for (Int_t yit = 0;
         yit < ETrueERec_Selected_ETrueNorm->GetYaxis()->GetNbins(); ++yit) {
      EProxyERec_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          EProxyERec_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueEProxy_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueEProxy_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERec_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERec_Selected->GetBinContent(bin_it + 1, yit + 1) / NETrueSlice);
      ETrueEmu_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueEmu_Selected->GetBinContent(bin_it + 1, yit + 1) / NETrueSlice);
      ETrueERecNeutron_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecNeutron_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecNeutron_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1, ETrueERecNeutron_timesep_Selected->GetBinContent(
                                   bin_it + 1, yit + 1) /
                                   NETrueSlice);
      ETrueERecProton_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecProton_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecProton_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecProton_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecPion_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecPion_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecPion_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecPion_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecOther_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecOther_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecOther_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecOther_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
    }
  }

  EProxyERec_Selected_ETrueNorm->SetDirectory(of);
  ETrueEProxy_Selected_ETrueNorm->SetDirectory(of);
  ETrueERec_Selected_ETrueNorm->SetDirectory(of);
  ETrueEmu_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecNeutron_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecNeutron_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecProton_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecProton_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecPion_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecPion_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecOther_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecOther_timesep_Selected_ETrueNorm->SetDirectory(of);


  TH1D *EventRates_Selected = nullptr;
  TH1D *ERec_Selected = nullptr;
  TH1D *EHadr_Selected = nullptr;
  TH2D *EvRateERec_Selected = nullptr;
  TH2D *EvRateEHadr_Selected = nullptr;
  TH2D *EvRateEProxy_Selected = nullptr;
  TH1D *EventRates_Corrected = nullptr;
  TH1D *ERec_Corrected = nullptr;
  TH1D *EHadr_Corrected = nullptr;
  TH2D *EvRateERec_Corrected = nullptr;
  TH2D *EvRateEHadr_Corrected = nullptr;
  TH2D *EvRateEProxy_Corrected = nullptr;

  std::vector<double> XRangeBins;
  if (CorrectInAbsPos) {
    for (Int_t bi_it = 0;
         bi_it < Eproxy_FV_abspos_seleff->GetYaxis()->GetNbins(); ++bi_it) {
      XRangeBins.push_back(
          Eproxy_FV_abspos_seleff->GetYaxis()->GetBinLowEdge(bi_it + 1));
    }
    XRangeBins.push_back(Eproxy_FV_abspos_seleff->GetYaxis()->GetBinUpEdge(
        Eproxy_FV_abspos_seleff->GetYaxis()->GetNbins()));

    EventRates_Selected = new TH1D("EventRates_Selected", ";Offset (cm);Count",
                                   (XRangeBins.size() - 1), XRangeBins.data());
    ERec_Selected =
        new TH1D("ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Selected = new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count",
                              (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Selected =
        new TH2D("EvRateERec_Selected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEHadr_Selected =
        new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEProxy_Selected =
        new TH2D("EvRateEProxy_Selected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EventRates_Corrected =
        new TH1D("EventRates_Corrected", ";Offset (cm);Count",
                 (XRangeBins.size() - 1), XRangeBins.data());
    ERec_Corrected =
        new TH1D("ERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Corrected = new TH1D("EHadr_Corrected", ";E_{Hadr} (GeV);Count",
                               (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Corrected =
        new TH2D("EvRateERec_Corrected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEHadr_Corrected =
        new TH2D("EvRateEHadr_Corrected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEProxy_Corrected =
        new TH2D("EvRateEProxy_Corrected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
  } else {
    EventRates_Selected =
        new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
    ERec_Selected =
        new TH1D("ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Selected = new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count",
                              (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Selected =
        new TH2D("EvRateERec_Selected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEHadr_Selected =
        new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEProxy_Selected =
        new TH2D("EvRateEProxy_Selected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EventRates_Corrected =
        new TH1D("EventRates_Corrected", ";Offset (cm);Count", 400, -250, 3750);
    ERec_Corrected =
        new TH1D("ERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Corrected = new TH1D("EHadr_Corrected", ";E_{Hadr} (GeV);Count",
                               (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Corrected =
        new TH2D("EvRateERec_Corrected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEHadr_Corrected =
        new TH2D("EvRateEHadr_Corrected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEProxy_Corrected =
        new TH2D("EvRateEProxy_Corrected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
  }

  TDirectory *oupD = of;

  TDirectory *sliceD = oupD->mkdir("SliceEventRates");
  sliceD->cd();

  std::vector<TH1D *> SliceERec_Selected;
  std::vector<TH1D *> SliceERec_Corrected;
  TH1D *SliceHelper = nullptr;
  if (CorrectInAbsPos) {
    SliceHelper =
        new TH1D("SliceHelper", "", (XRangeBins.size() - 1), XRangeBins.data());
    SliceHelper->SetDirectory(nullptr);

    for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
      SliceERec_Selected.push_back(new TH1D(
          (std::string("SliceERec_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{Rec} (GeV);Count", 100, 0, 10));
      SliceERec_Corrected.push_back(new TH1D(
          (std::string("SliceERec_Corrected") + to_str(bin_it - 1)).c_str(),
          ";E_{Rec} (GeV);Count", 100, 0, 10));
    }
  }

  TDirectory *stopD = oupD->mkdir("StopEventRates");
  stopD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEventRates_Corrected;

  std::vector<TH1D *> StopEProxy_Selected;
  std::vector<TH1D *> StopEProxy_Corrected;

  std::vector<TH1D *> StopERec_Selected;
  std::vector<TH1D *> StopERec_Corrected;

  for (Int_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_Corrected.push_back(new TH1D(
        (std::string("StopEventRates_Corrected_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));

    StopEProxy_Selected.push_back(new TH1D(
        (std::string("StopEProxy_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEProxy_Corrected.push_back(new TH1D(
        (std::string("StopEProxy_Corrected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));

    StopERec_Selected.push_back(new TH1D(
        (std::string("StopERec_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopERec_Corrected.push_back(new TH1D(
        (std::string("StopERec_Corrected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
  }
  oupD->cd();

  Int_t xb = SliceHelper->GetXaxis()->FindFixBin(-1 * edr.vtx[0]);
  if (xb && (xb <= Int_t(SliceERec_Selected.size()))) {
    SliceERec_Selected[xb - 1]->Fill(edr.PrimaryLep_4mom[3] +
                                         edr.TotalNonlep_Dep_FV +
                                         edr.TotalNonlep_Dep_veto,
                                     edr.stop_weight);
    SliceERec_Corrected[xb - 1]->Fill(edr.PrimaryLep_4mom[3] +
                                          edr.TotalNonlep_Dep_FV +
                                          edr.TotalNonlep_Dep_veto,
                                      edr.stop_weight * effweight);
  }
  }

  // Fill selected histos
  EventRates_Selected->Fill(-1 * edr.vtx[0], edr.stop_weight);
  ERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                        edr.TotalNonlep_Dep_veto,
                    edr.stop_weight);
  EHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                     edr.stop_weight);

  EvRateERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                              edr.TotalNonlep_Dep_veto,
                          -1 * edr.vtx[0], edr.stop_weight);
  EvRateEHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                           -1 * edr.vtx[0], edr.stop_weight);
  StopEventRates_Selected[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);
  EvRateEProxy_Selected->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                            edr.stop_weight);
  StopEProxy_Selected[edr.stop]->Fill(edr.ERecProxy_True, edr.stop_weight);
  StopERec_Selected[edr.stop]->Fill(edr.PrimaryLep_4mom[3] +
                                      edr.TotalNonlep_Dep_FV +
                                      edr.TotalNonlep_Dep_veto,
                                  edr.stop_weight);

  EventRates_Corrected->Fill(-1 * edr.vtx[0], edr.stop_weight * effweight);
  ERec_Corrected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                         edr.TotalNonlep_Dep_veto,
                     edr.stop_weight * effweight);
  EHadr_Corrected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                      edr.stop_weight * effweight);

  EvRateERec_Corrected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                               edr.TotalNonlep_Dep_veto,
                           -1 * edr.vtx[0], edr.stop_weight * effweight);
  EvRateEHadr_Corrected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                            -1 * edr.vtx[0], edr.stop_weight * effweight);
  StopEventRates_Corrected[edr.stop]->Fill(-1 * edr.vtx[0],
                                         edr.stop_weight * effweight);
  EvRateEProxy_Corrected->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                             edr.stop_weight * effweight);
  StopEProxy_Corrected[edr.stop]->Fill(edr.ERecProxy_True,
                                     edr.stop_weight * effweight);
  StopERec_Corrected[edr.stop]->Fill(edr.PrimaryLep_4mom[3] +
                                       edr.TotalNonlep_Dep_FV +
                                       edr.TotalNonlep_Dep_veto,
                                   edr.stop_weight * effweight);
