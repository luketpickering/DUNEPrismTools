{

  char buff[1000];

  TH2D * hDunePRISM = new TH2D("oa_enu_numode_numu","oa_enu_numode_numu", 1000, 0, 10, 130, 0, 6.45);
  
  for (double oa = 0; oa < 113; oa += 0.872665){

    double oaDegree = oa*180/3.14159265/1000.;
    
    sprintf(buff, "/mnt/home/f0003917/DunePRISM/offAxis/histos_OffAxis_%f_mrad.root", oa);
    std::cout << "Opening " << buff << std::endl;
    
    TFile * fIn = new TFile (buff);
    TH1D * slice = (TH1D*) fIn->Get("numu_flux_forplots");

    yBin = hDunePRISM->GetYaxis()->FindBin(oaDegree);
    for (int i = 1; i < slice->GetNbinsX(); i++){
      double xBinCenter = slice->GetBinCenter(i);
      int xBin = hDunePRISM->GetXaxis()->FindBin(xBinCenter);

      hDunePRISM->SetBinContent(xBin, yBin, slice->GetBinContent(i));
      hDunePRISM->SetBinError(xBin, yBin, slice->GetBinError(i));
    }
    delete slice;
    fIn->Close();
    delete fIn;
  }

  hDunePRISM->Draw("COLZ");

  TFile * outFile = new TFile("duneprism_spectra.root", "RECREATE");
  hDunePRISM->Write();
  outFile->Close();

}
