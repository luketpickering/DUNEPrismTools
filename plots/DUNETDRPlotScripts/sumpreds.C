void sumpreds() {
  TH1 *summer = nullptr;
  TCanvas c1("c1","");
  int num = 0;
  double max = 0;
  TFile *fin = new TFile("allfdnu_numu.root");
  for (size_t i = 0; i < 500; ++i) {
    std::stringstream n;
    n << "LBNF_numu_flux_Nom_22156546_" << i;

    TH1 *h = dynamic_cast<TH1 *>(fin->Get(n.str().c_str()));
    if (!h) {
      continue;
    }

    max = std::max(max,h->GetMaximum());

    if (!summer) {
      summer = dynamic_cast<TH1 *>(h->Clone("summer"));
      summer->SetDirectory(nullptr);
      h->GetYaxis()->SetRangeUser(0,7.47783e-14);
      h->Draw("HIST");
    } else {
      summer->Add(h);
      h->Draw("HISTSAME");
    }
    num++;
  }

  std::cout << "max: " << max << std::endl;

  summer->Scale(1.0/double(num));
  summer->SetLineColor(kRed);
  summer->SetLineWidth(4);
  summer->Draw("HISTSAME");

  c1.SaveAs("test.pdf");

  TFile *fout = new TFile("fdnumu_summer.root", "RECREATE");
  fout->WriteTObject(summer, "summer");
  fout->Close();
}
