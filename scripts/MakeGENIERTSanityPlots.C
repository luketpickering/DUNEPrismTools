{
  TTree *gRooTracker = dynamic_cast<TTree*>(_file0->Get("gRooTracker"));

  TCanvas *c1 = new TCanvas("c1");

  gRooTracker->Draw("StdHepP4[0][1]:StdHepP4[0][0]");
  c1->SaveAs("GENIEFlux.pdf[");
  c1->SaveAs("GENIEFlux.pdf");
  gRooTracker->Draw("StdHepP4[0][1]:StdHepP4[0][2]");
  c1->SaveAs("GENIEFlux.pdf");
  gRooTracker->Draw("StdHepP4[0][3]");
  c1->SaveAs("GENIEFlux.pdf");
  gRooTracker->Draw("EvtVtx[1]:EvtVtx[0]");
  c1->SaveAs("GENIEFlux.pdf");
  gRooTracker->Draw("EvtVtx[1]:EvtVtx[2]");
  c1->SaveAs("GENIEFlux.pdf");
  c1->SaveAs("GENIEFlux.pdf]");
}

