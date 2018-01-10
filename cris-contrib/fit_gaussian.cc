#include "stdlib.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TLegend.h"

// TH1D *sk_flux;
TH2D *oa_flux; //, *oa_flux_numub, *oa_flux_nue;
double nuprism_flux[1000];
double nuprism_error[1000];
double Elo,Ehi,enu, sigma, fluxNorm=1.;
TF1* gaussian;

//const int nBins = 60; // NuPRISM 
//const int binOffset = 22; // NuPRISM: Start at 1 degree OA

//const int NParams = 120; // DunePRISM: Go all the way to 6 degrees (is this sensible?)
//const int NParams = 70; // DunePRISM: Go all the way to 3.5 degrees: equivalent to 35 m wide detector
const int NParams = 66; // DunePRISM: Go all the way to 3.3 degrees: equivalent to 33 m wide detector
const int fOffset = 1; // DunePRISM: want to go down to 0 degrees


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   double chi2 = 0.;
    
   for(int j=0; j<1000; j++){
     nuprism_flux[j] = 0.;
     nuprism_error[j] = 0.;
   }

   //Make the predicted flux
   for(int i=0; i<NParams; i++)
     for(int j=0; j<500; j++){
       nuprism_flux[j] += par[i]*oa_flux->GetBinContent(j+1,i+fOffset);
       nuprism_error[j] += pow(par[i]*oa_flux->GetBinError(j+1,i+fOffset),2);
     }

   for(int i=0; i<NParams-1; i++) chi2 += pow((par[i]-par[i+1])/0.0002,2);
   //for(int i=0; i<60; i++) chi2 += pow((par[i]-0.0)/0.05,2);

   double max = 0.0;
   double maxe = 0.0;
   for(int j=0; j<500; j++){
     double enu = 0.01+(double)j*0.02;
     if(nuprism_flux[j]<0.) chi2 += pow((0.0-nuprism_flux[j])/(fluxNorm*0.00001),2);
     if(enu<Elo || enu>Ehi) continue;
     //Calculate chi2 with 2% error
     //600 MeV
     chi2 += pow(nuprism_flux[j]-gaussian->Eval(enu),2)/(pow(0.00005*fluxNorm,2)+pow(0.0001*gaussian->Eval(enu),2));
     //800 MeV
     //chi2 += pow(nuprism_flux[j]-gaussian->Eval(enu),2)/(pow(0.0004*fluxNorm*enu,2)+pow(0.0002*gaussian->Eval(enu),2));
     if(nuprism_flux[j]>max){
       max = nuprism_flux[j];
       maxe = enu;
     }  
   }


   f = chi2;//+pow((maxe-gaussian->GetParameter(1))/0.0001,2);
}
     

void fit_gaussian(double mean, double width, double elo=0.2, double ehi=9.0)
{

  Elo=elo;
  Ehi=ehi;

  enu=mean;
  sigma=width;

  gaussian = new TF1("gaussian","[0]*TMath::Exp(-0.5*pow(TMath::Abs((x-[1])/[2]),[3]))",0.4,2);
  gaussian->SetParameter(1,mean);
  gaussian->SetParameter(2,width);
  gaussian->SetParameter(3,1.7);

  //Read in the fluxes
  TFile *fi = new TFile("duneprism_spectra.root");

  
  
  oa_flux = (TH2D*)fi->Get("oa_enu_numode_numu");
  //  oa_flux_numub = (TH2D*)fi->Get("oa_enu_numode_numub");
  //  oa_flux_nue = (TH2D*)fi->Get("oa_enu_numode_nue");

  //  TH1D *sk_flux = (TH1D*)fi->Get("sk_nom_numode_numu");


  //Rebin by a factor of 2 so there are 0.1 degree bins
  //  oa_flux->RebinY(2);
  //oa_flux_numub->RebinY(2);
  //oa_flux_nue->RebinY(2);
  oa_flux->RebinX(2);
  //  oa_flux_numub->RebinX(2);
  //  oa_flux_nue->RebinX(2);

  for(int i=1; i<=oa_flux->GetNbinsX(); i++){
    double enu = oa_flux->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=oa_flux->GetNbinsY(); j++){
      oa_flux->SetBinContent(i,j,oa_flux->GetBinContent(i,j)*enu);
      //      oa_flux_numub->SetBinContent(i,j,oa_flux_numub->GetBinContent(i,j)*enu);
      //      oa_flux_nue->SetBinContent(i,j,oa_flux_nue->GetBinContent(i,j)*enu);
    }
  }
  
      
  //TH1D *oa_flux_comp = (TH1D*)oa_flux->ProjectionX("oa_flux_comp",34,35);
  // std::cout << "Off-axis comp: " <<  0.5*(oa_flux->GetYaxis()->GetBinCenter(34)+oa_flux->GetYaxis()->GetBinCenter(35)) << std::endl;
  //  double oaa =  0.5*(oa_flux->GetYaxis()->GetBinCenter(34)+oa_flux->GetYaxis()->GetBinCenter(35));
  //800 MeV
  //  TH1D *oa_flux_comp = (TH1D*)oa_flux->ProjectionX("oa_flux_comp",34,35);
  //  oa_flux_comp->Rebin(4);
  //  oa_flux_comp->Scale(0.25);


  //Divide by the area so the histogram is properly fluxNormalized
  oa_flux->Scale(1.0/1200.0/174.533);
  //  oa_flux_comp->Scale(1.0/1200.0/174.533);
  //  oa_flux_numub->Scale(1.0/1200.0/174.533);
  //  oa_flux_nue->Scale(1.0/1200.0/174.533);

  fluxNorm = oa_flux->GetMaximum();
  gaussian->SetParameter(0,oa_flux->GetMaximum());

  double bin_edges[NParams+1];
  for(int i=0; i<NParams+1; i++) bin_edges[i] = oa_flux->GetYaxis()->GetBinLowEdge(i+fOffset);


  //Set up the fitter
  TMinuit *gMinuit = new TMinuit(NParams); 
  gMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  //Set the parameter initial values
  for(int i=0; i<NParams; i++)
     gMinuit->mnparm(i, Form("oa%d",i), 1.0/30.0, 0.001, -10.0,10.0,ierflg);

  //Do the minimization
  arglist[0] = 500000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  std::cout << "Test 1" << std::endl;

  double par[NParams],perror[NParams]; 
   
  std::cout << "Test 2" << std::endl;

  //Graphs and histograms for the fit results
  TGraph *coeff = new TGraph();
  TH1D *fitted_flux = new TH1D("fitted_flux","",500,0,10.0);
 
  std::cout << "Test 3" << std::endl;

  TH1D *coefficients = new TH1D("Coefficient","",NParams,bin_edges);

  std::cout << "Test 4" << std::endl;

  //Iterate through the coefficients and fill this graphs and histograms with the fit results
  for(int i=0; i<NParams; i++){
    gMinuit->GetParameter(i,par[i],perror[i]);
    coefficients->SetBinContent(i+1,par[i]);
    coeff->SetPoint(i,(double)i*0.1+0.95,par[i]);
    for(int j=0; j<500; j++){
      fitted_flux->SetBinContent(j+1,fitted_flux->GetBinContent(j+1)+par[i]*oa_flux->GetBinContent(j+1,i+fOffset));
     }
  }

  //sk_flux->Rebin(2); 
  //sk_flux->Scale(fitted_flux->GetMaximum()/sk_flux->GetMaximum());
  //oa_flux_comp->Scale(fitted_flux->GetMaximum()/oa_flux_comp->GetMaximum());

  double max = fitted_flux->GetMaximum();
  for(int i=1; i<=fitted_flux->GetNbinsX(); i++) fitted_flux->SetBinError(i,max*0.05);

  gaussian->SetParameter(3,2.0);
  //gaussian->FixParameter(2,0.09);
  gaussian->SetParameter(2,0.1);
  gaussian->SetParLimits(2,0.01,1.0);
  //gaussian->FixParameter(1,0.8);
  gaussian->FixParameter(3,2.0);
  fitted_flux->Fit("gaussian","R");
  //gaussian->SetParameter(1,1.2);


  std::cout << "Test 5" << std::endl;
  //Draw the numu flux comparison
  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  gPad->SetRightMargin(0.05);
  fitted_flux->Rebin(2);
  fitted_flux->Scale(0.5);
  fitted_flux->SetLineWidth(2);
  fitted_flux->SetLineColor(2);
  fitted_flux->SetStats(false);
  fitted_flux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  fitted_flux->GetYaxis()->SetTitle("Arb. norm.");
  fitted_flux->GetXaxis()->SetRangeUser(0.0,10.0);
  fitted_flux->Draw("hist"); 
  
  gaussian->SetNpx(500);
  gaussian->SetLineWidth(1);
  gaussian->Draw("same");
  //  oa_flux_comp->Scale(fitted_flux->GetMaximum()/oa_flux_comp->GetMaximum());
  //  oa_flux_comp->SetLineColor(4);
  //  oa_flux_comp->Draw("hist same");
  TLegend *leg = new TLegend(0.5,0.5,0.8,0.8);
  leg->AddEntry(fitted_flux,"Linear Combination","l");
  //  leg->AddEntry(oa_flux_comp,Form("%0.1f#circ Off-axis Flux",oaa),"l");
  leg->AddEntry(gaussian,Form("Gaussian: Mean=%0.2f, RMS=%0.2f GeV",gaussian->GetParameter(1),gaussian->GetParameter(2)),"l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0.);
  leg->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",600,500);
  gPad->SetRightMargin(0.05);
  coefficients->SetMarkerStyle(20);
  coefficients->SetStats(false);
  coefficients->GetXaxis()->SetTitle("#theta_{OA} (#circ)");
  coefficients->GetYaxis()->SetTitle("Coefficient Value");
  coefficients->Draw("p");

  TFile fout("nuprism_coef.root","RECREATE");
  coefficients->Write("Coefficients");
  fitted_flux->Write("FittedFlux");
  fout.Close();


  return;


}
