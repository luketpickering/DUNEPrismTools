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
#include "TSystem.h"
#include "TApplication.h"

//#include "../VectorAnalysis/include/probNuPrism.h"

double flux_ratio;
TH1D *sk_flux; //, *sk_flux_numub, *sk_flux_nue;
TH2D *oa_flux; //, *oa_flux_numub, *oa_flux_nue;
double nuprism_flux[2000];
double nuprism_error[2000];
double Elo,Ehi;

// The number of nuPRISM slices to use in the fit
// const int NParams = 120; // Up to 6 degrees
const int NParams = 66; // Only up to 3.3 degrees / 33 m
const int fOffset = 1; // All the way down to 0 degrees

// There are 30 slices to span 3 degrees, so each slice is 0.1 degrees
// The off-axis angle flux histogram is actually binned in 0.05 degree bins,
// so in the fit fcn we loop over NParams*2 in steps of 2, adding up the number
// of entries from bins i and i+1.
//
// There is also a 22 bin offset applied.  This is because the flux histogram
// starts at -0.1 degrees off-axis, and so we need to get to 1 degree off axis, which is 22 bins

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double chi2 = 0.;

    for(int j=0; j<oa_flux->GetNbinsX(); j++){
        nuprism_flux[j] = 0.;
        nuprism_error[j] = 0.;
    }

    //Make the predicted flux
    for(int i=0; i<NParams*1; i+=1)
        for(int j=0; j<oa_flux->GetNbinsX(); j++){
            nuprism_flux[j] += flux_ratio*par[i/1]*(oa_flux->GetBinContent(j+1,i+fOffset)/*+oa_flux->GetBinContent(j+1,(i+1)+fOffset)*/);
	    nuprism_error[j] += pow(flux_ratio*par[i/1]*(oa_flux->GetBinContent(j+1,i+fOffset)/*+oa_flux->GetBinContent(j+1,(i+1)+fOffset)*/),2);
	    //	    nuprism_error[j] += pow(flux_ratio*par[i/1]*oa_flux->GetBinError(j+1,i+fOffset),2);
        }

    // Apply a strong penalty term if neighbouring parameters are not close to one another
    // This creates the smooth parameter distribution that reduces the statistical uncertainty
    for(int i=0; i<NParams; i++){
//        chi2 += pow(fabs(par[i] - 0.1)/2.,2);
      if(par[i] != 0 && par[i+1] != 0) chi2 += pow((par[i]-par[i+1])/(0.01),2);
    }

    for(int j=0; j<oa_flux->GetNbinsX(); j++){
        double enu = sk_flux->GetXaxis()->GetBinCenter(j+1);
        if(enu<Elo || enu>Ehi) continue;
        //Calculate chi2 with 2% error
        // chi2 += pow(nuprism_flux[j]-sk_flux->GetBinContent(j+1),2)/(pow(0.01*sk_flux->GetBinContent(j+1),2)+pow(sk_flux->GetBinError(j+1),2));
        chi2 += pow(nuprism_flux[j]-sk_flux->GetBinContent(j+1),2)/(pow(0.01*sk_flux->GetBinContent(j+1),2));


	//	std::cout << "bin="<<j<< " np_flux=" << nuprism_flux[j] << " sk_flux=" << sk_flux->GetBinContent(j+1) << std::endl;
    }

    f = chi2;
}


//int fit_spectrum(double theta, double dm2, bool appearance=false, double elo=0.2, double ehi=9.0)
int fit_spectrum(double elo=.65, double ehi=5., string oscFluxfile = "")
{

    // Set the low and high energy limits of the fit
    Elo=elo;
    Ehi=ehi;

    //    std::cout << "Fit for Theta 23 = " << theta << " and  Dm2 23 = " << dm2 << std::endl;

    //Approximate N/F ratio
    //    flux_ratio = pow(1.0/295.0,2);
    flux_ratio = pow(574/1300000., 2);
    // A more accurate ratio
    //flux_ratio = (0.96221*0.96221)/(295.016*295.016);

    //Read in the fluxes
    TFile *finup = new TFile("duneprism_spectra.root");

    //    TH1D * sk_flux = (TH1D*)fDuneFlux->Get("numu_fluxosc_forplots");
    //    sk_flux_numub = (TH1D*)finup->Get("sk_nom_numode_numub");
    //    sk_flux_nue = (TH1D*)finup->Get("sk_nom_numode_nue");

    //    sk_flux->Sumw2();
    //    sk_flux_numub->Sumw2();
    //    sk_flux_nue->Sumw2();

   
    //Rebin in X to make wider momentum bins (10MeV -> 5MeV)
    //  sk_flux->RebinX(10);
//    sk_flux_numub->RebinX(10);
//    sk_flux_nue->RebinX(10);

//    ProbNuPrism* prob = new ProbNuPrism(dm2,theta);

    //Apply the oscillation probability and scale by the bin energy to approximate the effects of the cross section
    //    for(int i=1; i<=sk_flux->GetNbinsX(); i++){

      //        double xsecweight_numu = 1;//sk_flux->GetXaxis()->GetBinCenter(i);
	//        double xsecweight_nue = 1;//sk_flux_nue->GetXaxis()->GetBinCenter(i);
	//        double xsecweight_numub = 1;//sk_flux_numub->GetXaxis()->GetBinCenter(i);

	//        sk_flux->SetBinContent(i,sk_flux->GetBinContent(i)*prob->GetProbNuMuNuMu(sk_flux->GetXaxis()->GetBinCenter(i))*xsecweight_numu);
	//        sk_flux_numub->SetBinContent(i,sk_flux_numub->GetBinContent(i)*prob->GetProbNuMuBarNuMuBar(sk_flux_numub->GetXaxis()->GetBinCenter(i))*xsecweight_numub);
	//        sk_flux_nue->SetBinContent(i,sk_flux_nue->GetBinContent(i)*prob->GetProbNuENuE(sk_flux_nue->GetXaxis()->GetBinCenter(i))*xsecweight_nue);
    //    }


    //    oa_flux = (TH2D*)finup->Get("oa_enu_numode_numu");
    //    oa_flux->RebinY(2);

    //    sk_flux->RebinX(4);
    //    oa_flux->RebinX(4);

    // Binning scheme for neutrino energy
    const int nXbins = 213;
    double xBins[nXbins+1];
    int binStops[] = {0, 75, 138, nXbins};
    float binWidths[] = {0.02, 0.04, 0.08};
    int binCounter = 0;
    for (int j = 0; j < sizeof(binWidths)/sizeof(binWidths[0]); j++){
      for (int i = binStops[j]; i <= binStops[j+1]; i++){
        if (!binCounter) xBins[binCounter] = 0.;
        else xBins[binCounter] = xBins[binCounter-1] + binWidths[j];
        std::cout << binCounter << " " << xBins[binCounter] << std::endl;
        binCounter++;
      }
    }
    
    
    //    sk_flux = (TH1D*) sk_flux->Rebin(nXbins, "osc_flux_rebin", xBins);
    //    TCanvas * ctemp1 = new TCanvas ("ctemp1", "ctemp1", 800, 600);
    //    sk_flux->Draw("");
    // ################################################################################
    // WARNING!!! REBINNING THIS WAY DOES NOT PROPAGATE BIN ERRORS. THIS IS ONLY VALID IF GETBINERROR NOT USED!!!
    // ################################################################################
    TFile * fDuneFlux;
    if (! oscFluxfile.size() ){
    //The SK fluxes
    //    TFile * fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_A.root"); Elo = 0.7; 
    //    TFile * fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_B.root"); Elo = 0.85;
    //    TFile * fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_C.root"); Elo = 0.6;
    //    TFile * fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_D.root"); Elo = 0.5;
    //    TFile * fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_E.root"); Elo = 0.95;
      fDuneFlux = new TFile("/home/cvilela/DunePRISM/offAxis/OscFlux/histos_OffAxis_0.0_mrad_F.root"); Elo = 0.5;
    } else {
      fDuneFlux = new TFile(oscFluxfile.c_str());
    }
    

    TH1D * sk_flux_fixed_binning = (TH1D*) fDuneFlux->Get("numu_fluxosc_forplots");
    
    TAxis *xaxis_sk = sk_flux_fixed_binning->GetXaxis();
    sk_flux = new TH1D("OscillatedFlux", "Oscillated Flux;E_{#nu} [GeV];Arbitrary", nXbins, xBins);
    for (int i = 1; i <= xaxis_sk->GetNbins(); i++){
      //      sk_flux->Fill(xaxis_sk->GetBinCenter(i), sk_flux_fixed_binning->GetBinContent(i));
      sk_flux->Fill(xaxis_sk->GetBinCenter(i),
                    sk_flux_fixed_binning->GetBinContent(i)*
                    xaxis_sk->GetBinWidth(i)/sk_flux->GetBinWidth(sk_flux->FindBin(xaxis_sk->GetBinCenter(i)))
                    );
    }

    TCanvas * ctemp1 = new TCanvas("ctemp1", "ctemp1", 800, 600);
    sk_flux->Draw("HIST");
    
    TH2D * oa_flux_fixed_binning =  (TH2D*)finup->Get("oa_enu_numode_numu");
    TAxis *xaxis = oa_flux_fixed_binning->GetXaxis();
    TAxis *yaxis = oa_flux_fixed_binning->GetYaxis();
    
    oa_flux = new TH2D ("oa_flux", "oa_flux", nXbins, xBins, yaxis->GetNbins(), yaxis->GetXmin(), yaxis->GetXmax());
    for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
        oa_flux->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),
                      oa_flux_fixed_binning->GetBinContent(i,j) *
                      xaxis->GetBinWidth(i)/oa_flux->GetXaxis()->GetBinWidth(oa_flux->GetXaxis()->FindBin(xaxis->GetBinCenter(i)))
                      );
      }
    }

    TCanvas * ctemp = new TCanvas("ctemp", "ctemp", 800, 600);
    oa_flux->Draw("COLZ");
    
//    oa_flux->GetXaxis()->Rebin(250, "oa_flux_rebin", xBins);
    
    //    oa_flux_numub = (TH2D*)finup->Get("oa_enu_numode_numub");
    //    oa_flux_nue = (TH2D*)finup->Get("oa_enu_numode_nue");

    //Rebin by a factor of 2 so there are 0.1 degree bins
//    oa_flux->RebinX(10);
//    oa_flux_numub->RebinX(10);
//    oa_flux_nue->RebinX(10);

    //Scale the off-axis fluxes by the energy
    for(int i=1; i<=oa_flux->GetNbinsX(); i++){
        for(int j=1; j<=oa_flux->GetNbinsY(); j++){
            
            double xsecweight_numu = 1;//oa_flux->GetXaxis()->GetBinCenter(i);
	    //            double xsecweight_nue = 1;//oa_flux_nue->GetXaxis()->GetBinCenter(i);
	    //            double xsecweight_numub = 1;//oa_flux_numub->GetXaxis()->GetBinCenter(i);

            oa_flux->SetBinContent(i,j, oa_flux->GetBinContent(i,j)*xsecweight_numu);
	    //            oa_flux_numub->SetBinContent(i,j, oa_flux_numub->GetBinContent(i,j)*xsecweight_numub);
	    //            oa_flux_nue->SetBinContent(i,j, oa_flux_nue->GetBinContent(i,j)*xsecweight_nue);
        }
    }

    //Divide by the area so the histogram is properly normalized
    //    oa_flux->Scale(1./1200.0/174.533);
    //    oa_flux_numub->Scale(1./1200.0/174.533);
    //    oa_flux_nue->Scale(1./1200.0/174.533);

    oa_flux->Scale(1e22/pow(1300000,2)); // POT?
    sk_flux->Scale(1e22/pow(1300000,2)); // POT?
    
    double bin_edges[NParams+1];
    for(int i=0; i<NParams+1; i++){
        //Print out the off-axis angle bin edges to check the calculations were correct
        std::cout << "oa_flux->GetYaxis()->GetBinLowEdge() = " << oa_flux->GetYaxis()->GetBinLowEdge(1*i+fOffset) << std::endl;
        bin_edges[i] = oa_flux->GetYaxis()->GetBinLowEdge(1*i+fOffset);
    }

    //Set up the fitter
    TMinuit *gMinuit = new TMinuit(NParams); 
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    //Set the parameter initial values
    for(int i=0; i<NParams; i++)
      gMinuit->mnparm(i, Form("oa%d",i), 1.0/double(NParams), 0.001, -10.0,10.0,ierflg);

    // Make gaps
    //    for (int i = 5; i < 15; i++){// 2.5 to 7.5 m
    //      gMinuit->mnparm(i, Form("oa%d",i), 0., 0.001, -10.0,10.0,ierflg);
    //      gMinuit->FixParameter(i);
    //    }
    //    for (int i = 25; i < 35; i++){// 12.5 to 17.5 m
    //      gMinuit->mnparm(i, Form("oa%d",i), 0., 0.001, -10.0,10.0,ierflg);
    //      gMinuit->FixParameter(i);
    //    }
    
    
    //Do the minimization
    arglist[0] = 500000;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    //    double par[100],perror[100];
    double par[NParams],perror[NParams]; 

    //Graphs and histograms for the fit results
    TGraph *coeff = new TGraph();
    TH1D *fitted_flux = (TH1D*)sk_flux->Clone("fitted_flux");
    fitted_flux->Reset();
    //  fitted_flux->SetBins(11, cov_energy_bins);
    //    TH1D *fitted_flux_numub = (TH1D*)sk_flux->Clone("fitted_flux_numub");
    //    fitted_flux_numub->Reset();
    //  fitted_flux_numub->SetBins(11, cov_energy_bins);
    //    TH1D *fitted_flux_nue = (TH1D*)sk_flux->Clone("fitted_flux_nue");
    //    fitted_flux_nue->Reset();
    //  fitted_flux_nue->SetBins(11, cov_energy_bins);
    TH1D *flux_dev = (TH1D*)sk_flux->Clone("flux_dev");
    flux_dev->Reset();
    //  flux_dev->SetBins(11, cov_energy_bins);
    //    TH1D *flux_dev_numub = (TH1D*)sk_flux->Clone("flux_dev_numub");
    //    flux_dev_numub->Reset();
    //  flux_dev_numub->SetBins(11, cov_energy_bins);
    //    TH1D *flux_dev_nue = (TH1D*)sk_flux->Clone("flux_dev_nue");
    //    flux_dev_nue->Reset();
    //  flux_dev_nue->SetBins(11, cov_energy_bins);

    TH1D *coefficients = new TH1D("Coefficient","",NParams,bin_edges);

    //Iterate through the coefficients and fill this graphs and histograms with the fit results
    //Again, the binning requires some odd indices
    for(int i=0; i<NParams*1; i+=1){
        gMinuit->GetParameter(i/1,par[i/1],perror[i/1]);

	coefficients->SetBinContent(i/1+1,par[i/1]);
        coeff->SetPoint(i,(double)(i+fOffset)*0.05 - 0.1,par[i]);
        for(int j=0; j<oa_flux->GetNbinsX(); j++){
            fitted_flux->SetBinContent(j+1,fitted_flux->GetBinContent(j+1)+flux_ratio*par[i/1]*(oa_flux->GetBinContent(j+1,i+fOffset)/* + oa_flux->GetBinContent(j+1,i+1+fOffset)*/));
	    //            fitted_flux_numub->SetBinContent(j+1,fitted_flux_numub->GetBinContent(j+1)+flux_ratio*par[i/1]*(oa_flux->GetBinContent(j+1,i+fOffset)/* + oa_flux->GetBinContent(j+1,i+1+fOffset)*/));
	    //            fitted_flux_nue->SetBinContent(j+1,fitted_flux_nue->GetBinContent(j+1)+flux_ratio*par[i/1]*(oa_flux->GetBinContent(j+1,i+fOffset)/* + oa_flux->GetBinContent(j+1,i+1+fOffset)*/));
	    flux_dev->SetBinContent(j+1,(fitted_flux->GetBinContent(j+1)-sk_flux->GetBinContent(j+1))/sk_flux->GetBinContent(j+1));
	    //            flux_dev_numub->SetBinContent(j+1,(fitted_flux_numub->GetBinContent(j+1)-sk_flux_numub->GetBinContent(j+1))/sk_flux_numub->GetBinContent(j+1));
	    // flux_dev_nue->SetBinContent(j+1,(fitted_flux_nue->GetBinContent(j+1)-sk_flux_nue->GetBinContent(j+1))/sk_flux_nue->GetBinContent(j+1));
        }
    }

    //Draw the numu flux comparison
    TCanvas *c1 = new TCanvas("c1","c1",600,500);
    fitted_flux->SetLineWidth(2);
    fitted_flux->SetLineColor(2);
    sk_flux->SetLineWidth(2);
    sk_flux->SetStats(false);
    sk_flux->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    sk_flux->GetYaxis()->SetTitle("Arbitrary");
    sk_flux->Draw("hist");
    fitted_flux->Draw("same hist"); 

    //Draw the numu-bar flux comparison
//    TCanvas *c4 = new TCanvas("c4","c4",600,500);
//    fitted_flux_numub->SetLineWidth(2);
//    fitted_flux_numub->SetLineColor(2);
//    sk_flux_numub->SetLineWidth(2);
//    sk_flux_numub->SetStats(false);
//    sk_flux_numub->GetXaxis()->SetTitle("E_{#nu} (GeV)");
//    sk_flux_numub->GetYaxis()->SetTitle("Flux/[cm^{2}#upoint 100 MEV #upoint 1e21 POT]");
//    sk_flux_numub->Draw();
//    fitted_flux_numub->Draw("same hist"); 

    //Draw the coefficents from the fit 
//    TCanvas *c2 = new TCanvas("c2","c2",600,500);
//    coeff->SetMarkerStyle(20);
//    coeff->SetMarkerSize(1.0);
//    coeff->GetXaxis()->SetTitle("Off-axis Angle (degrees)");
//    coeff->GetYaxis()->SetTitle("Fitted Coefficient");
//    coeff->Draw("AP");

    //Draw the deviation of the numu flux
    TCanvas *c3 = new TCanvas("c3","c3",600,500);
    flux_dev->SetStats(false);
    flux_dev->SetLineWidth(2);
    flux_dev->SetLineColor(2);
    flux_dev->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    flux_dev->GetYaxis()->SetTitle("Fractional Deviation");
    flux_dev->Draw("hist");

    //Draw the deviation of the numu-bar flux 
    //    TCanvas *c5 = new TCanvas("c5","c5",600,500);
    //    flux_dev_numub->SetStats(false);
    //    flux_dev_numub->SetLineWidth(2);
    //    flux_dev_numub->SetLineColor(2);
    //    flux_dev_numub->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    //    flux_dev_numub->GetYaxis()->SetTitle("Fractional Deviation");
    //    flux_dev_numub->Draw("hist");

    //Draw the nue flux
    //    TCanvas *c6 = new TCanvas("c6","c6",600,500);
    //    fitted_flux_nue->SetLineWidth(2);
    //    fitted_flux_nue->SetLineColor(2);
    //    sk_flux_nue->SetLineWidth(2);
    //    sk_flux_nue->SetStats(false);
    //    sk_flux_nue->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    //    sk_flux_nue->GetYaxis()->SetTitle("Flux/[cm^{2}#upoint 100 MEV #upoint 1e21 POT]");
    //    sk_flux_nue->Draw();
    //    fitted_flux_nue->Draw("same hist"); 

    //    char* outname = Form("nuprism_coef_%d_%d.root",int(round(dm2*100000)), int(round(theta*100000)));
    char* outname = Form("nuprism_coef_oscSpectrum.root");
    TFile fout(outname,"RECREATE");
    sk_flux->Write();
    fitted_flux->Write();
    coefficients->Write("Coefficients");
    fout.Close();

//    gApplication->Terminate();

    return 0;
}
