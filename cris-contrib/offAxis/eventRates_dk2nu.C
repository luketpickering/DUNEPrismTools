#define eventRates_dk2nu_cxx
#include "eventRates_dk2nu.h"
#include "OscLib/OscCalculator.cxx"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <TGaxis.h>
#include <TParameter.h>

// EQUIVALENT PDG CODES for DK2NU!
const int    nue     =  12;
const int    nuebar  =  -12;
const int    numu    =  14;
const int    numubar =  -14;
const int    nutau    =  16; // is this right?
const int    nutaubar =  -16; // is this right?
const int    muplus  =  -13;
const int    muminus =  13;
const int    piplus = 211;
const int    piminus = -211;
const int    kplus = 321;
const int    kminus = -321;
const int    kzero = 130;

using namespace std;

int iread = 0;

void eventRates_dk2nu::Loop()
{

  //   In a ROOT session, you can do:
  //      Root > .L eventRates.C++  // The ++ is important because the code runs MUCH faster when compiled
  //      Root > eventRates t
  //      Root > t.Loop();       // Loop on all entries
  //

   //
   // Define the Reference pot 
   //

   double refpot          = 1;
   std::string potref_str = eventRates_dk2nu::GetPOTAsString(refpot);

   //
   // Set Histogram Binning
   //

   // simple binning for some alignment plots
   int nbins   = 1000;
   double xmin = 0.0;
   double xmax = 10.0;

   // more complex binning for passing to fastMC
   std::vector< Double_t > fastmc_bins;
   // 0.125 GeV bins up to 8 GeV
   for(int i = 0; i<8/0.125; i++)
     fastmc_bins.push_back(i*.125);
   // 0.5 GeV bins up to 20 GeV
   for(int i = 0; i<(20-8)/0.5; i++)
     fastmc_bins.push_back(8.0+i*.5);
  // 2.0 GeV bins up to 120 GeV
   for(int i = 0; i<(120-20)/2.0; i++)
     fastmc_bins.push_back(20.0+i*2.0);

   //
   // Declare all the histograms
   //

   TH1D *fhNuMuFlux    = new TH1D("numu_flux_forplots",  
				  "numu_flux_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarFlux = new TH1D("numubar_flux_forplots", 
				  "numubar_flux_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEFlux     = new TH1D("nue_flux_forplots",
				  "nue_flux_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarFlux  = new TH1D("nuebar_flux_forplots",
				  "nuebar_flux_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauFlux     = new TH1D("nutau_flux_forplots",
				    "nutau_flux_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarFlux  = new TH1D("nutaubar_flux_forplots",
				    "nutaubar_flux_forplots", nbins/2,xmin,xmax);
   
   TH1D *fhNuMuCCEventRate    = new TH1D("numu_cceventrate_forplots",  
				  "numu_cceventrate_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarCCEventRate = new TH1D("numubar_cceventrate_forplots", 
				  "numubar_cceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuECCEventRate     = new TH1D("nue_cceventrate_forplots",
				  "nue_cceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarCCEventRate  = new TH1D("nuebar_cceventrate_forplots",
				  "nuebar_cceventrate_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauCCEventRate     = new TH1D("nutau_cceventrate_forplots",
				    "nutau_cceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarCCEventRate  = new TH1D("nutaubar_cceventrate_forplots",
				    "nutaubar_cceventrate_forplots", nbins/2,xmin,xmax);

   TH1D *fhNuMuNCEventRate    = new TH1D("numu_nceventrate_forplots",  
				  "numu_nceventrate_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarNCEventRate = new TH1D("numubar_nceventrate_forplots", 
				  "numubar_nceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuENCEventRate     = new TH1D("nue_nceventrate_forplots",
				  "nue_nceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarNCEventRate  = new TH1D("nuebar_nceventrate_forplots",
				  "nuebar_nceventrate_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauNCEventRate     = new TH1D("nutau_nceventrate_forplots",
				    "nutau_nceventrate_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarNCEventRate  = new TH1D("nutaubar_nceventrate_forplots",
				    "nutaubar_nceventrate_forplots", nbins/2,xmin,xmax);


   TH1D *fhNuMuFlux_FastMC    = new TH1D("numu_flux",  
					 "numu_flux",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarFlux_FastMC = new TH1D("numubar_flux", 
					 "numubar_flux", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEFlux_FastMC     = new TH1D("nue_flux",
					 "nue_flux", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarFlux_FastMC  = new TH1D("nuebar_flux",
					 "nuebar_flux", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauFlux_FastMC     = new TH1D("nutau_flux",
					   "nutau_flux", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarFlux_FastMC  = new TH1D("nutaubar_flux",
					   "nutaubar_flux", fastmc_bins.size()-1,&fastmc_bins[0]);

   TH1D *fhNuMuCCEventRate_FastMC    = new TH1D("numu_cceventrate",  
					 "numu_cceventrate",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarCCEventRate_FastMC = new TH1D("numubar_cceventrate", 
					 "numubar_cceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuECCEventRate_FastMC     = new TH1D("nue_cceventrate",
					 "nue_cceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarCCEventRate_FastMC  = new TH1D("nuebar_cceventrate",
					 "nuebar_cceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauCCEventRate_FastMC     = new TH1D("nutau_cceventrate",
					   "nutau_cceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarCCEventRate_FastMC  = new TH1D("nutaubar_cceventrate",
					   "nutaubar_cceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);

   TH1D *fhNuMuNCEventRate_FastMC    = new TH1D("numu_nceventrate",  
					 "numu_nceventrate",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarNCEventRate_FastMC = new TH1D("numubar_nceventrate", 
					 "numubar_nceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuENCEventRate_FastMC     = new TH1D("nue_nceventrate",
					 "nue_nceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarNCEventRate_FastMC  = new TH1D("nuebar_nceventrate",
					 "nuebar_nceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauNCEventRate_FastMC     = new TH1D("nutau_nceventrate",
					   "nutau_nceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarNCEventRate_FastMC  = new TH1D("nutaubar_nceventrate",
					   "nutaubar_nceventrate", fastmc_bins.size()-1,&fastmc_bins[0]);

   double globes_low = 0;
   double globes_high = 125.25;
   double globes_binwidth = 0.25;
   int globes_nbins = ( globes_high - globes_low ) / globes_binwidth;

   TH1D *fhNuMuFlux_Globes    = new TH1D("numu_flux_globes",  
					 "numu_flux_globes", globes_nbins,globes_low,globes_high);
   TH1D *fhNuMuBarFlux_Globes = new TH1D("numubar_flux_globes", 
				  "numubar_flux_globes", globes_nbins, globes_low, globes_high);
   TH1D *fhNuEFlux_Globes     = new TH1D("nue_flux_globes",
				  "nue_flux_globes", globes_nbins, globes_low, globes_high);
   TH1D *fhNuEBarFlux_Globes  = new TH1D("nuebar_flux_globes",
				  "nuebar_flux_globes", globes_nbins, globes_low, globes_high);
 
   TH1D *fhNuTauFlux_Globes     = new TH1D("nutau_flux_globes",
				    "nutau_flux_globes", globes_nbins, globes_low, globes_high);
   TH1D *fhNuTauBarFlux_Globes  = new TH1D("nutaubar_flux_globes",
					   "nutau_flux_globes", globes_nbins, globes_low, globes_high);

   TH1D *fhNuMuFluxOsc    = new TH1D("numu_fluxosc_forplots",  
				  "numu_fluxosc_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarFluxOsc = new TH1D("numubar_fluxosc_forplots", 
				  "numubar_fluxosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEFluxOsc     = new TH1D("nue_fluxosc_forplots",
				  "nue_fluxosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarFluxOsc  = new TH1D("nuebar_fluxosc_forplots",
				  "nuebar_fluxosc_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauFluxOsc     = new TH1D("nutau_fluxosc_forplots",
				    "nutau_fluxosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarFluxOsc  = new TH1D("nutaubar_fluxosc_forplots",
				    "nutaubar_fluxosc_forplots", nbins/2,xmin,xmax);   

   TH1D *fhNuMuCCEventRateOsc    = new TH1D("numu_cceventrateosc_forplots",  
				  "numu_cceventrateosc_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarCCEventRateOsc = new TH1D("numubar_cceventrateosc_forplots", 
				  "numubar_cceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuECCEventRateOsc     = new TH1D("nue_cceventrateosc_forplots",
				  "nue_cceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarCCEventRateOsc  = new TH1D("nuebar_cceventrateosc_forplots",
				  "nuebar_cceventrateosc_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauCCEventRateOsc     = new TH1D("nutau_cceventrateosc_forplots",
				    "nutau_cceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarCCEventRateOsc  = new TH1D("nutaubar_cceventrateosc_forplots",
				    "nutaubar_cceventrateosc_forplots", nbins/2,xmin,xmax);

   TH1D *fhNuMuNCEventRateOsc    = new TH1D("numu_nceventrateosc_forplots",  
				  "numu_nceventrateosc_forplots", nbins,xmin,xmax);
   TH1D *fhNuMuBarNCEventRateOsc = new TH1D("numubar_nceventrateosc_forplots", 
				  "numubar_nceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuENCEventRateOsc     = new TH1D("nue_nceventrateosc_forplots",
				  "nue_nceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuEBarNCEventRateOsc  = new TH1D("nuebar_nceventrateosc_forplots",
				  "nuebar_nceventrateosc_forplots", nbins/2,xmin,xmax);
 
   TH1D *fhNuTauNCEventRateOsc     = new TH1D("nutau_nceventrateosc_forplots",
				    "nutau_nceventrateosc_forplots", nbins/2,xmin,xmax);
   TH1D *fhNuTauBarNCEventRateOsc  = new TH1D("nutaubar_nceventrateosc_forplots",
				    "nutaubar_nceventrateosc_forplots", nbins/2,xmin,xmax);

   TH1D *fhNuMuFluxOsc_FastMC    = new TH1D("numu_fluxosc",  
					 "numu_fluxosc",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarFluxOsc_FastMC = new TH1D("numubar_fluxosc", 
					 "numubar_fluxosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEFluxOsc_FastMC     = new TH1D("nue_fluxosc",
					 "nue_fluxosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarFluxOsc_FastMC  = new TH1D("nuebar_fluxosc",
					 "nuebar_fluxosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauFluxOsc_FastMC     = new TH1D("nutau_fluxosc",
					   "nutau_fluxosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarFluxOsc_FastMC  = new TH1D("nutaubar_fluxosc",
					   "nutaubar_fluxosc", fastmc_bins.size()-1,&fastmc_bins[0]);

   TH1D *fhNuMuCCEventRateOsc_FastMC    = new TH1D("numu_cceventrateosc",  
					 "numu_cceventrateosc",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarCCEventRateOsc_FastMC = new TH1D("numubar_cceventrateosc", 
					 "numubar_cceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuECCEventRateOsc_FastMC     = new TH1D("nue_cceventrateosc",
					 "nue_cceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarCCEventRateOsc_FastMC  = new TH1D("nuebar_cceventrateosc",
					 "nuebar_cceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauCCEventRateOsc_FastMC     = new TH1D("nutau_cceventrateosc",
					   "nutau_cceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarCCEventRateOsc_FastMC  = new TH1D("nutaubar_cceventrateosc",
					   "nutaubar_cceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);

   TH1D *fhNuMuNCEventRateOsc_FastMC    = new TH1D("numu_nceventrateosc",  
					 "numu_nceventrateosc",fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuMuBarNCEventRateOsc_FastMC = new TH1D("numubar_nceventrateosc", 
					 "numubar_nceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuENCEventRateOsc_FastMC     = new TH1D("nue_nceventrateosc",
					 "nue_nceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuEBarNCEventRateOsc_FastMC  = new TH1D("nuebar_nceventrateosc",
					 "nuebar_nceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauNCEventRateOsc_FastMC     = new TH1D("nutau_nceventrateosc",
					   "nutau_nceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);
   TH1D *fhNuTauBarNCEventRateOsc_FastMC  = new TH1D("nutaubar_nceventrateosc",
					   "nutaubar_nceventrateosc", fastmc_bins.size()-1,&fastmc_bins[0]);

   //
   // Call Sumw2 to make sure histogram errors are propagated correctly
   //
 
   fhNuMuFlux->Sumw2();   
   fhNuMuBarFlux->Sumw2();   
   fhNuEFlux->Sumw2();   
   fhNuEBarFlux->Sumw2();   
   fhNuTauFlux->Sumw2();   
   fhNuTauBarFlux->Sumw2();   

   fhNuMuFlux_FastMC->Sumw2();   
   fhNuMuBarFlux_FastMC->Sumw2();   
   fhNuEFlux_FastMC->Sumw2();   
   fhNuEBarFlux_FastMC->Sumw2();   
   fhNuTauFlux_FastMC->Sumw2();   
   fhNuTauBarFlux_FastMC->Sumw2();   

   fhNuMuFlux_Globes->Sumw2();   
   fhNuMuBarFlux_Globes->Sumw2();   
   fhNuEFlux_Globes->Sumw2();   
   fhNuEBarFlux_Globes->Sumw2();   
   fhNuTauFlux_Globes->Sumw2();   
   fhNuTauBarFlux_Globes->Sumw2();   

   fhNuMuFluxOsc->Sumw2();   
   fhNuMuBarFluxOsc->Sumw2();   
   fhNuEFluxOsc->Sumw2();   
   fhNuEBarFluxOsc->Sumw2();   
   fhNuTauFluxOsc->Sumw2();   
   fhNuTauBarFluxOsc->Sumw2();   

   fhNuMuFluxOsc_FastMC->Sumw2();   
   fhNuMuBarFluxOsc_FastMC->Sumw2();   
   fhNuEFluxOsc_FastMC->Sumw2();   
   fhNuEBarFluxOsc_FastMC->Sumw2();   
   fhNuTauFluxOsc_FastMC->Sumw2();   
   fhNuTauBarFluxOsc_FastMC->Sumw2();   

   fhNuMuCCEventRate->Sumw2();   
   fhNuMuBarCCEventRate->Sumw2();   
   fhNuECCEventRate->Sumw2();   
   fhNuEBarCCEventRate->Sumw2();   
   fhNuTauCCEventRate->Sumw2();   
   fhNuTauBarCCEventRate->Sumw2();   
   fhNuMuNCEventRate->Sumw2();   
   fhNuMuBarNCEventRate->Sumw2();   
   fhNuENCEventRate->Sumw2();   
   fhNuEBarNCEventRate->Sumw2();   
   fhNuTauNCEventRate->Sumw2();   
   fhNuTauBarNCEventRate->Sumw2();   

   fhNuMuCCEventRate_FastMC->Sumw2();   
   fhNuMuBarCCEventRate_FastMC->Sumw2();   
   fhNuECCEventRate_FastMC->Sumw2();   
   fhNuEBarCCEventRate_FastMC->Sumw2();   
   fhNuTauCCEventRate_FastMC->Sumw2();   
   fhNuTauBarCCEventRate_FastMC->Sumw2();   
   fhNuMuNCEventRate_FastMC->Sumw2();   
   fhNuMuBarNCEventRate_FastMC->Sumw2();   
   fhNuENCEventRate_FastMC->Sumw2();   
   fhNuEBarNCEventRate_FastMC->Sumw2();   
   fhNuTauNCEventRate_FastMC->Sumw2();   
   fhNuTauBarNCEventRate_FastMC->Sumw2();   

   fhNuMuCCEventRateOsc->Sumw2();   
   fhNuMuBarCCEventRateOsc->Sumw2();   
   fhNuECCEventRateOsc->Sumw2();   
   fhNuEBarCCEventRateOsc->Sumw2();   
   fhNuTauCCEventRateOsc->Sumw2();   
   fhNuTauBarCCEventRateOsc->Sumw2();   
   fhNuMuNCEventRateOsc->Sumw2();   
   fhNuMuBarNCEventRateOsc->Sumw2();   
   fhNuENCEventRateOsc->Sumw2();   
   fhNuEBarNCEventRateOsc->Sumw2();   
   fhNuTauNCEventRateOsc->Sumw2();   
   fhNuTauBarNCEventRateOsc->Sumw2();   

   fhNuMuCCEventRateOsc_FastMC->Sumw2();   
   fhNuMuBarCCEventRateOsc_FastMC->Sumw2();   
   fhNuECCEventRateOsc_FastMC->Sumw2();   
   fhNuEBarCCEventRateOsc_FastMC->Sumw2();   
   fhNuTauCCEventRateOsc_FastMC->Sumw2();   
   fhNuTauBarCCEventRateOsc_FastMC->Sumw2();   
   fhNuMuNCEventRateOsc_FastMC->Sumw2();   
   fhNuMuBarNCEventRateOsc_FastMC->Sumw2();   
   fhNuENCEventRateOsc_FastMC->Sumw2();   
   fhNuEBarNCEventRateOsc_FastMC->Sumw2();   
   fhNuTauNCEventRateOsc_FastMC->Sumw2();   
   fhNuTauBarNCEventRateOsc_FastMC->Sumw2();   


   //
   // Set histogram titles
   //

   std::string fluxtitle      = "Neutrinos / GeV / m^{2} / POT";
   std::string oscfluxtitle      = "Oscillated Neutrinos / GeV / m^{2} / POT";
   std::string cceventratetitle      = "CC Events / POT";
   std::string nceventratetitle      = "CC Events / POT";

   SetTitles(fhNuMuFlux,         "#nu_{#mu} Energy (GeV) ", "Unosc #nu_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuMuBarFlux,      "#bar{#nu}_{#mu} Energy (GeV)", "Unosc #bar{#nu}_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuEFlux,          "#nu_{e} Energy (GeV)", "Unosc #nu_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuEBarFlux,       "#bar{#nu}_{e} Energy (GeV)", "Unosc #bar{#nu}_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauFlux,          "#nu_{#tau} Energy (GeV)", "Unosc #nu_{#tau}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauBarFlux,       "#bar{#nu}_{#tau} Energy (GeV)", "Unosc #bar{#nu}_{#tau}s / GeV / m^{2} / POT");

   SetTitles(fhNuMuFlux_FastMC,         "#nu_{#mu} Energy (GeV) ", "Unosc #nu_{#mu}s / m^{2} / POT");
   SetTitles(fhNuMuBarFlux_FastMC,      "#bar{#nu}_{#mu} Energy (GeV)", "Unosc #bar{#nu}_{#mu}s / m^{2} / POT");
   SetTitles(fhNuEFlux_FastMC,          "#nu_{e} Energy (GeV)", "Unosc #nu_{e}s / m^{2} / POT");
   SetTitles(fhNuEBarFlux_FastMC,       "#bar{#nu}_{e} Energy (GeV)", "Unosc #bar{#nu}_{e}s / m^{2} / POT");
   SetTitles(fhNuTauFlux_FastMC,          "#nu_{#tau} Energy (GeV)", "Unosc #nu_{#tau}s / m^{2} / POT");
   SetTitles(fhNuTauBarFlux_FastMC,       "#bar{#nu}_{#tau} Energy (GeV)", "Unosc #bar{#nu}_{#tau}s / m^{2} / POT");

   SetTitles(fhNuMuFlux_Globes,         "#nu_{#mu} Energy (GeV)", "Unosc #nu_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuMuBarFlux_Globes,      "#bar{#nu}_{#mu} Energy (GeV)", "Unosc #bar{#nu}_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuEFlux_Globes,          "#nu_{e} Energy (GeV)", "Unosc #nu_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuEBarFlux_Globes,       "#bar{#nu}_{e} Energy (GeV)", "Unosc #bar{#nu}_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauFlux_Globes,          "#nu_{#tau} Energy (GeV)", "Unosc #nu_{#tau}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauBarFlux_Globes,       "#bar{#nu}_{#tau} Energy (GeV)", "Unosc #bar{#nu}_{#tau}s / GeV / m^{2} / POT");

   SetTitles(fhNuMuFluxOsc,         "Energy (GeV)", "Oscillated #nu_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuMuBarFluxOsc,      "Energy (GeV)", "Oscillated #bar{#nu}_{#mu}s / GeV / m^{2} / POT");
   SetTitles(fhNuEFluxOsc,          "Energy (GeV)", "Oscillated #nu_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuEBarFluxOsc,       "Energy (GeV)", "Oscillated #bar{#nu}_{e}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauFluxOsc,          "Energy (GeV)", "Oscillated #nu_{#tau}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauBarFluxOsc,       "Energy (GeV)", "Oscillated #bar{#nu}_{#tau}s / GeV / m^{2} / POT");

   SetTitles(fhNuMuFluxOsc_FastMC,         "Energy (GeV)", "Oscillated #nu_{#mu}s / m^{2} / POT");
   SetTitles(fhNuMuBarFluxOsc_FastMC,      "Energy (GeV)", "Oscillated #bar{#nu}_{#mu}s / m^{2} / POT");
   SetTitles(fhNuEFluxOsc_FastMC,          "Energy (GeV)", "Oscillated #nu_{e}s / m^{2} / POT");
   SetTitles(fhNuEBarFluxOsc_FastMC,       "Energy (GeV)", "Oscillated #bar{#nu}_{e}s / m^{2} / POT");
   SetTitles(fhNuTauFluxOsc_FastMC,          "Energy (GeV)", "Oscillated #nu_{#tau}s / GeV / m^{2} / POT");
   SetTitles(fhNuTauBarFluxOsc_FastMC,       "Energy (GeV)", "Oscillated #bar{#nu}_{#tau}s / GeV / m^{2} / POT");


   SetTitles(fhNuMuCCEventRate,         "Energy (GeV)", "#nu_{#mu} CC Events / GeV / kTon / POT");
   SetTitles(fhNuMuBarCCEventRate,      "Energy (GeV)", "#bar{#nu}_{#mu} CC Events / kTon / POT");
   SetTitles(fhNuECCEventRate,          "Energy (GeV)", "#nu_{e} CC Events / GeV / kTon / POT");
   SetTitles(fhNuEBarCCEventRate,       "Energy (GeV)", "#bar{#nu}_{e} CC Events / kTon / POT");
   SetTitles(fhNuTauCCEventRate,          "Energy (GeV)", "#nu_{#tau} CC Events / kTon / POT");
   SetTitles(fhNuTauBarCCEventRate,       "Energy (GeV)", "#bar{#nu}_{#tau} CC Events / kTon /  POT");

   SetTitles(fhNuMuCCEventRate_FastMC,         "Energy (GeV)", "#nu_{#mu} CC Events / kTon / POT");
   SetTitles(fhNuMuBarCCEventRate_FastMC,      "Energy (GeV)", "#bar{#nu}_{#mu} CC Events / kTon / POT");
   SetTitles(fhNuECCEventRate_FastMC,          "Energy (GeV)", "#nu_{e} CC Events / kTon / POT");
   SetTitles(fhNuEBarCCEventRate_FastMC,       "Energy (GeV)", "#bar{#nu}_{e} CC Events / kTon / POT");
   SetTitles(fhNuTauCCEventRate_FastMC,          "Energy (GeV)", "#nu_{#tau} CC Events / kTon / POT");
   SetTitles(fhNuTauBarCCEventRate_FastMC,       "Energy (GeV)", "#bar{#nu}_{#tau} CC Events / kTon /  POT");

   SetTitles(fhNuMuCCEventRateOsc,         "Energy (GeV)", "Oscillated #nu_{#mu} CC Events / GeV / kTon / POT");
   SetTitles(fhNuMuBarCCEventRateOsc,      "Energy (GeV)", "Oscillated #bar{#nu}_{#mu} CC Events / GeV / kTon / POT");
   SetTitles(fhNuECCEventRateOsc,          "Energy (GeV)", "Oscillated #nu_{e} CC Events / GeV / kTon / POT");
   SetTitles(fhNuEBarCCEventRateOsc,       "Energy (GeV)", "Oscillated #bar{#nu}_{e} CC Events / GeV / kTon / POT");
   SetTitles(fhNuTauCCEventRateOsc,          "Energy (GeV)", "Oscillated #nu_{#tau} CC Events / GeV / kTon / POT");
   SetTitles(fhNuTauBarCCEventRateOsc,       "Energy (GeV)", "Oscillated #bar{#nu}_{#tau} CC Events / GeV / kTon /  POT");

   SetTitles(fhNuMuCCEventRateOsc_FastMC,         "Energy (GeV)", "Oscillated #nu_{#mu} CC Events / kTon / POT");
   SetTitles(fhNuMuBarCCEventRateOsc_FastMC,      "Energy (GeV)", "Oscillated #bar{#nu}_{#mu} CC Events / kTon / POT");
   SetTitles(fhNuECCEventRateOsc_FastMC,          "Energy (GeV)", "Oscillated #nu_{e} CC Events / kTon / POT");
   SetTitles(fhNuEBarCCEventRateOsc_FastMC,       "Energy (GeV)", "Oscillated #bar{#nu}_{e} CC Events / kTon / POT");
   SetTitles(fhNuTauCCEventRateOsc_FastMC,          "Energy (GeV)", "Oscillated #nu_{#tau} CC Events / kTon / POT");
   SetTitles(fhNuTauBarCCEventRateOsc_FastMC,       "Energy (GeV)", "Oscillated #bar{#nu}_{#tau} CC Events / kTon /  POT");


   SetTitles(fhNuMuNCEventRate,         "Energy (GeV)", "#nu_{#mu} NC Events / GeV / kTon / POT");
   SetTitles(fhNuMuBarNCEventRate,      "Energy (GeV)","#bar{#nu}_{#mu} NC Events / GeV / kTon / POT");
   SetTitles(fhNuENCEventRate,          "Energy (GeV)", "#nu_{e} NC Events / GeV / kTon / POT");
   SetTitles(fhNuEBarNCEventRate,       "Energy (GeV)", "#bar{#nu}_{e} NC Events / GeV / kTon / POT");
   SetTitles(fhNuTauNCEventRate,          "Energy (GeV)", "#nu_{#tau} NC Events / GeV / kTon / POT");
   SetTitles(fhNuTauBarNCEventRate,       "Energy (GeV)", "#bar{#nu}_{#tau} NC Events / GeV / kTon / POT");

   SetTitles(fhNuMuNCEventRate_FastMC,         "Energy (GeV)", "#nu_{#mu} NC Events / kTon / POT");
   SetTitles(fhNuMuBarNCEventRate_FastMC,      "#bar{#nu}_{#mu} Energy (GeV)","#bar{#nu}_{#mu} NC Events / kTon / POT");
   SetTitles(fhNuENCEventRate_FastMC,          "#nu_{e} Energy (GeV)", "#nu_{e} NC Events / kTon / POT");
   SetTitles(fhNuEBarNCEventRate_FastMC,       "#bar{#nu}_{e} Energy (GeV)", "#bar{#nu}_{e} NC Events / kTon / POT");
   SetTitles(fhNuTauNCEventRate_FastMC,          "#nu_{#tau} Energy (GeV)", "#nu_{#tau} NC Events / kTon / POT");
   SetTitles(fhNuTauBarNCEventRate_FastMC,       "#bar{#nu}_{#tau} Energy (GeV)", "#bar{#nu}_{#tau} NC Events / kTon / POT");

   SetTitles(fhNuMuNCEventRateOsc,         "#nu_{#mu} Energy (GeV)", "Oscillated #nu_{#mu} NC Events / GeV / kTon / POT");
   SetTitles(fhNuMuBarNCEventRateOsc,      "#bar{#nu}_{#mu} Energy (GeV)","Oscillated #bar{#nu}_{#mu} NC Events / GeV / kTon / POT");
   SetTitles(fhNuENCEventRateOsc,          "#nu_{e} Energy (GeV)", "Oscillated #nu_{e} NC Events / GeV / kTon / POT");
   SetTitles(fhNuEBarNCEventRateOsc,       "#bar{#nu}_{e} Energy (GeV)", "Oscillated #bar{#nu}_{e} NC Events / GeV / kTon / POT");
   SetTitles(fhNuTauNCEventRateOsc,          "#nu_{#tau} Energy (GeV)", "Oscillated #nu_{#tau} NC Events / GeV / kTon / POT");
   SetTitles(fhNuTauBarNCEventRateOsc,       "#bar{#nu}_{#tau} Energy (GeV)", "Oscillated #bar{#nu}_{#tau} NC Events / GeV / kTon / POT");

   SetTitles(fhNuMuNCEventRateOsc_FastMC,         "#nu_{#mu} Energy (GeV)", "Oscillated #nu_{#mu} NC Events / kTon / POT");
   SetTitles(fhNuMuBarNCEventRateOsc_FastMC,      "#bar{#nu}_{#mu} Energy (GeV)","Oscillated #bar{#nu}_{#mu} NC Events / kTon / POT");
   SetTitles(fhNuENCEventRateOsc_FastMC,          "#nu_{e} Energy (GeV)", "Oscillated #nu_{e} NC Events / kTon / POT");
   SetTitles(fhNuEBarNCEventRateOsc_FastMC,       "#bar{#nu}_{e} Energy (GeV)", "Oscillated #bar{#nu}_{e} NC Events / kTon / POT");
   SetTitles(fhNuTauNCEventRateOsc_FastMC,          "#nu_{#tau} Energy (GeV)", "Oscillated #nu_{#tau} NC Events / kTon / POT");
   SetTitles(fhNuTauBarNCEventRateOsc_FastMC,       "#bar{#nu}_{#tau} Energy (GeV)", "Oscillated #bar{#nu}_{#tau} NC Events / kTon / POT");

   //
   //start loop over entries in ntuple
   //

   Long64_t nentries = fChain->GetEntries();
   std::cout << "Total number of Entries = " << nentries << std::endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   //for (Long64_t jentry=0; jentry<1000;jentry++) //fast -- for testing
   {
     //     std::cout << "jentry " << jentry << " Loading tree" << std::endl;
     Long64_t ientry = LoadTree(jentry);
     //     std::cout << "jentry " << jentry << " Tree loaded" << std::endl;
     nb = fChain->GetEntry(jentry);  nbytes += nb;
     //     std::cout << "jentry " << jentry << " Got Entry" << std::endl;
     if (ientry < 0) break;

      ++iread;
      
      if(iread % 10000 == 0)
      {
	 std::cout << "Reading Entry " << iread << std::endl;
      }


      //
      // Compute the detector location weight
      // This calculation needs to know the coordinates
      // of the detector where you want to plot the flux
      // which is set by the eventRates constructor 
      //

      double nuenergyatsomedet     = -999.0;
      double detectorwghtatsomedet = -999.0;
      std::vector<double> detvec;
      detvec.push_back(detx);
      detvec.push_back(dety);
      detvec.push_back(detz);

      //this function computes the location weight and neutrino energy at detx,y,z
      eventRates_dk2nu::GetWeight(detvec, detectorwghtatsomedet, nuenergyatsomedet);
      
      // Calculate the total weight, including location weight, importance weight
      // and POT scale factor

      double fluxwghtsomedet = (detectorwghtatsomedet*decay_nimpwt/3.1415)*(refpot/fTotalPOT);
      
      // Multiply flux weights by cross sections to get event rate std

      string current_string = "CC";
      double cceventratewghtsomedet       = fluxwghtsomedet * GetXSec((double)decay_ntype,nuenergyatsomedet,current_string);
      
      current_string = "NC";
      double nceventratewghtsomedet       = fluxwghtsomedet * GetXSec(decay_ntype,nuenergyatsomedet,current_string);

      // Get Oscillated Flavor 
      int decay_ntype_osc = GetOscillatedNeutrinoType(nuenergyatsomedet);
      ///
      //Now fill the histograms
      //by neutrino type
      //
      if(decay_ntype == numu)
      {
	 fhNuMuFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuCCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuNCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuMuFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuCCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuNCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuMuFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }
      if(decay_ntype == numubar)
      {
	 fhNuMuBarFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuBarCCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuBarNCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuMuBarFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuBarCCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuBarNCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuMuBarFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }
      if(decay_ntype == nue)
      {
	 fhNuEFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuECCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuENCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuEFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuECCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuENCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuEFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }
      if(decay_ntype == nuebar)
      {
	 fhNuEBarFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuEBarCCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuEBarNCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuEBarFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuEBarCCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuEBarNCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuEBarFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }
      
      if(decay_ntype == nutau)
      {
	 fhNuTauFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauCCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauNCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuTauFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauCCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauNCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuTauFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }
      if(decay_ntype == 59)
      {
	 fhNuTauBarFlux -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauBarCCEventRate -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauBarNCEventRate -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

	 fhNuTauBarFlux_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauBarCCEventRate_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauBarNCEventRate_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuTauBarFlux_Globes-> Fill(nuenergyatsomedet, fluxwghtsomedet);
      }

      if(decay_ntype_osc == numu)
      {
	 fhNuMuFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuCCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuNCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuMuFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuCCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuNCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
      }
      if(decay_ntype_osc == numubar)
      {
	 fhNuMuBarFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuBarCCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuBarNCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuMuBarFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuMuBarCCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuMuBarNCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

      }
      if(decay_ntype_osc == nue)
      {
	 fhNuEFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuECCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuENCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuEFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuECCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuENCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
      }
      if(decay_ntype_osc == nuebar)
      {
	 fhNuEBarFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuEBarCCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuEBarNCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuEBarFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuEBarCCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuEBarNCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

      }
      
      if(decay_ntype_osc == nutau)
      {
	 fhNuTauFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauCCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauNCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuTauFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauCCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauNCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);

      }
      if(decay_ntype_osc == 59)
      {
	 fhNuTauBarFluxOsc -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauBarCCEventRateOsc -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauBarNCEventRateOsc -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
	 fhNuTauBarFluxOsc_FastMC -> Fill(nuenergyatsomedet, fluxwghtsomedet);
	 fhNuTauBarCCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, cceventratewghtsomedet);
	 fhNuTauBarNCEventRateOsc_FastMC -> Fill(nuenergyatsomedet, nceventratewghtsomedet);
      }
   }//end loop over entries

   // normalize by bin width for plots but not for fastmc
   fhNuMuFlux->Scale(1.0,"width");
   fhNuMuCCEventRate->Scale(1.0,"width");
   fhNuMuNCEventRate->Scale(1.0,"width");

   fhNuMuBarFlux->Scale(1.0,"width");
   fhNuMuBarCCEventRate->Scale(1.0,"width");
   fhNuMuBarNCEventRate->Scale(1.0,"width");

   fhNuEFlux->Scale(1.0,"width");
   fhNuECCEventRate->Scale(1.0,"width");
   fhNuENCEventRate->Scale(1.0,"width");

   fhNuEBarFlux->Scale(1.0,"width");
   fhNuEBarCCEventRate->Scale(1.0,"width");
   fhNuEBarNCEventRate->Scale(1.0,"width");

   fhNuTauFlux->Scale(1.0,"width");
   fhNuTauCCEventRate->Scale(1.0,"width");
   fhNuTauNCEventRate->Scale(1.0,"width");

   fhNuTauBarFlux->Scale(1.0,"width");
   fhNuTauBarCCEventRate->Scale(1.0,"width");
   fhNuTauBarNCEventRate->Scale(1.0,"width");

   fhNuMuFluxOsc->Scale(1.0,"width");
   fhNuMuCCEventRateOsc->Scale(1.0,"width");
   fhNuMuNCEventRateOsc->Scale(1.0,"width");

   fhNuMuBarFluxOsc->Scale(1.0,"width");
   fhNuMuBarCCEventRateOsc->Scale(1.0,"width");
   fhNuMuBarNCEventRateOsc->Scale(1.0,"width");

   fhNuEFluxOsc->Scale(1.0,"width");
   fhNuECCEventRateOsc->Scale(1.0,"width");
   fhNuENCEventRateOsc->Scale(1.0,"width");

   fhNuEBarFluxOsc->Scale(1.0,"width");
   fhNuEBarCCEventRateOsc->Scale(1.0,"width");
   fhNuEBarNCEventRateOsc->Scale(1.0,"width");

   fhNuTauFluxOsc->Scale(1.0,"width");
   fhNuTauCCEventRateOsc->Scale(1.0,"width");
   fhNuTauNCEventRateOsc->Scale(1.0,"width");

   fhNuTauBarFluxOsc->Scale(1.0,"width");
   fhNuTauBarCCEventRateOsc->Scale(1.0,"width");
   fhNuTauBarNCEventRateOsc->Scale(1.0,"width");

   fhNuMuFlux_Globes->Scale(1.0,"width");
   fhNuMuBarFlux_Globes->Scale(1.0,"width");
   fhNuEFlux_Globes->Scale(1.0,"width");
   fhNuEBarFlux_Globes->Scale(1.0,"width");
   fhNuTauFlux_Globes->Scale(1.0,"width");
   fhNuTauBarFlux_Globes->Scale(1.0,"width");

   //
   // Style histograms
   //

   TGaxis::SetMaxDigits(2);
   fhNuMuFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuMuBarFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuEFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuEBarFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauBarFlux_FastMC->GetYaxis()->SetTitleOffset(1.4);

   fhNuMuCCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuMuBarCCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuECCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuEBarCCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauCCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauBarCCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
 
   fhNuMuNCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuMuBarNCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuENCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuEBarNCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauNCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);
   fhNuTauBarNCEventRate_FastMC->GetYaxis()->SetTitleOffset(1.4);


   // put the location of the detector in the file
   TParameter<double> det_x("det_x",detx);
   TParameter<double> det_y("det_y",dety);
   TParameter<double> det_z("det_z",detz);

   // Save histograms to a root file
   std::cout<<"writing "+ffilename+".root"<<std::endl;
   TFile f((ffilename+".root").c_str(),"recreate");
   det_x.Write();
   det_y.Write();
   det_z.Write();
   fhNuMuFlux->Write();
   fhNuMuBarFlux->Write();
   fhNuEFlux->Write();
   fhNuEBarFlux->Write();
   fhNuTauFlux->Write();
   fhNuTauBarFlux->Write();
   fhNuMuCCEventRate->Write();
   fhNuMuBarCCEventRate->Write();
   fhNuECCEventRate->Write();
   fhNuEBarCCEventRate->Write();
   fhNuTauCCEventRate->Write();
   fhNuTauBarCCEventRate->Write();
   fhNuMuNCEventRate->Write();
   fhNuMuBarNCEventRate->Write();
   fhNuENCEventRate->Write();
   fhNuEBarNCEventRate->Write();
   fhNuTauNCEventRate->Write();
   fhNuTauBarNCEventRate->Write();
   fhNuMuFluxOsc->Write();
   fhNuMuBarFluxOsc->Write();
   fhNuEFluxOsc->Write();
   fhNuEBarFluxOsc->Write();
   fhNuTauFluxOsc->Write();
   fhNuTauBarFluxOsc->Write();
   fhNuMuCCEventRateOsc->Write();
   fhNuMuBarCCEventRateOsc->Write();
   fhNuECCEventRateOsc->Write();
   fhNuEBarCCEventRateOsc->Write();
   fhNuTauCCEventRateOsc->Write();
   fhNuTauBarCCEventRateOsc->Write();
   fhNuMuNCEventRateOsc->Write();
   fhNuMuBarNCEventRateOsc->Write();
   fhNuENCEventRateOsc->Write();
   fhNuEBarNCEventRateOsc->Write();
   fhNuTauNCEventRateOsc->Write();
   fhNuTauBarNCEventRateOsc->Write();
   f.Close();

   TFile g((ffilename+"_fastmc.root").c_str(),"recreate");
   det_x.Write();
   det_y.Write();
   det_z.Write();
   fhNuMuFlux_FastMC->Write();
   fhNuMuBarFlux_FastMC->Write();
   fhNuEFlux_FastMC->Write();
   fhNuEBarFlux_FastMC->Write();
   fhNuTauFlux_FastMC->Write();
   fhNuTauBarFlux_FastMC->Write();
   fhNuMuCCEventRate_FastMC->Write();
   fhNuMuBarCCEventRate_FastMC->Write();
   fhNuECCEventRate_FastMC->Write();
   fhNuEBarCCEventRate_FastMC->Write();
   fhNuTauCCEventRate_FastMC->Write();
   fhNuTauBarCCEventRate_FastMC->Write();
   fhNuMuNCEventRate_FastMC->Write();
   fhNuMuBarNCEventRate_FastMC->Write();
   fhNuENCEventRate_FastMC->Write();
   fhNuEBarNCEventRate_FastMC->Write();
   fhNuTauNCEventRate_FastMC->Write();
   fhNuTauBarNCEventRate_FastMC->Write();

   fhNuMuFluxOsc_FastMC->Write();
   fhNuMuBarFluxOsc_FastMC->Write();
   fhNuEFluxOsc_FastMC->Write();
   fhNuEBarFluxOsc_FastMC->Write();
   fhNuTauFluxOsc_FastMC->Write();
   fhNuTauBarFluxOsc_FastMC->Write();
   fhNuMuCCEventRateOsc_FastMC->Write();
   fhNuMuBarCCEventRateOsc_FastMC->Write();
   fhNuECCEventRateOsc_FastMC->Write();
   fhNuEBarCCEventRateOsc_FastMC->Write();
   fhNuTauCCEventRateOsc_FastMC->Write();
   fhNuTauBarCCEventRateOsc_FastMC->Write();
   fhNuMuNCEventRateOsc_FastMC->Write();
   fhNuMuBarNCEventRateOsc_FastMC->Write();
   fhNuENCEventRateOsc_FastMC->Write();
   fhNuEBarNCEventRateOsc_FastMC->Write();
   fhNuTauNCEventRateOsc_FastMC->Write();
   fhNuTauBarNCEventRateOsc_FastMC->Write();
   g.Close();

   // Write Globes files
   ofstream myfile;
   myfile.open((ffilename+"_globes_flux.txt").c_str());
   for(int i = 0; i<globes_nbins; i++) {
     myfile<<fhNuMuFlux_Globes->GetBinCenter(i+1)<<" "<<
       fhNuEFlux_Globes->GetBinContent(i+1)<<" "<<
       fhNuMuFlux_Globes->GetBinContent(i+1)<<" "<<
       fhNuTauFlux_Globes->GetBinContent(i+1)<<" "<<
       fhNuEBarFlux_Globes->GetBinContent(i+1)<<" "<<
       fhNuMuBarFlux_Globes->GetBinContent(i+1)<<" "<<
       fhNuTauBarFlux_Globes->GetBinContent(i+1)<<" "<<std::endl;
   }
   
}

//-------------------------------------------------------------------------------------
std::string eventRates_dk2nu::GetPOTAsString(const double dpot)
{

   std::stringstream potstrm;
   potstrm << scientific << dpot;

   string potstr = potstrm.str();

   //
   //get base
   //
   size_t baselength;
   if(potstr.find("e",0) != string::npos)
   {
      baselength = potstr.find("e",0);
   }
   else if(potstr.find("E",0) != string::npos)
   {
      baselength = potstr.find("E",0);
   }
   else
   {
      cout << "eventRates::GetPOTAsString - PROBLEM: pot is not in scientific notation" << endl;
      return "Problem";
   }
   
   string base = potstr.substr(0, baselength);

   //
   //get exp
   //
   size_t exppos;
   if(potstr.find("+",baselength) != string::npos)
   {
      exppos = potstr.find("+",baselength);
   }
   else if(potstr.find("-",baselength) != string::npos)
   {
      exppos = potstr.find("-",baselength);
   }
   else
   {
      cout << "eventRates::GetPOTAsString - PROBLEM: pot is not in scientific notation" << endl;
      return "Problem";
   }
   
   string exp = potstr.substr(exppos);


   //
   //modify base string if needed
   //
   string baseNumber = base;   

   size_t baseDecimalpos = base.find(".",0);
   if(baseDecimalpos != string::npos)
   {
      size_t baseNotZeropos = base.find_last_not_of("0", string::npos);
      if(baseNotZeropos != string::npos)
      {
         if(baseNotZeropos > baseDecimalpos)
         {
            baseNumber = base.substr(0,baseNotZeropos+1);
         }
         else
         {
            baseNumber = base.substr(0,baseDecimalpos+2);
         }

      }      
   }
   else
   {
      baseNumber = baseNumber + ".0";
   }

   //
   //modify exp string if needed
   //

   string expSign   = exp.substr(0, 1);
   string expNumber = exp.substr(1, string::npos);

   size_t expNotZeropos = expNumber.find_first_not_of("0",0);
   if(expNotZeropos != string::npos)
   {
      expNumber = expNumber.substr(expNotZeropos, string::npos);
   }



   string potfinalstr;

   //
   //put base and exp together
   //

   if(baseNumber.empty() && expNumber.empty())
   {
      cout << "eventRates::GetPOTAsString - PROBLEM: base number and exp number are both empty" << endl;
      return "Problem";
   }
   
   if(baseNumber == "1.0")
   {
      if(expSign == "-")
         potfinalstr = "10^{" + expSign + expNumber + "}";
      else
         potfinalstr = "10^{" + expNumber + "}";
   }
   else
   {
      if(expSign == "-")
         potfinalstr = baseNumber + "#times10^{" + expSign + expNumber + "}";
      else
         potfinalstr = baseNumber + "#times10^{" + expNumber + "}";
   }

   
   /*
   cout << "pot str = " << potstr << endl
        << " base = " << base << endl
        << " baseNumber = " << baseNumber << endl
        << " expSign = " << expSign << endl
        << " expNumber = " << expNumber << endl
        << " potfinalstr = " << potfinalstr << endl;
   */
   
   
  return potfinalstr;
}


//---------------------------------------------------------------------------------------------
void eventRates_dk2nu::SetTitles(TH1* h, const std::string &xtitle, const std::string &ytitle)
{
   if(!ytitle.empty())
   {
      h -> GetYaxis() -> SetTitle(ytitle.c_str());
      h -> GetYaxis() -> CenterTitle();
   }
   if(!xtitle.empty())
   {
      h -> GetXaxis() -> SetTitle(xtitle.c_str());
      h -> GetXaxis() -> CenterTitle();
   }
}

//---------------------------------------------------------------------------------------------
double eventRates_dk2nu::GetWeight(const std::vector<double> xdet,
			 double& nu_wght, 
			 double& nu_energy)
{

//   if(iread > 60000)
//      std::cout << "start iread = " << iread;

   //assumes units are GeV and cm

   const double rdet    = 100.0; //in cm
   const double pimass  =    0.13957; //in GeV
   const double kmass   =    0.49368;
   const double k0mass  =    0.49767;
   const double mumass  =    0.10389;
   const double taumass =    1.77682;

   //these are geant codes not PDG
   //   const int    nue     =  53;
   //   const int    nuebar  =  nuebar;
   //   const int    numu    =  numu;
   //   const int    numubar =  numubar;
   //   const int    nutau    =  nutau; // is this right?
   //   const int    nutaubar =  59; // is this right?
   //   const int    muplus  =  5;
   //   const int    muminus =  6;
   // DECLARE THESE GLOBALLY!!!!!!!
   // EQUIVALENT PDG CODES for DK2NU!
//   const int    nue     =  12;
//   const int    nuebar  =  -12;
//   const int    numu    =  14;
//   const int    numubar =  -14;
//   const int    nutau    =  16; // is this right?
//   const int    nutaubar =  -16; // is this right?
//   const int    muplus  =  -13;
//   const int    muminus =  13;


   double parent_mass=0.;
   if      (decay_ptype == piplus  || decay_ptype == piminus)  parent_mass = pimass;
   else if (decay_ptype == kplus || decay_ptype == kminus) parent_mass = kmass;
   else if (decay_ptype == kzero)                parent_mass = k0mass;
   else if (decay_ptype == muplus  || decay_ptype == muminus)  parent_mass = mumass;
   else 
   {
      cout <<"eventRates::GetWeight - Wrong parent type!! "<< decay_ptype << " = "
	   << decay_ptype << " Decay code = " << decay_ndecay <<endl;
      
     return -999;
   }
   
   double parent_energy = sqrt(decay_pdpx*decay_pdpx +
			       decay_pdpy*decay_pdpy +
			       decay_pdpz*decay_pdpz + 
			       parent_mass*parent_mass);
   double gamma = parent_energy / parent_mass;
   double gamma_sqr = gamma*gamma;
   double beta_mag = sqrt((gamma_sqr-1.)/gamma_sqr);
   
   double enuzr = decay_necm;
   
   double rad = sqrt((xdet[0]-decay_vx)*(xdet[0]-decay_vx) +
		     (xdet[1]-decay_vy)*(xdet[1]-decay_vy) +
		     (xdet[2]-decay_vz)*(xdet[2]-decay_vz));
   
   double parentp = sqrt((decay_pdpx*decay_pdpx)+
			 (decay_pdpy*decay_pdpy)+
			 (decay_pdpz*decay_pdpz));
   double costh_pardet = (decay_pdpx*(xdet[0]-decay_vx) +
			 decay_pdpy*(xdet[1]-decay_vy) +
			  decay_pdpz*(xdet[2]-decay_vz))/(parentp*rad);

  if (costh_pardet>1.) costh_pardet = 1.;
  else if (costh_pardet<-1.) costh_pardet = -1.;
  double theta_pardet = acos(costh_pardet);

  double emrat = 1./(gamma*(1. - beta_mag * cos(theta_pardet)));

  nu_energy = emrat*enuzr;

  double sangdet = (rdet*rdet /(rad*rad)/ 4.); 

  nu_wght = sangdet * emrat * emrat;

  //done for all except polarized muon
  // in which case need to modify weight
  if (decay_ptype==muplus || decay_ptype==muminus)
  {
     //boost new neutrino to mu decay cm
     double beta[3];
     double p_nu[3]; //nu momentum
     beta[0]=decay_pdpx / parent_energy;
     beta[1]=decay_pdpy / parent_energy;
     beta[2]=decay_pdpz / parent_energy;
     
     p_nu[0] = (xdet[0]- decay_vx) * nu_energy / rad;
     p_nu[1] = (xdet[1]- decay_vy) * nu_energy / rad;
     p_nu[2] = (xdet[2]- decay_vz) * nu_energy / rad;

     double partial = gamma*(beta[0]*p_nu[0]+
			     beta[1]*p_nu[1]+
			     beta[2]*p_nu[2]);
     partial = nu_energy-partial / (gamma+1.);
     double p_dcm_nu[4];
     for (int i=0;i<3;i++) p_dcm_nu[i]=p_nu[i]-beta[i]*gamma*partial;
     p_dcm_nu[3]=0.;
     for (int i=0;i<3;i++) p_dcm_nu[3]+=p_dcm_nu[i]*p_dcm_nu[i];
     p_dcm_nu[3]=sqrt(p_dcm_nu[3]);
     
     //boost parent of mu to mu production cm
     gamma=decay_ppenergy / parent_mass;
     beta[0] = decay_ppdxdz * decay_pppz / decay_ppenergy;
     beta[1] = decay_ppdydz * decay_pppz / decay_ppenergy;
     beta[2] =                  decay_pppz / decay_ppenergy;
     partial = gamma*(beta[0]*decay_muparpx+
		      beta[1]*decay_muparpy+
		      beta[2]*decay_muparpz);
     partial = decay_mupare - partial / (gamma+1.);
     double p_pcm_mp[4];
     p_pcm_mp[0]=decay_muparpx-beta[0]*gamma*partial;
     p_pcm_mp[1]=decay_muparpy-beta[1]*gamma*partial;
     p_pcm_mp[2]=decay_muparpz-beta[2]*gamma*partial;
     p_pcm_mp[3]=0.;
     for (int i=0;i<3;i++) p_pcm_mp[3]+=p_pcm_mp[i]*p_pcm_mp[i];
     p_pcm_mp[3]=sqrt(p_pcm_mp[3]);
     
     double wt_ratio = 1.;
     //have to check p_pcm_mp
     //it can be 0 if mupar..=0. (I guess muons created in target??)
     if (p_pcm_mp[3] != 0. ) {
	//calc new decay angle w.r.t. (anti)spin direction
	double costh = (p_dcm_nu[0]*p_pcm_mp[0]+ 
			p_dcm_nu[1]*p_pcm_mp[1]+ 
			p_dcm_nu[2]*p_pcm_mp[2])/(p_dcm_nu[3]*p_pcm_mp[3]);
	
	if (costh>1.) costh = 1.;
	else if (costh<-1.) costh = -1.;
	
	//calc relative weight due to angle difference
	if (decay_ntype == nue || decay_ntype == nuebar)
	{
	   wt_ratio = 1.-costh;
	}
	else if (decay_ntype == numu || decay_ntype == numubar) 
	{
	   double xnu = 2.* enuzr / mumass;
	   wt_ratio = ( (3.-2.*xnu) - (1.-2.*xnu)*costh ) / (3.-2.*xnu);
	}
	else if (decay_ntype == nutau || decay_ntype == nutaubar) 
	{
	   double xnu = 2.* enuzr / taumass;
	   wt_ratio = ( (3.-2.*xnu) - (1.-2.*xnu)*costh ) / (3.-2.*xnu);
	   std::cout << "calculating weight for tau neutrino; this may not be correct" << std::endl;
	}
	else 
	{
	   std::cout << "eventRates:: Bad neutrino type = " << decay_ntype << std::endl;
	}
     }
     nu_wght *= wt_ratio;
  }
  
//   if(iread > 60000)
//      std::cout << " end iread = " << iread << std::endl;
  
  return nu_wght;
}

double eventRates_dk2nu::GetXSec( double nu_type, 
			    double nu_energy, 
			    std::string current)

{

  if( current != "NC" && current != "CC") {
    cout <<" eventRates::GetXSec: Current other than NC or CC specified... I don't know what to do." << endl;
    return -999;
  }

  int file_index=0;
  if (nu_type == nuebar) file_index = 0;
  if (nu_type == nue) file_index = 1;
  if (nu_type == numubar) file_index = 2;
  if (nu_type == numu) file_index = 3;

  int current_index = 0;
  if( current == "NC")
    current_index = 1;

  // calculate cross section
  double thexsec = 0.;

  // scale factor
  double scale_factor = 6.026e-10;
  //      xseccc = xseccc * 4.09e9
  //     10**-38 cm2 * 10**6 kg/kton * 3.8e20 POT/year *1iron/numu/1.66e-27kg
  //     * 10**-4 m2/cm2 = 4.09e9 conversion factor
  //     since table was 26*sigma(nu-proton) + 30*sigma(nu-neutron)
  //     
  //     10**-38 cm2 * 10**6 kg/kton *1iron/numu/1.66e-27kg
  //     * 10**-4 m2/cm2 = 1.076e-11 conversion factor
  //	
  //     10**-38 cm2 * 10**6 kg/kton * 1 nucleon / 1.66e-27 kg
  //     * 10**-4 m2/cm2 = 6.026e-10

  // if energy is higher than any available xsection point
  // return xsection of highest available energy point
  if ( nu_energy > f_e_arr[fnlines-1][file_index][current_index] ) {
    thexsec = f_xsec_arr[fnlines-1][file_index][current_index]*f_e_arr[fnlines-1][file_index][current_index]*scale_factor;
  }
  // if energy is lower than any available xsection point
  // return xsection of lowest available energy poitn  
  else if ( nu_energy < f_e_arr[0][file_index][current_index] ) {
    thexsec = f_xsec_arr[0][file_index][current_index]*f_e_arr[0][file_index][current_index]*scale_factor;
  }
  else {
    // if not, find the xsections for energy values immendiately above and 
    // below the requested energy
    int energy_index = 0;
    for(int i = 0; i< fnlines-1; i++)
      if( nu_energy > f_e_arr[i][file_index][current_index] && 
	  nu_energy < f_e_arr[i+1][file_index][current_index]) {
	energy_index = i;
	break;
      }

    double sig1 = f_xsec_arr[energy_index][file_index][current_index];
    double sig2 = f_xsec_arr[energy_index+1][file_index][current_index];
	
    thexsec = sig1 + ((sig2 - sig1)/(f_e_arr[energy_index+1][file_index][current_index]-f_e_arr[energy_index][file_index][current_index]))*(nu_energy - f_e_arr[energy_index][file_index][current_index]);
    thexsec = thexsec * nu_energy;

    thexsec = thexsec * scale_factor;
  }

  //std::cout<<"nu_type "<<nu_type<<" current "<<current<<" energy "<<nu_energy <<" xsec "<<thexsec/nu_energy/scale_factor<<std::endl;

  return thexsec;



}


void eventRates_dk2nu::ReadXSecsFromFiles() {

  fnbins = 1500;
  
  std::string base("data/argon_genie2.8.4/");
   
  for(int current = 1; current <= 2; current++) {

  string charge = "cc";

  if ( current == 2 ) {
      charge = "nc";
    }

    const int narr=4;
    string suffix[narr];
    if (current == 1) {
      suffix[0] = "_nuebar.dat";
      suffix[1] = "_nue.dat";
      suffix[2] = "_numubar.dat";
      suffix[3] = "_numu.dat";
    } else if (current==2) {
      suffix[0] = "_numubar.dat";
      suffix[1] = "_numu.dat";
      suffix[2] = "_numubar.dat";
      suffix[3] = "_numu.dat";
    }
    string filename[narr];
    for (int i=0; i<narr; i++) {
      filename[i] = "xsec_"+charge+suffix[i];
      string tmpfilename = base + filename[i];
      
      fdat_file[i].open(tmpfilename.c_str());
      if (fdat_file[i].fail()) {
	cout << " File not found: " << filename[i] << endl;
	//assert(0);
      }
      else {
	cout << " Opened "<<filename[i] << endl;
      }
      double row[2];
      fnlines = 0;
      
      while ( fdat_file[i] >> row[0] >> row[1] ) {
	
	if (fnlines >= fnbins) {
	  cout << " length of data file exceed array size. Fix me. " << filename << endl;
	  assert(0);
	}

	f_e_arr[fnlines][i][current-1] = row[0];
	f_xsec_arr[fnlines][i][current-1] = row[1];
	fnlines++;
      }
      fdat_file[i].close();
    }
  }
}

int eventRates_dk2nu::GetOscillatedNeutrinoType(double E) {

  
  int flavbefore = 0;
//  const int    nue     =  53;
//  const int    nuebar  =  52;
//  const int    numu    =  56;
//  const int    numubar =  55;
//  const int    nutau    =  58; 
//  const int    nutaubar =  59; 
  if (decay_ntype == nuebar) flavbefore = -12;
  else if (decay_ntype == numubar) flavbefore = -14;
  else if (decay_ntype == nutaubar) flavbefore = -16;
  else if (decay_ntype == nue) flavbefore = 12;
  else if (decay_ntype == numu) flavbefore = 14;
  else if (decay_ntype == nutau) flavbefore = 16;
  else {
    std::cout<<"WARNING: unrecognized neutrino type: "<<decay_ntype<<std::endl;
  }

  osc::OscCalculator calculator;
  calculator.SetRho(2.8); // from LBNE Snowmass Document
  calculator.SetL(sqrt(detx*detx+dety*dety+detz*detz)/100/1000); // cm to km
  //  calculator.SetDmsq21(7.50e-5); // PDG 2013
  //  calculator.SetDmsq32(0.00232); // PDG 2013
  calculator.SetDmsq32(0.0025); calculator.SetTh23(0.785); // Set A
  //calculator.SetDmsq32(0.0028); calculator.SetTh23(0.938); // Set B
  //calculator.SetDmsq32(0.0022); calculator.SetTh23(0.938); // Set C
  //calculator.SetDmsq32(0.0025); calculator.SetTh23(0.938); // Set D
  //calculator.SetDmsq32(0.0028); calculator.SetTh23(0.785); // Set E
  //  calculator.SetDmsq32(0.0022); calculator.SetTh23(0.785); // Set F
  
  //  calculator.SetDmsq32(0.0029); // Max
  //  calculator.SetDmsq32(0.0026); // Med
  calculator.SetTh12(0.591); // PDG 2013
  calculator.SetTh13(.16); // PDG 2013
  //  calculator.SetTh23(0.785); // Maximal
  calculator.SetdCP(0);

  double P_e = 0;
  double P_t = 0;
  double P_m = 0;
  
  if(flavbefore<0) {
    P_e = calculator.P(flavbefore, -12, E);
    P_m = calculator.P(flavbefore, -14, E);
    P_t = calculator.P(flavbefore, -16, E);
  }
  else {
    P_e = calculator.P(flavbefore, 12, E);
    P_m = calculator.P(flavbefore, 14, E);
    P_t = calculator.P(flavbefore, 16, E);
  }
  
  // throw a random number to determine which neutrino is produced     
  double random_number = rand3->Rndm();
  
  int ntype_oscillated;
  if(flavbefore <0) {
    if(random_number < P_e)
      ntype_oscillated = nuebar;
    else if(random_number < P_e + P_t)
      ntype_oscillated = nutaubar;
    else
      ntype_oscillated = numubar;
  }
  else {
    if(random_number < P_e)
      ntype_oscillated = nue;
    else if(random_number < P_e + P_t)
      ntype_oscillated = nutau;
    else
      ntype_oscillated = numu;
  }

  return ntype_oscillated;

}

