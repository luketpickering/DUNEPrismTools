#!/usr/bin/python

import numpy as np
from array import array
import ROOT

g = ROOT.TGraphErrors()


h2d = ROOT.TH2D("h2d", "h2d", 56, -0.25, 7.-0.25, 250, 0., 10.)

for mu in np.arange(0.5, 6.25, 0.25) :
    fIn = ROOT.TFile("nuprism_coef_"+str(mu)+".root")
#    fIn.FittedFlux.Fit("gaus", "QM")
    func = fIn.FittedFlux.Fit("gaus", "SQM", "", 0.1, 9.0)
    
    mean = func.GetParams()[1]
    sigma = func.GetParams()[2]

    print mu, mean, sigma
    n = g.GetN()
    g.SetPoint(n, mu, mean)
    g.SetPointError(n, 0., sigma)

    nBinHist = h2d.GetXaxis().FindBin(mu)
    for i in range(1, 251) :
        h2d.SetBinContent(nBinHist, i, fIn.FittedFlux.GetBinContent(i))
        

cgraph = ROOT.TCanvas("cgraph")
g.SetLineWidth(2)
g.SetLineColor(ROOT.kRed)
g.Draw("ALP")

ROOT.gStyle.SetPalette(ROOT.kSunset)

NRGBs = 2
NCont = 99
stops = [ 0.00,  1.00 ]
red   = [ 1.,     0.0 ]
green = [ 1.,     0.0 ]
blue  = [ 1.,     1.0 ]
stopsArray = array('d', stops)
redArray   = array('d', red)
greenArray = array('d', green)
blueArray  = array('d', blue)
ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)
ROOT.gStyle.SetNumberContours(NCont)


h2d.Draw("SAMECOLZ")
h2d.SetMinimum(0.)
h2d.RebinY(2)

lTargetMean = ROOT.TLine(0.0, 0.0, 6., 6.)
lTargetm1 = ROOT.TLine(0.0, 0.0, 6., 5.4)
lTargetp1 = ROOT.TLine(0.0, 0.0, 6., 6.6)

lTargetMean.SetLineColor(ROOT.kBlack)
lTargetm1.SetLineColor(ROOT.kBlack)
lTargetp1.SetLineColor(ROOT.kBlack)

lTargetm1.SetLineStyle(2)
lTargetp1.SetLineStyle(2)

lTargetMean.SetLineWidth(2)
lTargetm1.SetLineWidth(2)
lTargetp1.SetLineWidth(2)

lTargetMean.Draw()
lTargetm1.Draw()
lTargetp1.Draw()

g.Draw("LP")



raw_input()
