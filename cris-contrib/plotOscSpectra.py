import ROOT
import sys, os


if len(sys.argv) != 2 :
    print "No filename given"
    exit(-1)
else :
    print sys.argv[1]


fIn = ROOT.TFile(sys.argv[1])

fIn.OscillatedFlux.SetLineColor(ROOT.kBlue)
fIn.OscillatedFlux.SetTitle("")

c1 = ROOT.TCanvas()

fIn.OscillatedFlux.Draw("HIST")
fIn.fitted_flux.Draw("SAMEHIST")

outFname = os.path.splitext(sys.argv[1])

c1.SaveAs(outFname[0]+".eps")
c1.SaveAs(outFname[0]+".png")
c1.SaveAs(outFname[0]+".pdf")
