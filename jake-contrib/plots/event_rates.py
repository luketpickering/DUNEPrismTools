from ROOT import *
#from os import listdir
import os 
from array import *
from arg_parser import init_parser
import xml.etree.ElementTree
from doEff import doEff
from stops_parser import parse_stops

args = init_parser().parse_args()
         
det_file = args.config


det_configs = parse_stops(det_file)

print det_configs    
datadir = '/home/calcuttj/DUNEPrismSim/'
setting = args.setting
veto = args.veto


DIR = args.DIR
print "DIR:\n\t", DIR
output = os.listdir(datadir + setting + "/" +DIR)

if not os.path.isdir("OAA_events/"+setting+"/"+DIR):
  os.mkdir("OAA_events/"+setting+"/"+DIR)


dc = det_configs[0]
default_size = str(dc['x']/100) + "x" + str(dc['y']/100) + "x" + str(dc['z']/100) +"m"
default_x = dc['x']/100

stops = dict()
for dc in det_configs:
  if dc['shift'] == 0:
    stops[dc['stop']] = "On Axis"
  else:
    stops[dc['stop']] = str(abs(dc['shift']/100)) + "m Off Axis"
chain = TChain("EDeps")

for f in output:
  if len(f.split('.')) == 7:
    chain.Add(datadir + setting + "/" + DIR + "/" + f)

total = chain.GetEntries()
print total 
print "Loaded ",len(output), "files"

out = args.o
fout = TFile(out,"RECREATE")

gStyle.SetOptStat(0)
gStyle.SetPalette(54)
c1 = TCanvas("c1","c1")


hCCEvents = dict()
hCCTotal = dict()
hNCTotal = dict()
hCCEventsHadAcc = dict()
hNCEvents = dict()
hNCEventsHadAcc = dict()
hKE = dict()

fout.cd()

flavors = [12,-12,14,-14]

flavnames = {12:"Nue",-12:"Nuebar",14:"Numu",-14:"Numubar"}
flavcutsCC =  {14:"LepPDG == 13 && NuPDG == 14",
               -14:"LepPDG == -13 && NuPDG == -14",
               12:"LepPDG == 11 && NuPDG == 12",
               -12:"LepPDG == -11 && NuPDG == -12"}
containcut = {14:" && !(flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow || flagLepExitXHigh || flagLepExitXLow)",
              -14:" && !(flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow || flagLepExitXHigh || flagLepExitXLow)",
              12:"&& eElectronShowerDepOutside/(eElectronShowerDepInside + eElectronShowerDepOutside) < 0.05",
              -12:"&& eElectronShowerDepOutside/(eElectronShowerDepInside + eElectronShowerDepOutside) < 0.05"}

exitcut = " && (flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow || flagLepExitXHigh || flagLepExitXLow) && sqrt(lepExitingMomX*lepExitingMomX + lepExitingMomY*lepExitingMomY + lepExitingMomZ*lepExitingMomZ) > 0.114"
vetocut = " && (TotalNonlep_Dep_veto) <= " +str(veto)
flavcutsNC =  {14:"LepPDG == 14 && NuPDG == 14",
               -14:"LepPDG == -14 && NuPDG == -14",
               12:"LepPDG == 12 && NuPDG == 12",
               -12:"LepPDG == -12 && NuPDG == -12"}
#for stop in stops:
for stop in [0]:
  stopcut = "stop == " + str(stop)
  for flavor in flavors:
    print stopcut, flavor


    chain.Draw("Enu>>hCC" + flavnames[flavor] + "Total" + str(stop) + "(100,0,40)",flavcutsCC[flavor] + " && " + stopcut)
    hCCTotal[stop] = gDirectory.Get("hCC" + flavnames[flavor] + "Total" + str(stop))
    hCCTotal[stop].SetTitle(str(hCCTotal[stop].Integral()))
    hCCTotal[stop].GetXaxis().SetTitle("Enu (GeV)")
    hCCTotal[stop].Write()

    chain.Draw("Enu>>hNC" + flavnames[flavor] + "Total" + str(stop) + "(100,0,40)",flavcutsNC[flavor] + " && " + stopcut)
    hNCTotal[stop] = gDirectory.Get("hNC" + flavnames[flavor] + "Total" + str(stop))
    hNCTotal[stop].SetTitle(str(hNCTotal[stop].Integral()))
    hNCTotal[stop].GetXaxis().SetTitle("Enu (GeV)")
    hNCTotal[stop].Write()

    chain.Draw("Enu>>hCC" + flavnames[flavor] + "EventsMuonContained" + str(stop) + "(100,0,40)",flavcutsCC[flavor] + containcut[flavor] + vetocut + " && " + stopcut)
    hCCEvents[stop] = gDirectory.Get("hCC" + flavnames[flavor] + "EventsMuonContained" + str(stop))
    hCCEvents[stop].SetTitle(str(hCCEvents[stop].Integral()))
    hCCEvents[stop].GetXaxis().SetTitle("Enu (GeV)")
    hCCEvents[stop].Write()

    chain.Draw("Enu>>hCC" + flavnames[flavor] + "EventsMuonExit" + str(stop) + "(100,0,40)",flavcutsCC[flavor] + exitcut + vetocut + " && " + stopcut)
    hCCEventsHadAcc[stop] = gDirectory.Get("hCC" + flavnames[flavor] + "EventsMuonExit" + str(stop))
    hCCEventsHadAcc[stop].SetTitle(str(hCCEventsHadAcc[stop].Integral()))
    hCCEventsHadAcc[stop].GetXaxis().SetTitle("Enu (GeV)")
    hCCEventsHadAcc[stop].Write()

    chain.Draw("Enu>>hNC" + flavnames[flavor] + "Events" + str(stop) + "(100,0,40)",flavcutsNC[flavor] + vetocut + " && " + stopcut)
    hNCEvents[stop] = gDirectory.Get("hNC" + flavnames[flavor] + "Events" + str(stop))
    hNCEvents[stop].SetTitle(str(hNCEvents[stop].Integral()))
    hNCEvents[stop].GetXaxis().SetTitle("Enu (GeV)")
    hNCEvents[stop].Write()

    chain.Draw("sqrt(lepExitingMomX*lepExitingMomX + lepExitingMomY*lepExitingMomY + lepExitingMomZ*lepExitingMomZ + .105*.105) - .105>>hKE"+flavnames[flavor]+str(stop)+"(50,0,5)","flagLepExitBack && " + stopcut)
    hKE[stop] = gDirectory.Get("hKE"+flavnames[flavor]+str(stop))
    hKE[stop].Write()

fout.Close()
