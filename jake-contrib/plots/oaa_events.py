from ROOT import *
#from os import listdir
import os 
from array import *
from arg_parser import init_parser
import xml.etree.ElementTree
from doEff import doEff
from stops_parser import parse_stops
from load_chain import load_chain
import numpy as np

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


#chain = TChain("EDeps")
#for f in output:
#  if len(f.split('.')) == 7:
#    chain.Add(datadir + setting + "/" + DIR + "/" + f)
chain = load_chain(output, "EDeps",datadir,setting,DIR)
total = chain.GetEntries()
print total 
print "Loaded ",len(output), "files"

gStyle.SetOptStat(0)
gStyle.SetPalette(54)
c1 = TCanvas("c1","c1")

bins = dict()
full_bins = []
#lowerbins = []
#higherbins = []
for dc in det_configs:
  #bins[dc['stop']] = [-1*dc['shift']/100. - dc['x']/200., -1*dc['shift']/100. + dc['x']/200. ]
  bins[dc['stop']] = [-1*dc['shift']/100. - dc['x']/200.,
                      -1*dc['shift']/100. - dc['x']/200. + .1,
                      -1*dc['shift']/100. - dc['x']/200. + .3,
                      -1*dc['shift']/100. - dc['x']/200. + .7,
                      -1*dc['shift']/100. - dc['x']/200. + 1.1,
                      -1*dc['shift']/100. - dc['x']/200. + 1.5,
                      -1*dc['shift']/100. - dc['x']/200. + 1.9,
                      -1*dc['shift']/100. - dc['x']/200. + 2.3,
                      -1*dc['shift']/100. - dc['x']/200. + 2.7,
                      -1*dc['shift']/100. - dc['x']/200. + 2.9,
                      -1*dc['shift']/100. + dc['x']/200.]
  for b in bins[dc['stop']]:
    if b not in full_bins:
      full_bins.append(b)
print "full: ", full_bins                                        
                        
  #lowerbins.append(-1*dc['shift']/100. - dc['x']/200.)
  #higherbins.append(-1*dc['shift']/100. + dc['x']/200.)


hOAAEvents = dict()
hOAAEventsAcc = dict()
hOAAEventsEff = dict()
hOAAEventsMuHadAcc = dict()
hOAAEventsMuHadEff = dict()
hOAAEventsExitAcc = dict()
hOAAEventsExitEff = dict()
hOAAEventsExitMuHadAcc = dict()
hOAAEventsExitMuHadEff = dict()

#start_bin = (min(lowerbins))
#end_bin = (max(higherbins))

#n_bins = int((end_bin - start_bin)*10)
n_bins = len(full_bins) - 1
#print n_bins,start_bin,end_bin,default_x*10

#hOAAEventsFull = TH1D("hOAAEventsFull","",n_bins,start_bin,end_bin)
#hOAAEventsEffFull = TH1D("hOAAEventsEffFull","",n_bins,start_bin,end_bin)
#hOAAEventsMuHadEffFull = TH1D("hOAAEventsMuHadEffFull","",n_bins,start_bin,end_bin)
#hOAAEventsExitEffFull = TH1D("hOAAEventsExitEffFull","",n_bins,start_bin,end_bin)
#hOAAEventsExitMuHadEffFull = TH1D("hOAAEventsExitMuHadEffFull","",n_bins,start_bin,end_bin)
hOAAEventsFull = TH1D("hOAAEventsFull","",n_bins,array('d',full_bins))
hOAAEventsEffFull = TH1D("hOAAEventsEffFull","",n_bins,array('d',full_bins))
hOAAEventsMuHadEffFull = TH1D("hOAAEventsMuHadEffFull","",n_bins,array('d',full_bins))
hOAAEventsExitEffFull = TH1D("hOAAEventsExitEffFull","",n_bins,array('d',full_bins))
hOAAEventsExitMuHadEffFull = TH1D("hOAAEventsExitMuHadEffFull","",n_bins,array('d',full_bins))

exitcut = "(flagLepExitBack || flagLepExitXHigh || flagLepExitXLow || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow) && sqrt(lepExitingMomX*lepExitingMomX + lepExitingMomY*lepExitingMomY + lepExitingMomZ*lepExitingMomZ) > 0.114"

flavcut = {"FHC":"LepPDG == 13 && NuPDG == 14",
           "RHC":"LepPDG == -13 && NuPDG == -14"}

for stop in stops:
#for stop in [0]:
  
  stopcut = "stop == " + str(stop)
  print stopcut
  print array('d',bins[stop])
  hOAAEvents[stop] = TH1D("hOAAEvents" + str(stop), "", len(bins[stop]) - 1,array('d',bins[stop]))
  chain.Draw("vtx[0]/-1.e2>>hOAAEvents" + str(stop),flavcut[setting] + " && " + stopcut)
  c1.SaveAs("try" + str(stop) + ".png")

  hOAAEventsAcc[stop] = TH1D("hOAAEventsAcc" + str(stop),"",len(bins[stop]) - 1,array('d',bins[stop]))
  chain.Draw("vtx[0]/-1.e2>>hOAAEventsAcc"+str(stop),"!(flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow || flagLepExitXHigh || flagLepExitXLow) && " + flavcut[setting] + " && " + stopcut)

  print stop, "Accepted", hOAAEventsAcc[stop].Integral()

  hOAAEventsExitAcc[stop] = TH1D("hOAAEventsExitAcc" + str(stop),"",len(bins[stop]) - 1,array('d',bins[stop]))
  chain.Draw("vtx[0]/-1.e2>>hOAAEventsExitAcc"+str(stop), exitcut + " && " + flavcut[setting] + " && " + stopcut)

  hOAAEventsMuHadAcc[stop] = TH1D("hOAAEventsMuHadAcc" + str(stop),"",len(bins[stop]) - 1,array('d',bins[stop]))
  chain.Draw("vtx[0]/-1.e2>>hOAAEventsMuHadAcc"+str(stop),"!(flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow || flagLepExitXHigh || flagLepExitXLow) && " + flavcut[setting] + " && (TotalNonlep_Dep_veto) <= " + str(veto) + " && " + stopcut)

  hOAAEventsExitMuHadAcc[stop] = TH1D("hOAAEventsExitMuHadAcc" + str(stop),"",len(bins[stop]) - 1,array('d',bins[stop]))
  chain.Draw("vtx[0]/-1.e2>>hOAAEventsExitMuHadAcc"+str(stop),exitcut + " && " + flavcut[setting] + " && (TotalNonlep_Dep_veto) <= " + str(veto) + " && " + stopcut)

  hOAAEventsEff[stop] = doEff(hOAAEvents[stop],hOAAEventsAcc[stop],"hOAAEventsEff"+str(stop))
  hOAAEventsMuHadEff[stop] = doEff(hOAAEvents[stop],hOAAEventsMuHadAcc[stop],"hOAAEventsMuHadEff"+str(stop))
  hOAAEventsExitEff[stop] = doEff(hOAAEvents[stop],hOAAEventsExitAcc[stop],"hOAAEventsExitEff"+str(stop))
  hOAAEventsExitMuHadEff[stop] = doEff(hOAAEvents[stop],hOAAEventsExitMuHadAcc[stop],"hOAAEventsExitMuHadEff"+str(stop))

for stop in stops:
#for stop in [0]:
  for i in range(1,hOAAEvents[stop].GetNbinsX()+1):
    full_bin = hOAAEventsFull.GetXaxis().FindBin(hOAAEvents[stop].GetBinCenter(i))

    hOAAEventsFull.SetBinContent(full_bin,hOAAEvents[stop].GetBinContent(i))
    hOAAEventsEffFull.SetBinContent(full_bin,hOAAEventsEff[stop].GetBinContent(i))    
    hOAAEventsMuHadEffFull.SetBinContent(full_bin,hOAAEventsMuHadEff[stop].GetBinContent(i))
    hOAAEventsExitEffFull.SetBinContent(full_bin,hOAAEventsExitEff[stop].GetBinContent(i))    
    hOAAEventsExitMuHadEffFull.SetBinContent(full_bin,hOAAEventsExitMuHadEff[stop].GetBinContent(i))

gPad.SetTicks(1,1)

hOAAEventsFull.SetTitle("Events - " + default_size + " - FV (50cm veto)")
hOAAEventsFull.SetXTitle("Off-axis position (m)")
hOAAEventsFull.GetXaxis().SetTitleSize(.06)
hOAAEventsFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_full.png")

hOAAEventsEffFull.SetTitle("Mu containment - " + default_size + " - FV (50cm veto)")
hOAAEventsEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsEffFull.SetYTitle("Efficiency")
hOAAEventsEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsEffFull.SetMaximum(1.0)
hOAAEventsEffFull.SetMinimum(0.0)
hOAAEventsEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_lepcontainedeff_full.png")

hOAAEventsMuHadEffFull.SetTitle("Mu containment + <" + str(veto*1000) + "MeV in veto - " + default_size + " - FV (50cm veto)")
hOAAEventsMuHadEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsMuHadEffFull.SetYTitle("Efficiency")
hOAAEventsMuHadEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsMuHadEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsMuHadEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsMuHadEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsMuHadEffFull.SetMaximum(1.0)
hOAAEventsMuHadEffFull.SetMinimum(0.0)
hOAAEventsMuHadEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_muHad_lepcontainedeff_full.png")

hOAAEventsExitEffFull.SetTitle("Mu exits KE > 50 MeV - " + default_size + " - FV (50cm veto)")
hOAAEventsExitEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsExitEffFull.SetYTitle("Efficiency")
hOAAEventsExitEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsExitEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsExitEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsExitEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsExitEffFull.SetMaximum(1.0)
hOAAEventsExitEffFull.SetMinimum(0.0)
hOAAEventsExitEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_lepexiteff_full.png")

hOAAEventsExitMuHadEffFull.SetTitle("Mu exits KE > 50 MeV + <" + str(veto*1000) +"MeV in veto - " + default_size + " - FV (50cm veto)")
hOAAEventsExitMuHadEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsExitMuHadEffFull.SetYTitle("Efficiency")
hOAAEventsExitMuHadEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsExitMuHadEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsExitMuHadEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsExitMuHadEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsExitMuHadEffFull.SetMaximum(1.0)
hOAAEventsExitMuHadEffFull.SetMinimum(0.0)
hOAAEventsExitMuHadEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_muhad_lepexiteff_full.png")

hOAAEnuEvents = dict()
hOAAEnuEventsAcc = dict()
hOAAEnuEventsEff = dict()
#hOAAEnuEventsFull = TH2D("hOAAEnuEventsFull","",50,0,5,n_bins,start_bin,end_bin)
#hOAAEnuEventsEffFull = TH2D("hOAAEnuEventsEffFull","",50,0,5,n_bins,start_bin,end_bin)
hOAAEnuEventsFull = TH2D("hOAAEnuEventsFull","",50,np.arange(0,5.1,.1),n_bins,array('d',full_bins))
hOAAEnuEventsEffFull = TH2D("hOAAEnuEventsEffFull","",50,np.arange(0,5.1,.1),n_bins,array('d',full_bins))

gStyle.SetPalette(54)

for stop in stops:
#for stop in [0]:

  stopcut = "stop == " + str(stop)

  hOAAEnuEvents[stop] = TH2D("hOAAEnuEvents" + str(stop), "",50,np.arange(0,5.1,.1), len(bins[stop]) - 1,array('d',bins[stop])) 
  chain.Draw("vtx[0]/-1.e2:Enu>>hOAAEnuEvents"+str(stop),flavcut[setting],"colz")

  hOAAEnuEventsAcc[stop] = TH2D("hOAAEnuEventsAcc" + str(stop), "",50,np.arange(0,5.1,.1), len(bins[stop]) - 1,array('d',bins[stop])) 
  chain.Draw("vtx[0]/-1.e2:Enu>>hOAAEnuEventsAcc"+str(stop),flavcut[setting] + " && !( flagLepExitBack || flagLepExitFront || flagLepExitYHigh || flagLepExitYLow|| flagLepExitXLow || flagLepExitXHigh)","colz")

  hOAAEnuEventsEff[stop] = doEff(hOAAEnuEvents[stop],hOAAEnuEventsAcc[stop],"hOAAEnuEventsEff"+str(stop))

for stop in stops:
#for stop in [0]:
  for i in range(1,hOAAEnuEvents[stop].GetNbinsX()+1):
    for j in range(1,hOAAEnuEvents[stop].GetNbinsY()+1):
#      full_bin = int((stop - min(stops.keys()))/10 + j)
      full_bin = hOAAEventsFull.GetXaxis().FindBin(hOAAEvents[stop].GetBinCenter(j))
      hOAAEnuEventsFull.SetBinContent(i,full_bin,hOAAEnuEvents[stop].GetBinContent(i,j))

      hOAAEnuEventsEffFull.SetBinContent(i,full_bin,hOAAEnuEventsEff[stop].GetBinContent(i,j))

hOAAEnuEventsFull.SetTitle("Events - " + default_size)
hOAAEnuEventsFull.SetXTitle("Enu (GeV)")
hOAAEnuEventsFull.SetYTitle("Off-axis position (m)")
hOAAEnuEventsFull.GetXaxis().SetTitleSize(.06)
hOAAEnuEventsFull.GetXaxis().SetTitleOffset(.6)
hOAAEnuEventsFull.GetYaxis().SetTitleSize(.06)
hOAAEnuEventsFull.GetYaxis().SetTitleOffset(.6)
hOAAEnuEventsFull.Draw("colz")
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_Enu_events_full.png")                     

hOAAEnuEventsEffFull.SetTitle("Mu containment - " + default_size)
hOAAEnuEventsEffFull.SetXTitle("Enu (GeV)")
hOAAEnuEventsEffFull.SetYTitle("Off-axis position (m)")
hOAAEnuEventsEffFull.GetXaxis().SetTitleSize(.06)
hOAAEnuEventsEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEnuEventsEffFull.GetYaxis().SetTitleSize(.06)
hOAAEnuEventsEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEnuEventsEffFull.Draw("colz")
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_Enu_events_eff_full.png")

#hOAAEnuExit = dict()
#for stop in stops:
#  hOAAEnuExit[stop] = TH2D("hOAAEnuExit" + str(stop), "", default_x*10,bins[stop][0],bins[stop][1]) 
#  chains[stop].Draw("vtx[0]/1.e2>>hOAAEnuexitLow"+str(stop),"LepPDG==13 && flagExitXLow","colz")
#  c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/exitingleplow"+str(stop)+".png")
#  chains[stop].Draw("vtx[0]/1.e2>>hOAAEnuexitHigh"+str(stop),"LepPDG==13 && flagExitXHigh","colz")
#  c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/exitinglephigh"+str(stop)+".png")
