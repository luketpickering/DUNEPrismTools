from ROOT import *
#from os import listdir
import os 
from array import *
from arg_parser import init_parser
import xml.etree.ElementTree
from doEff import doEff


args = init_parser().parse_args()
         
det_file = args.config

det_configs = []

print "Config file:\n\t",det_file 
config = xml.etree.ElementTree.parse(det_file)
for rp in config.findall('RunPlan'):
  stops = rp.findall('Stops')[0]
  det = rp.findall('Detector')[0]
  for stop in stops.findall('Stop'):
    det_configs.append(
      {'shift':int(stop.get('LateralOffset_m'))*-100,
      'x':int(float(det.get('DetectorFiducialWidth_m')))*100,
      'y':int(float(det.get('DetectorFiducialHeight_m')))*100,
      'z':int(float(det.get('DetectorFiducialDepth_m')))*100}
    )

print det_configs    
datadir = '/home/calcuttj/DUNEPrismSim/'
setting = args.setting
veto = args.veto


DIR = args.DIR
print "DIR:\n\t", DIR
output = os.listdir(datadir + setting + "/" +DIR)

if not os.path.isdir("OAA_events/"+setting+"/"+DIR):
  os.mkdir("OAA_events/"+setting+"/"+DIR)

trees = dict()
for dc in det_configs:
#  trees[dc['shift']] = "events_" + str(dc['x']) + "x" + str(dc['y']) + "x" + str(dc['z']) + "_50x50x50_"+str(dc['shift'])
   trees[dc['shift']] = "EDep_Stop" + str(dc['shift']/-100) + "_m"

default_size = str(dc['x']/100) + "x" + str(dc['y']/100) + "x" + str(dc['z']/100) +"m"
default_x = dc['x']/100

print trees

chains = dict()
stops = dict()
for dc in det_configs:
  chains[dc['shift']] = TChain(trees[dc['shift']])
  if dc['shift'] == 0:
    stops[dc['shift']] = "On Axis"
  else:
    stops[dc['shift']] = str(abs(dc['shift']/100)) + "m Off Axis"
print chains

for f in output:
  if len(f.split('.')) == 7:
    #print f
    for stop in chains:
      chains[stop].Add(datadir + setting + "/" + DIR + "/" + f)

total = 0
for dc in det_configs:
  total = total + chains[dc['shift']].GetEntries()
print total 
print "Loaded ",len(output), "files"

gStyle.SetOptStat(0)
gStyle.SetPalette(54)
c1 = TCanvas("c1","c1")

bins = dict()
for dc in det_configs:
  bins[-1*dc['shift']] = [-1*dc['shift']/100. - dc['x']/200., -1*dc['shift']/100. + dc['x']/200. ]

hOAAEvents = dict()
hOAAEventsAcc = dict()
hOAAEventsEff = dict()
hOAAEventsMuHadAcc = dict()
hOAAEventsMuHadEff = dict()
hOAAEventsExitAcc = dict()
hOAAEventsExitEff = dict()
hOAAEventsExitMuHadAcc = dict()
hOAAEventsExitMuHadEff = dict()

start_bin = (min(list(bins.keys()))/100. - default_x/2.)
end_bin = (max(list(bins.keys()))/100. + default_x/2.)

n_bins = int((end_bin - start_bin)*10)

print n_bins,start_bin,end_bin,default_x*10

hOAAEventsFull = TH1D("hOAAEventsFull","",n_bins,start_bin,end_bin)
hOAAEventsEffFull = TH1D("hOAAEventsEffFull","",n_bins,start_bin,end_bin)
hOAAEventsMuHadEffFull = TH1D("hOAAEventsMuHadEffFull","",n_bins,start_bin,end_bin)
hOAAEventsExitEffFull = TH1D("hOAAEventsExitEffFull","",n_bins,start_bin,end_bin)
hOAAEventsExitMuHadEffFull = TH1D("hOAAEventsExitMuHadEffFull","",n_bins,start_bin,end_bin)

exitcut = "(flagLepExitBack || flagLepExitXHigh || flagLepExitXLow || flagLepExitFront || flagLepExitY) && sqrt(lepExitingMomX*lepExitingMomX + lepExitingMomY*lepExitingMomY + lepExitingMomZ*lepExitingMomZ) > 0.114"

for stop in stops:
  hOAAEvents[stop] = TH1D("hOAAEvents" + str(stop), "", default_x*10,bins[-1*stop][0],bins[-1*stop][1])
  chains[stop].Draw("vtx[0]/-1.e2>>hOAAEvents" + str(stop),"LepPDG == 13 && NuPDG == 14 && !((flagLepExitXHigh || flagLepExitXLow) && sqrt(lepExitingMomX*lepExitingMomX + lepExitingMomY*lepExitingMomY + lepExitingMomZ*lepExitingMomZ) == 0)")

  hOAAEventsAcc[stop] = TH1D("hOAAEventsAcc" + str(stop),"",default_x*10,bins[-1*stop][0],bins[-1*stop][1])
  chains[stop].Draw("vtx[0]/-1.e2>>hOAAEventsAcc"+str(stop),"!(flagLepExitBack || flagLepExitFront || flagLepExitY || flagLepExitXHigh || flagLepExitXLow) && LepPDG == 13 && NuPDG == 14")

  hOAAEventsExitAcc[stop] = TH1D("hOAAEventsExitAcc" + str(stop),"",default_x*10,bins[-1*stop][0],bins[-1*stop][1])
  chains[stop].Draw("vtx[0]/-1.e2>>hOAAEventsExitAcc"+str(stop), exitcut + " && LepPDG == 13 && NuPDG == 14")

  hOAAEventsMuHadAcc[stop] = TH1D("hOAAEventsMuHadAcc" + str(stop),"",default_x*10,bins[-1*stop][0],bins[-1*stop][1])
  chains[stop].Draw("vtx[0]/-1.e2>>hOAAEventsMuHadAcc"+str(stop),"!(flagLepExitBack || flagLepExitFront || flagLepExitY || flagLepExitXHigh || flagLepExitXLow) && LepPDG == 13 && NuPDG == 14 && (TotalDep_veto + Pi0Dep_veto) <= " + str(veto) )

  hOAAEventsExitMuHadAcc[stop] = TH1D("hOAAEventsExitMuHadAcc" + str(stop),"",default_x*10,bins[-1*stop][0],bins[-1*stop][1])
  chains[stop].Draw("vtx[0]/-1.e2>>hOAAEventsExitMuHadAcc"+str(stop),exitcut + " && LepPDG == 13 && NuPDG == 14 && (TotalDep_veto + Pi0Dep_veto) <= " + str(veto) )

  hOAAEventsEff[stop] = doEff(hOAAEvents[stop],hOAAEventsAcc[stop],"hOAAEventsEff"+str(stop))
  hOAAEventsMuHadEff[stop] = doEff(hOAAEvents[stop],hOAAEventsMuHadAcc[stop],"hOAAEventsMuHadEff"+str(stop))
  hOAAEventsExitEff[stop] = doEff(hOAAEvents[stop],hOAAEventsExitAcc[stop],"hOAAEventsExitEff"+str(stop))
  hOAAEventsExitMuHadEff[stop] = doEff(hOAAEvents[stop],hOAAEventsExitMuHadAcc[stop],"hOAAEventsExitMuHadEff"+str(stop))

for stop in stops:
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
hOAAEventsEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_lepcontainedeff_full.png")

hOAAEventsMuHadEffFull.SetTitle("Mu containement + <" + str(veto*1000) + "MeV in veto - " + default_size + " - FV (50cm veto)")
hOAAEventsMuHadEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsMuHadEffFull.SetYTitle("Efficiency")
hOAAEventsMuHadEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsMuHadEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsMuHadEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsMuHadEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsMuHadEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_muHad_lepcontainedeff_full.png")

hOAAEventsExitEffFull.SetTitle("Mu exits KE > 50 MeV - " + default_size + " - FV (50cm veto)")
hOAAEventsExitEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsExitEffFull.SetYTitle("Efficiency")
hOAAEventsExitEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsExitEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsExitEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsExitEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsExitEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_lepexiteff_full.png")

hOAAEventsExitMuHadEffFull.SetTitle("Mu exits KE > 50 MeV + <" + str(veto*1000) +"MeV in veto - " + default_size + " - FV (50cm veto)")
hOAAEventsExitMuHadEffFull.SetXTitle("Off-axis position (m)")
hOAAEventsExitMuHadEffFull.SetYTitle("Efficiency")
hOAAEventsExitMuHadEffFull.GetXaxis().SetTitleSize(.06)
hOAAEventsExitMuHadEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEventsExitMuHadEffFull.GetYaxis().SetTitleSize(.06)
hOAAEventsExitMuHadEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEventsExitMuHadEffFull.Draw()
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_events_muhad_lepexiteff_full.png")

hOAAEnuEvents = dict()
hOAAEnuEventsAcc = dict()
hOAAEnuEventsEff = dict()
hOAAEnuEventsFull = TH2D("hOAAEnuEventsFull","",50,0,5,n_bins,start_bin,end_bin)
hOAAEnuEventsEffFull = TH2D("hOAAEnuEventsEffFull","",50,0,5,n_bins,start_bin,end_bin)

gStyle.SetPalette(54)

for stop in stops:
  hOAAEnuEvents[stop] = TH2D("hOAAEnuEvents" + str(stop), "",50,0.,5., default_x*10,bins[-1*stop][0],bins[-1*stop][1]) 
  chains[stop].Draw("vtx[0]/-1.e2:Enu>>hOAAEnuEvents"+str(stop),"LepPDG==13","colz")

  hOAAEnuEventsAcc[stop] = TH2D("hOAAEnuEventsAcc" + str(stop), "",50,0.,5., default_x*10,bins[-1*stop][0],bins[-1*stop][1]) 
  chains[stop].Draw("vtx[0]/-1.e2:Enu>>hOAAEnuEventsAcc"+str(stop),"LepPDG==13 && !( flagLepExitBack || flagLepExitFront || flagLepExitY|| flagLepExitXLow || flagLepExitXHigh)","colz")

  hOAAEnuEventsEff[stop] = doEff(hOAAEnuEvents[stop],hOAAEnuEventsAcc[stop],"hOAAEnuEventsEff"+str(stop))

for stop in stops:
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
