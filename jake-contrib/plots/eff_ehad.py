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
  for det in rp.findall('Detector'):
    det_configs.append(
      {'shift':int(det.get('shift'))*100,
      'x':int(float(det.get('detectorSizeX')))*100,
      'y':int(float(det.get('detectorSizeY')))*100,
      'z':int(float(det.get('detectorSizeZ')))*100}
    )

print det_configs    
datadir = '/home/calcuttj/DUNEPrismSim/'
setting = args.setting



DIR = args.DIR
output = os.listdir(datadir + setting + "/"+DIR)

print "DIR:\n\t", DIR
if not os.path.isdir("OAA_events/"+setting+"/"+DIR):
  os.mkdir("OAA_events/"+setting+"/"+DIR)

trees = dict()
for dc in det_configs:
  trees[dc['shift']] = "events_" + str(dc['x']) + "x" + str(dc['y']) + "x" + str(dc['z']) + "_50x50x50_"+str(dc['shift'])

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

for f in output:
  if len(f.split('.')) == 7:
    print f
    for stop in chains:
      chains[stop].Add(datadir + setting + "/" +DIR+"/"+ f)

print chains[0].GetEntries()

gStyle.SetOptStat(0)
gStyle.SetPalette(54)
c1 = TCanvas("c1","c1")

bins = dict()
for dc in det_configs:
  bins[dc['shift']] = [dc['shift']/100. - dc['x']/200., dc['shift']/100. + dc['x']/200. ]

hOAAEvents = dict()
hOAAEventsAcc = dict()
hOAAEventsEff = dict()
hOAAEventsMuHadAcc = dict()
hOAAEventsMuHadEff = dict()

start_bin = (min(list(bins.keys()))/100. - default_x/2.)
end_bin = (max(list(bins.keys()))/100. + default_x/2.)

n_bins = int((end_bin - start_bin)*10)

print n_bins,start_bin,end_bin,default_x*10

hEHadEvents = dict()
hEHadEventsAcc = dict()
hEHadEventsHadAcc = dict()
hEHadEventsEff = dict()
hEHadEventsHadEff = dict()
hOAAEHadEvents = dict()
hOAAEHadEventsAcc = dict()
hOAAEHadEventsEff = dict()
hOAAEHadEventsHadAcc = dict()
hOAAEHadEventsHadEff = dict()
hOAAEHadEventsFull = TH2D("hOAAEHadEventsFull","",50,0,5,n_bins,start_bin,end_bin)
hOAAEHadEventsEffFull = TH2D("hOAAEHadEventsEffFull","",50,0,5,n_bins,start_bin,end_bin)
hOAAEHadEventsHadEffFull = TH2D("hOAAEHadEventsHadEffFull","",50,0,5,n_bins,start_bin,end_bin)

gStyle.SetPalette(54)

for stop in stops:
  #1D Section
  hEHadEvents[stop] = TH1D("hEHadEvents" + str(stop), "",50,0.,5.) 
  chains[stop].Draw("eHadTrueTotal>>hEHadEvents"+str(stop),"lepPDG==13","colz")

  hEHadEventsAcc[stop] = TH1D("hEHadEventsAcc" + str(stop), "",50,0.,5.) 
  chains[stop].Draw("eHadTrueTotal>>hEHadEventsAcc"+str(stop),"lepPDG==13 && flagLepContained","colz")

  hEHadEventsHadAcc[stop] = TH1D("hEHadEventsHadAcc" + str(stop), "",50,0.,5.) 
  chains[stop].Draw("eHadTrueTotal>>hEHadEventsHadAcc"+str(stop),"lepPDG==13 && flagLepContained && (eHadPrimaryDepIn + eHadSecondaryDepIn)/(eHadPrimaryDepIn + eHadSecondaryDepIn + eHadPrimaryDepOut + eHadSecondaryDepOut) >= .95","colz")

  hEHadEventsEff[stop] = doEff(hEHadEvents[stop],hEHadEventsAcc[stop],"hEHadEventsEff"+str(stop))
  hEHadEventsHadEff[stop] = doEff(hEHadEvents[stop],hEHadEventsHadAcc[stop],"hEHadEventsHadEff"+str(stop))

  hEHadEventsEff[stop].SetTitle("Mu containment - " + default_size)
  hEHadEventsEff[stop].SetXTitle("EHad (GeV)")
  hEHadEventsEff[stop].GetXaxis().SetTitleSize(.06)
  hEHadEventsEff[stop].GetXaxis().SetTitleOffset(.6)
  hEHadEventsEff[stop].Draw("colz")
#  c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/EHad_events_eff_" + str(stop) + ".png")

  hEHadEventsHadEff[stop].SetTitle("Mu + .95 EHad containment - " + default_size)
  hEHadEventsHadEff[stop].SetXTitle("EHad (GeV)")
  hEHadEventsHadEff[stop].GetXaxis().SetTitleSize(.06)
  hEHadEventsHadEff[stop].GetXaxis().SetTitleOffset(.6)
  hEHadEventsHadEff[stop].Draw("colz")
#  c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/EHad_events_had_eff_" + str(stop) + ".png")

  #2D Section
  hOAAEHadEvents[stop] = TH2D("hOAAEHadEvents" + str(stop), "",50,0.,5., default_x*10,bins[stop][0],bins[stop][1]) 
  chains[stop].Draw("vtx_X/1.e2:eHadTrueTotal>>hOAAEHadEvents"+str(stop),"lepPDG==13","colz")

  hOAAEHadEventsAcc[stop] = TH2D("hOAAEHadEventsAcc" + str(stop), "",50,0.,5., default_x*10,bins[stop][0],bins[stop][1]) 
  chains[stop].Draw("vtx_X/1.e2:eHadTrueTotal>>hOAAEHadEventsAcc"+str(stop),"lepPDG==13 && flagLepContained","colz")

  hOAAEHadEventsHadAcc[stop] = TH2D("hOAAEHadEventsHadAcc" + str(stop), "",50,0.,5., default_x*10,bins[stop][0],bins[stop][1]) 
  chains[stop].Draw("vtx_X/1.e2:eHadTrueTotal>>hOAAEHadEventsHadAcc"+str(stop),"lepPDG==13 && flagLepContained && (eHadPrimaryDepIn + eHadSecondaryDepIn)/(eHadPrimaryDepIn + eHadSecondaryDepIn + eHadPrimaryDepOut + eHadSecondaryDepOut) >= .95","colz")

  hOAAEHadEventsEff[stop] = doEff(hOAAEHadEvents[stop],hOAAEHadEventsAcc[stop],"hOAAEHadEventsEff"+str(stop))
  hOAAEHadEventsHadEff[stop] = doEff(hOAAEHadEvents[stop],hOAAEHadEventsHadAcc[stop],"hOAAEHadEventsHadEff"+str(stop))

for stop in stops:
  for i in range(1,hOAAEHadEvents[stop].GetNbinsX()+1):
    for j in range(1,hOAAEHadEvents[stop].GetNbinsY()+1):
      full_bin = int((stop - min(stops.keys()))/10 + j)
      hOAAEHadEventsFull.SetBinContent(i,full_bin,hOAAEHadEvents[stop].GetBinContent(i,j))
      hOAAEHadEventsEffFull.SetBinContent(i,full_bin,hOAAEHadEventsEff[stop].GetBinContent(i,j))
      hOAAEHadEventsHadEffFull.SetBinContent(i,full_bin,hOAAEHadEventsHadEff[stop].GetBinContent(i,j))

hOAAEHadEventsFull.SetTitle("Events - " + default_size)
hOAAEHadEventsFull.SetXTitle("EHad (GeV)")
hOAAEHadEventsFull.SetYTitle("Off-axis position (m)")
hOAAEHadEventsFull.GetXaxis().SetTitleSize(.06)
hOAAEHadEventsFull.GetXaxis().SetTitleOffset(.6)
hOAAEHadEventsFull.GetYaxis().SetTitleSize(.06)
hOAAEHadEventsFull.GetYaxis().SetTitleOffset(.6)
hOAAEHadEventsFull.Draw("colz")
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_EHad_events_full.png")

hOAAEHadEventsEffFull.SetTitle("Mu containment - " + default_size)
hOAAEHadEventsEffFull.SetXTitle("EHad (GeV)")
hOAAEHadEventsEffFull.SetYTitle("Off-axis position (m)")
hOAAEHadEventsEffFull.GetXaxis().SetTitleSize(.06)
hOAAEHadEventsEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEHadEventsEffFull.GetYaxis().SetTitleSize(.06)
hOAAEHadEventsEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEHadEventsEffFull.Draw("colz")
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_EHad_events_eff_full.png")

hOAAEHadEventsHadEffFull.SetTitle("Mu + .95 EHad containment - " + default_size)
hOAAEHadEventsHadEffFull.SetXTitle("EHad (GeV)")
hOAAEHadEventsHadEffFull.SetYTitle("Off-axis position (m)")
hOAAEHadEventsHadEffFull.GetXaxis().SetTitleSize(.06)
hOAAEHadEventsHadEffFull.GetXaxis().SetTitleOffset(.6)
hOAAEHadEventsHadEffFull.GetYaxis().SetTitleSize(.06)
hOAAEHadEventsHadEffFull.GetYaxis().SetTitleOffset(.6)
hOAAEHadEventsHadEffFull.Draw("colz")
c1.SaveAs("OAA_events/"+setting+"/"+DIR+"/OAA_EHad_events_had_eff_full.png")

