from ROOT import *
#from os import listdir, mkdir, path.isdir
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

DIR = args.DIR
output = os.listdir(datadir + setting + "/"+DIR)

print "DIR:\n\t", DIR
if not os.path.isdir("ElastEnuEff/"+setting+"/"+DIR):
  os.mkdir("ElastEnuEff/"+setting+"/"+DIR)

trees = dict()
for dc in det_configs:
#  trees[dc['shift']] = "events_" + str(dc['x']) + "x" + str(dc['y']) + "x" + str(dc['z']) + "_50x50x50_"+str(dc['shift'])
   trees[dc['shift']] = "EDep_Stop" + str(dc['shift']/-100) + "_m"

default_size = str(dc['x']/100) + "x" + str(dc['y']/100) + "x" + str(dc['z']/100) +"m"

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
#    print f
    for stop in chains:
      chains[stop].Add(datadir + setting + "/" +DIR+"/"+ f)

print chains[0].GetEntries()

gStyle.SetOptStat(0)
gStyle.SetPalette(54)
c1 = TCanvas("c1","c1")
c1.SetLogz()

a = 0
hAcc = {} 
hTot = {} 
hEff = {} 

c1.SetLogz()


nXBins = 11
xBins = array('d',[0.,1.,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0])
nYBins = 10
yBins = array('d',[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0])

hElastEnuTot = dict();
hElastEnuAcc = dict();
hElastEnuEff = dict();

for stop in stops:
  c1.SetLogz()
  hElastEnuTot[stop] = TH2D("hElastEnuTot_" + str(stop),"",nXBins,xBins,nYBins,yBins)
  chains[stop].Draw("(1-yTrue):Enu>>hElastEnuTot_"+str(stop),"","colz")
  hElastEnuTot[stop].SetTitle("Events - " + default_size + " - " +stops[stop])
  hElastEnuTot[stop].SetXTitle("Enu (GeV)")
  hElastEnuTot[stop].SetYTitle("Elasticity (1-y)")
  hElastEnuTot[stop].GetXaxis().SetTitleSize(.06)
  hElastEnuTot[stop].GetYaxis().SetTitleSize(.06)
  hElastEnuTot[stop].GetXaxis().SetTitleOffset(.6)
  hElastEnuTot[stop].GetYaxis().SetTitleOffset(.6)
  hElastEnuTot[stop].Draw("colz")
  c1.SaveAs("ElastEnuEff/"+setting+ "/"+DIR+ "/Tot"+str(stop)+".png")


  hElastEnuAcc[stop] = TH2D("hElastEnuAcc_" + str(stop),"",nXBins,xBins,nYBins,yBins)
  chains[stop].Draw("(1-yTrue):Enu>>hElastEnuAcc_"+str(stop),"TotalDep_veto < 0.05","colz")
  c1.SaveAs("ElastEnuEff/"+setting+ "/"+DIR+ "/Acc"+str(stop)+".png")

  c1.SetLogz(0)
  hElastEnuEff[stop] = doEff(hElastEnuTot[stop],hElastEnuAcc[stop],"hElastEnuEff_"+str(stop)) 
  hElastEnuEff[stop].SetMaximum(1.0)
  hElastEnuEff[stop].SetMinimum(0.)
  hElastEnuEff[stop].SetTitle(" Had Containment - " + default_size + " - " + stops[stop])
  hElastEnuEff[stop].SetXTitle("Enu (GeV)")
  hElastEnuEff[stop].SetYTitle("Elasticity (1-y)")
  hElastEnuEff[stop].GetXaxis().SetTitleSize(.06)
  hElastEnuEff[stop].GetYaxis().SetTitleSize(.06)
  hElastEnuEff[stop].GetXaxis().SetTitleOffset(.6)
  hElastEnuEff[stop].GetYaxis().SetTitleOffset(.6)
  hElastEnuEff[stop].Draw("colz")
  c1.SaveAs("ElastEnuEff/"+setting+ "/"+DIR+"/Eff"+str(stop)+".png")

