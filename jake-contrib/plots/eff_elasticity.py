from ROOT import *
#from os import listdir, mkdir, path.isdir
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
DIR = DIR + args.extra

flavcut = {"FHC":"LepPDG == 13 && NuPDG == 14",
           "RHC":"LepPDG == -13 && NuPDG == -14"}
extracut = "1"
for stop in stops:

  stopcut = "stop == " + str(stop)

  c1.SetLogz()
  hElastEnuTot[stop] = TH2D("hElastEnuTot_" + str(stop),"",nXBins,xBins,nYBins,yBins)
  chain.Draw("(1-yTrue):Enu>>hElastEnuTot_"+str(stop),stopcut + " && " + flavcut[setting] ,"colz")
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
  chain.Draw("(1-yTrue):Enu>>hElastEnuAcc_"+str(stop),stopcut + " && " + extracut + " && (TotalNonlep_Dep_veto) <= " + str(veto),"colz")

  c1.SetLogz(0)
  hElastEnuEff[stop] = doEff(hElastEnuTot[stop],hElastEnuAcc[stop],"hElastEnuEff_"+str(stop)) 
  hElastEnuEff[stop].SetMaximum(1.0)
  hElastEnuEff[stop].SetMinimum(0.)
  hElastEnuEff[stop].SetTitle(str(veto) + " had energy veto - " + default_size + " - " + stops[stop])
  hElastEnuEff[stop].SetXTitle("Enu (GeV)")
  hElastEnuEff[stop].SetYTitle("Elasticity (1-y)")
  hElastEnuEff[stop].GetXaxis().SetTitleSize(.06)
  hElastEnuEff[stop].GetYaxis().SetTitleSize(.06)
  hElastEnuEff[stop].GetXaxis().SetTitleOffset(.6)
  hElastEnuEff[stop].GetYaxis().SetTitleOffset(.6)
  hElastEnuEff[stop].Draw("colz")
  c1.SaveAs("ElastEnuEff/"+setting+ "/"+DIR+"/Eff"+str(stop)+".png")

