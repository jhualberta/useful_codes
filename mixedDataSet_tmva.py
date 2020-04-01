import ROOT
import sys, os, getopt
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
path = os.getcwd()
file_list = os.listdir(path)
### signals
inputSignal = "MCsolar_200004to206391.root"
f0 = TFile(inputSignal)
treeSolar = f0.Get("TsolarMC")
totalEvents = treeSolar.GetEntries()
posx = array('d',[0]);posy = array('d',[0]);posz = array('d',[0]);dirx = array('d',[0]);diry = array('d',[0]);dirz = array('d',[0])
itr = array('d',[0]);
posFOM = array('d',[0]);posFOM2 = array('d',[0])
nhits = array('d',[0]);beta14 = array('d',[0]);thetaij = array('d',[0])
treeSolar.SetBranchAddress("posx", posx)
treeSolar.SetBranchAddress("posy",posy)
treeSolar.SetBranchAddress("posz", posz)
treeSolar.SetBranchAddress("dirx", dirx)
treeSolar.SetBranchAddress("diry", diry)
treeSolar.SetBranchAddress("dirz", dirz)

treeSolar.SetBranchAddress("itr", itr)
treeSolar.SetBranchAddress("posFOM", posFOM)
treeSolar.SetBranchAddress("posFOM2",posFOM2)
treeSolar.SetBranchAddress("nhits", nhits)
treeSolar.SetBranchAddress("beta14", beta14)
treeSolar.SetBranchAddress("thetaij", thetaij)
### backgrounds
inputBkg = "MCTl208_500files_200004to200893.root"
f1 = TFile(inputBkg)
treeTl208 = f1.Get("TsolarMC")
totalEvents = treeSolar.GetEntries()
treeTl208.SetBranchAddress("posx", posx)
treeTl208.SetBranchAddress("posy",posy)
treeTl208.SetBranchAddress("posz", posz)
treeTl208.SetBranchAddress("dirx", dirx)
treeTl208.SetBranchAddress("diry", diry)
treeTl208.SetBranchAddress("dirz", dirz)

treeTl208.SetBranchAddress("itr", itr)
treeTl208.SetBranchAddress("posFOM", posFOM)
treeTl208.SetBranchAddress("posFOM2",posFOM2)
treeTl208.SetBranchAddress("nhits", nhits)
treeTl208.SetBranchAddress("beta14", beta14)
treeTl208.SetBranchAddress("thetaij", thetaij)

data_posx = []
data_posy = []
data_posz = []
data_dirx = []
data_diry = []
data_dirz = []
data_itr = []
data_beta14 = []
data_thetaij = []
data_nhits = []
data_fom = []
data_fom2 = []

### Must create a root file, otherwise it complains ntuple filling
fff = TFile("MixedDataSet_ntuple.root","recreate")

# create a TNtuple
ntSolar = ROOT.TNtuple("ntupleMixed","solar+bkg","itr:beta14:thetaij:nhits:udotR:logL:scaleLogL:signal")
#
countSignal = 0
countBkg = 0

# generate 'signal' and 'background' distributions
# !!!!! signals
for i in range(0, totalEvents):
    treeSolar.GetEntry(i)
    #data_posx.append(posx[0])
    #data_posy.append(posy[0]) 
    #data_posz.append(posz[0])
    #data_dirx.append(dirx[0])
    #data_diry.append(diry[0])
    #data_dirz.append(dirz[0])
    #data_beta14.append(beta14[0])
    #data_thetaij.append(thetaij[0])
    #data_nhits.append(nhits[0])
    #data_fom.append(posFOM[0])
    #data_fom2.append(posFOM2[0])
    udotR = (posx[0]*dirx[0]+posy[0]*diry[0]+(posz[0]-108)*dirz[0])/(posx[0]*posx[0]+posy[0]*posy[0]+(posz[0]-108)*(posz[0]-108))
    scaleLogL = posFOM[0]/posFOM2[0]
    ntSolar.Fill(itr[0], beta14[0],thetaij[0],nhits[0], udotR, posFOM[0], scaleLogL, 1)
    countSignal += 1

totalBkg = treeTl208.GetEntries()
for i in range(0,totalBkg): 
    treeTl208.GetEntry(i)
    #data_posx.append(posx[0])
    #data_posy.append(posy[0])
    #data_posz.append(posz[0])
    #data_dirx.append(dirx[0])
    #data_diry.append(diry[0])
    #data_dirz.append(dirz[0])
    #data_beta14.append(beta14[0])
    #data_thetaij.append(thetaij[0])
    #data_nhits.append(nhits[0])
    #data_fom.append(posFOM[0])
    #data_fom2.append(posFOM2[0])
    udotR = (posx[0]*dirx[0]+posy[0]*diry[0]+(posz[0]-108)*dirz[0])/(posx[0]*posx[0]+posy[0]*posy[0]+(posz[0]-108)*(posz[0]-108))
    scaleLogL = posFOM[0]/posFOM2[0]
    ntSolar.Fill(itr[0], beta14[0],thetaij[0],nhits[0], udotR, posFOM[0], scaleLogL, -1)
    countBkg += 1


print "signals", countSignal, "backgrounds", countBkg

fff.cd()

ntSolar.Write()

fff.Close()
