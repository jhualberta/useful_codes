import ROOT
import sys, os, getopt
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
from numpy import sqrt
import csv

from array import array
path = os.getcwd()
file_list = os.listdir(path)
### signals
inputSignal = "Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to200658_s0_p0.root"
f0 = TFile(inputSignal)
treeSolar = f0.Get("T")
totalEvents = treeSolar.GetEntries()
posx = array('d',[0]);posy = array('d',[0]);posz = array('d',[0]);dirx = array('d',[0]);diry = array('d',[0]);dirz = array('d',[0])
itr = array('d',[0]);
posFOM = array('d',[0]);posFOM2 = array('d',[0])
nhits = array('I',[0]);beta14 = array('d',[0]);thetaij = array('d',[0])
energy = array('d',[0]) # use corrected energy
Utest = array('d',[0]); Gtest = array('d',[0]); scaleLogL = array('d',[0])
medianDevHit = array('d',[0]);medianDev = array('d',[0]);medianProbHit = array('d',[0]);medianProb = array('d',[0]);

posxA = array('d',[0]);posyA = array('d',[0]);poszA = array('d',[0]);dirxA = array('d',[0]);diryA = array('d',[0]);dirzA = array('d',[0])
itrA = array('d',[0]);
posFOMA = array('d',[0]);posFOM2A = array('d',[0])
nhitsA = array('I',[0]);beta14A = array('d',[0]);thetaijA = array('d',[0])
energyA = array('d',[0]) # use corrected energy
UtestA = array('d',[0]); GtestA = array('d',[0]); scaleLogLA = array('d',[0])
medianDevHitA = array('d',[0]);medianDevA = array('d',[0]);medianProbHitA = array('d',[0]);medianProbA= array('d',[0]);

posxB = array('d',[0]);posyB = array('d',[0]);poszB = array('d',[0]);dirxB = array('d',[0]);diryB = array('d',[0]);dirzB = array('d',[0])
itrB = array('d',[0]);
posFOMB = array('d',[0]);posFOM2B = array('d',[0])
nhitsB = array('I',[0]);beta14B = array('d',[0]);thetaijB = array('d',[0])
energyB = array('d',[0]) # use corrected energy
UtestB = array('d',[0]); GtestB = array('d',[0]); scaleLogLB = array('d',[0])
medianDevHitB = array('d',[0]);medianDevB = array('d',[0]);medianProbHitB = array('d',[0]);medianProbB= array('d',[0]);

treeSolar.SetBranchAddress("posx", posx)
treeSolar.SetBranchAddress("posy",posy)
treeSolar.SetBranchAddress("posz", posz)
treeSolar.SetBranchAddress("dirx", dirx)
treeSolar.SetBranchAddress("diry", diry)
treeSolar.SetBranchAddress("dirz", dirz)
treeSolar.SetBranchAddress("itr", itr)
treeSolar.SetBranchAddress("nhits", nhits)
treeSolar.SetBranchAddress("beta14", beta14)
treeSolar.SetBranchAddress("thetaij", thetaij)
treeSolar.SetBranchAddress("energy", energy) # use corrected energy
treeSolar.SetBranchAddress("Gtest", Gtest) # use corrected energy
treeSolar.SetBranchAddress("Utest", Utest) # use corrected energy
treeSolar.SetBranchAddress("scaleLogL", scaleLogL) # use corrected energy
treeSolar.SetBranchAddress("medianDevHit", medianDevHit)
treeSolar.SetBranchAddress("medianDev", medianDev)
treeSolar.SetBranchAddress("medianProbHit", medianProbHit)
treeSolar.SetBranchAddress("medianProb", medianProb)

### backgrounds, Bi214 internal
inputBkg = "Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214Run_r200004to200658_s0_p0.root"
f1 = TFile(inputBkg)
treeBi214 = f1.Get("T")
totalEvents = treeSolar.GetEntries()
treeBi214.SetBranchAddress("posx", posxA)
treeBi214.SetBranchAddress("posy",posyA)
treeBi214.SetBranchAddress("posz", poszA)
treeBi214.SetBranchAddress("dirx", dirxA)
treeBi214.SetBranchAddress("diry", diryA)
treeBi214.SetBranchAddress("dirz", dirzA)
treeBi214.SetBranchAddress("itr", itrA)
treeBi214.SetBranchAddress("nhits", nhitsA)
treeBi214.SetBranchAddress("beta14", beta14A)
treeBi214.SetBranchAddress("thetaij", thetaijA)
treeBi214.SetBranchAddress("energy", energyA) # use corrected energy
treeBi214.SetBranchAddress("Gtest", GtestA) # use corrected energy
treeBi214.SetBranchAddress("Utest", UtestA) # use corrected energy
treeBi214.SetBranchAddress("scaleLogL", scaleLogLA) # use corrected energy
treeBi214.SetBranchAddress("medianDevHit", medianDevHitA)
treeBi214.SetBranchAddress("medianDev", medianDevA)
treeBi214.SetBranchAddress("medianProbHit", medianProbHitA)
treeBi214.SetBranchAddress("medianProb", medianProbA)
### backgrounds, Tl208 internal
inputBkg1 = "Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208Run_r200004to200658_s0_p0.root"
f2 = TFile(inputBkg1)
treeTl208 = f2.Get("T")
totalEvents = treeSolar.GetEntries()
treeTl208.SetBranchAddress("posx", posxB)
treeTl208.SetBranchAddress("posy",posyB)
treeTl208.SetBranchAddress("posz", poszB)
treeTl208.SetBranchAddress("dirx", dirxB)
treeTl208.SetBranchAddress("diry", diryB)
treeTl208.SetBranchAddress("dirz", dirzB)
treeTl208.SetBranchAddress("itr", itrB)
treeTl208.SetBranchAddress("nhits", nhitsB)
treeTl208.SetBranchAddress("beta14", beta14B)
treeTl208.SetBranchAddress("thetaij", thetaijB)
treeTl208.SetBranchAddress("energy", energyB) # use corrected energy
treeTl208.SetBranchAddress("Gtest", GtestB) # use corrected energy
treeTl208.SetBranchAddress("Utest", UtestB) # use corrected energy
treeTl208.SetBranchAddress("scaleLogL", scaleLogLB) # use corrected energy
treeTl208.SetBranchAddress("medianDevHit", medianDevHitB)
treeTl208.SetBranchAddress("medianDev", medianDevB)
treeTl208.SetBranchAddress("medianProbHit", medianProbHitB)
treeTl208.SetBranchAddress("medianProb", medianProbB)
### Must create a root file, otherwise it complains ntuple filling
fff = TFile("MixedDataSet_MP_r200004to200658_4to5MeV.root","recreate")

# create a TNtuple
ntSolarMix = ROOT.TNtuple("ntupleMixed","solar+bkg","itr:beta14:thetaij:nhits:energy:udotR:Gtest:Utest:scaleLogL:zFactor:signal")
#
countSignal = 0
countBkg = 0

# generate 'signal' and 'background' distributions
# !!!!! signals
for i in range(0, totalEvents):
    treeSolar.GetEntry(i)
    if sqrt(posx[0]*posx[0]+posy[0]*posy[0]+(posz[0]-108)*(posz[0]-108))<6000 and nhits[0]>14 and energy[0]<5 and energy[0]>4:
      udotR = (posx[0]*dirx[0]+posy[0]*diry[0]+(posz[0]-108)*dirz[0])/(posx[0]*posx[0]+posy[0]*posy[0]+(posz[0]-108)*(posz[0]-108))
      zFactor = 1#1-3*(medianDevHit[0]+medianDev[0])/(medianProbHit[0]-medianProb[0])
      ntSolarMix.Fill(itr[0], beta14[0],thetaij[0],nhits[0], energy[0], udotR, Gtest[0], Utest[0], scaleLogL[0], zFactor, 1)
      countSignal += 1

totalBkg1 = treeBi214.GetEntries()
for i in range(0,totalBkg1):
    treeBi214.GetEntry(i)
    if sqrt(posxA[0]*posxA[0]+posyA[0]*posyA[0]+(poszA[0]-108)*(poszA[0]-108))<6000 and nhitsA[0]>14 and energyA[0]<5 and energyA[0]>4:
      udotR = (posxA[0]*dirxA[0]+posyA[0]*diryA[0]+(poszA[0]-108)*dirzA[0])/(posxA[0]*posxA[0]+posyA[0]*posyA[0]+(poszA[0]-108)*(poszA[0]-108))
      zFactor = 1#-3*(medianDevHitA[0]+medianDevA[0])/(medianProbHitA[0]-medianProbA[0])
      ntSolarMix.Fill(itrA[0], beta14A[0],thetaijA[0],nhitsA[0], energyA[0], udotR, GtestA[0], UtestA[0], scaleLogLA[0], zFactor,-1)
      countBkg += 1

totalBkg2 = treeTl208.GetEntries()
for i in range(0,totalBkg2): 
    treeTl208.GetEntry(i)
    if sqrt(posxB[0]*posxB[0]+posyB[0]*posyB[0]+(poszB[0]-108)*(poszB[0]-108))<6000 and nhitsB[0]>14 and energyB[0]<5 and energyB[0]>4:
      udotR = (posxB[0]*dirxB[0]+posyB[0]*diryB[0]+(poszB[0]-108)*dirzB[0])/(posxB[0]*posxB[0]+posyB[0]*posyB[0]+(poszB[0]-108)*(poszB[0]-108))
      zFactor = 1#-3*(medianDevHitB[0]+medianDevB[0])/(medianProbHitB[0]-medianProbB[0])
      ntSolarMix.Fill(itrB[0], beta14B[0],thetaijB[0],nhitsB[0], energyB[0], udotR, GtestB[0], UtestB[0], scaleLogLB[0], zFactor,-1)
      countBkg += 1

print "signals", countSignal, "backgrounds", countBkg

fff.cd()

ntSolarMix.Write()

fff.Close()
