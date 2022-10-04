import ROOT
import sys, os, getopt
from ROOT import *
from rat import *
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
r = 850.0 # AV radius

checkID = 1155751 
#fileList = ["Mergedtree_030708.root"]#Mergedtree_030747.root"]
fileList = ["MergedNew_9runs.root"]
#fileList = ["Mergedtree_030708.root"]
#fileList = ["Mergedtree_030717.root"]
#fileList = ["Mergedtree_030726.root"]
#fileList = ["Mergedtree_030741.root"]
#fileList = ["Mergedtree_030747.root"]
#fileList = ["Mergedtree_030756.root"]
#fileList = ["Mergedtree_030765.root"]
#fileList = ["Mergedtree_030774.root"]

#for files in fileList:
fname = fileList[0]
#fnew = TFile("saveClusters_"+fname[12:17]+".root","recreate")
fnew = TFile("saveClusters_all.root","recreate")
# ==============================
# 2D phi,theta; phi, theta Q weight; phi, theta T weight
hPhiCosTheta_singleCluster = TH2F("hPhiCosTheta_singleCluster",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Qweighted_singleCluster = TH2F("hPhiCosTheta_Qweighted_singleCluster", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Tweighted_singleCluster = TH2F("hPhiCosTheta_Tweighted_singleCluster", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

# 1D phi ; theta; phi; theta Q weight; phi; theta T weight
hPhi_singleCluster = TH1F("hPhi_singleCluster",";#phi(rad);" ,100, 0, 6)
hCosTheta_singleCluster = TH1F("hCosTheta_singleCluster",";Cos(#theta)" , 200, -1, 1)
hPhi_Qweighted_singleCluster = TH1F("hPhi_Qweighted_singleCluster", "charge;#phi(rad);" ,100, 0, 6)
hCosTheta_Qweighted_singleCluster = TH1F("hCosTheta_Qweighted_singleCluster", "charge;Cos(#theta);" , 200, -1, 1)
hPhi_Tweighted_singleCluster = TH1F("hPhi_Tweighted_singleCluster", "hitTime;#phi(rad);" ,100, 0, 6)
hCosTheta_Tweighted_singleCluster = TH1F("hCosTheta_Tweighted_singleCluster", "hitTime;Cos(#theta);" , 200, -1, 1)

# average case
hAvePhiCosTheta_singleCluster = TH2F("hAvePhiCosTheta_singleCluster",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Qweighted_singleCluster = TH2F("hAvePhiCosTheta_Qweighted_singleCluster","charge;#phi(rad);" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tweighted_singleCluster = TH2F("hAvePhiCosTheta_Tweighted_singleCluster","hitTime;Cos(#theta);" ,100, 0, 6, 200, -1, 1)

hAvePhiCosTheta_Qfraction_singleCluster = TH2F("hAvePhiCosTheta_Qfraction_single",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tfraction_singleCluster = TH2F("hAvePhiCosTheta_Tfraction_single",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

hAvePhi_singleCluster = TH1F("hAvePhi_singleCluster",";#phi(rad);" ,100, 0, 6)
hAveCosTheta_singleCluster = TH1F("hAveCosTheta_singleCluster",";Cos(#theta);" , 200, -1, 1)
hAvePhi_Qweighted_singleCluster = TH1F("hAvePhi_Qweighted_singleCluster","charge;#phi(rad);" ,100, 0, 6)
hAvePhi_Tweighted_singleCluster = TH1F("hAvePhi_Tweighted_singleCluster","hitTime;#phi(rad);" ,100, 0, 6)
hAveCosTheta_Qweighted_singleCluster = TH1F("hAveCosTheta_Qweighted_singleCluster","charge;Cos(#theta);" , 200, -1, 1)
hAveCosTheta_Tweighted_singleCluster = TH1F("hAveCosTheta_Tweighted_singleCluster","hitTime;Cos(#theta);" , 200, -1, 1)
# ==============================
# 2D phi,theta; phi, theta Q weight; phi, theta T weight
hPhiCosTheta_doubleA = TH2F("hPhiCosTheta_doubleA",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Qweighted_doubleA = TH2F("hPhiCosTheta_Qweighted_doubleA", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Tweighted_doubleA = TH2F("hPhiCosTheta_Tweighted_doubleA", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

# 1D phi ; theta; phi; theta Q weight; phi; theta T weight
hPhi_doubleA = TH1F("hPhi_doubleA",";#phi(rad);" ,100, 0, 6)
hCosTheta_doubleA = TH1F("hCosTheta_doubleA",";Cos(#theta)" , 200, -1, 1)
hPhi_Qweighted_doubleA = TH1F("hPhi_Qweighted_doubleA", "charge;#phi(rad);" ,100, 0, 6)
hCosTheta_Qweighted_doubleA = TH1F("hCosTheta_Qweighted_doubleA", "charge;Cos(#theta);" , 200, -1, 1)
hPhi_Tweighted_doubleA = TH1F("hPhi_Tweighted_doubleA", "hitTime;#phi(rad);" ,100, 0, 6)
hCosTheta_Tweighted_doubleA = TH1F("hCosTheta_Tweighted_doubleA", "hitTime;Cos(#theta);" , 200, -1, 1)

# average case
hAvePhiCosTheta_doubleA = TH2F("hAvePhiCosTheta_doubleA",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Qweighted_doubleA = TH2F("hAvePhiCosTheta_Qweighted_doubleA","charge;#phi(rad);" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tweighted_doubleA = TH2F("hAvePhiCosTheta_Tweighted_doubleA","hitTime;Cos(#theta);" ,100, 0, 6, 200, -1, 1)

hAvePhiCosTheta_Qfraction_doubleA = TH2F("hAvePhiCosTheta_Qfraction_double",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tfraction_doubleA = TH2F("hAvePhiCosTheta_Tfraction_double",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

hAvePhi_doubleA = TH1F("hAvePhi_doubleA",";#phi(rad);" ,100, 0, 6)
hAveCosTheta_doubleA = TH1F("hAveCosTheta_doubleA",";Cos(#theta);" , 200, -1, 1)
hAvePhi_Qweighted_doubleA = TH1F("hAvePhi_Qweighted_doubleA","charge;#phi(rad);" ,100, 0, 6)
hAvePhi_Tweighted_doubleA = TH1F("hAvePhi_Tweighted_doubleA","hitTime;#phi(rad);" ,100, 0, 6)
hAveCosTheta_Qweighted_doubleA = TH1F("hAveCosTheta_Qweighted_doubleA","charge;Cos(#theta);" , 200, -1, 1)
hAveCosTheta_Tweighted_doubleA = TH1F("hAveCosTheta_Tweighted_doubleA","hitTime;Cos(#theta);" , 200, -1, 1)
# ==============================
# 2D phi,theta; phi, theta Q weight; phi, theta T weight
hPhiCosTheta_doubleB = TH2F("hPhiCosTheta_doubleB",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Qweighted_doubleB = TH2F("hPhiCosTheta_Qweighted_doubleB", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hPhiCosTheta_Tweighted_doubleB = TH2F("hPhiCosTheta_Tweighted_doubleB", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

# 1D phi ; theta; phi; theta Q weight; phi; theta T weight
hPhi_doubleB = TH1F("hPhi_doubleB",";#phi(rad);" ,100, 0, 6)
hCosTheta_doubleB = TH1F("hCosTheta_doubleB",";Cos(#theta)" , 200, -1, 1)
hPhi_Qweighted_doubleB = TH1F("hPhi_Qweighted_doubleB", "charge;#phi(rad);" ,100, 0, 6)
hCosTheta_Qweighted_doubleB = TH1F("hCosTheta_Qweighted_doubleB", "charge;Cos(#theta);" , 200, -1, 1)
hPhi_Tweighted_doubleB = TH1F("hPhi_Tweighted_doubleB", "hitTime;#phi(rad);" ,100, 0, 6)
hCosTheta_Tweighted_doubleB = TH1F("hCosTheta_Tweighted_doubleB", "hitTime;Cos(#theta);" , 200, -1, 1)

# average case
hAvePhiCosTheta_doubleB = TH2F("hAvePhiCosTheta_doubleB",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Qweighted_doubleB = TH2F("hAvePhiCosTheta_Qweighted_doubleB","charge;#phi(rad);" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tweighted_doubleB = TH2F("hAvePhiCosTheta_Tweighted_doubleB","hitTime;Cos(#theta);" ,100, 0, 6, 200, -1, 1)

hAvePhiCosTheta_Qfraction_doubleB = TH2F("hAvePhiCosTheta_Qfraction_double",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)
hAvePhiCosTheta_Tfraction_doubleB = TH2F("hAvePhiCosTheta_Tfraction_double",";#phi(rad);Cos(#theta)" ,100, 0, 6, 200, -1, 1)

hAvePhi_doubleB = TH1F("hAvePhi_doubleB",";#phi(rad);" ,100, 0, 6)
hAveCosTheta_doubleB = TH1F("hAveCosTheta_doubleB",";Cos(#theta);" , 200, -1, 1)
hAvePhi_Qweighted_doubleB = TH1F("hAvePhi_Qweighted_doubleB","charge;#phi(rad);" ,100, 0, 6)
hAvePhi_Tweighted_doubleB = TH1F("hAvePhi_Tweighted_doubleB","hitTime;#phi(rad);" ,100, 0, 6)
hAveCosTheta_Qweighted_doubleB = TH1F("hAveCosTheta_Qweighted_doubleB","charge;Cos(#theta);" , 200, -1, 1)
hAveCosTheta_Tweighted_doubleB = TH1F("hAveCosTheta_Tweighted_doubleB","hitTime;Cos(#theta);" , 200, -1, 1)

#### save PMT info, average calculations into single-cluster, double-cluster, triple-cluster trees
treeI = TTree("TI","single_cluster")
evtID1 = array('l',[0])
sumQ1 = array('f',[0]) # sum charge of cluster
aveQ1 = array('f',[0]) # average charge of cluster
maxQ1 = array('f',[0]) # max charge of cluster
aveTime1 = array('f',[0]) # average PMT hitTime 
modeTime1 = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
qwTime1 = array('f',[0]) # charge-weighted time
earlyTime1 = array('f',[0]) # min(t)
lateTime1 = array('f',[0]) # max(t)
time_maxQ1 = array('f',[0]) # pmt time of the maxQ pmt 
evtx1 = array('f',[0])
evty1 = array('f',[0])
evtz1 = array('f',[0])
avePhi1 = array('f',[0]) # average of PMT phi-position
aveCos1 = array('f',[0]) # average of PMT cosTheta-position
aveQPhi1 = array('f',[0]) # Q-average of PMT phi-position
aveQCos1 = array('f',[0]) # Q-average of PMT cosTheta-position
aveTPhi1 = array('f',[0]) # T-average of PMT phi-position
aveTCos1 = array('f',[0]) # T-average of PMT cosTheta-position

runID1 = array('l',[0])
subrunID1 = array('l',[0])

treeI.Branch("runID",runID1,"runID/L")
treeI.Branch("subrunID",subrunID1,"subrunID/L")
treeI.Branch("eventID",evtID1,"eventID/L")
treeI.Branch("sumQ",sumQ1,"sumQ/F") 
treeI.Branch("aveQ",aveQ1,"aveQ/F") #ULong64_
treeI.Branch("maxQ", maxQ1,"maxQ/F")
treeI.Branch("aveTime", aveTime1,"aveTime/F")
treeI.Branch("modeTime", modeTime1,"modeTime/F")
treeI.Branch("qwTime", qwTime1, "qwTime/F")
treeI.Branch("earlyTime", earlyTime1,"earlyTime/F")
treeI.Branch("lateTime", lateTime1,"lateTime/F")
treeI.Branch("time_maxQ", time_maxQ1, "time_maxQ/F")
treeI.Branch("avePhi", avePhi1,"avePhi/F")
treeI.Branch("aveCos", aveCos1,"aveCos/F")
treeI.Branch("aveQPhi", aveQPhi1,"aveQPhi/F")
treeI.Branch("aveQCos", aveQCos1,"aveQCos/F")
treeI.Branch("aveTPhi", aveTPhi1,"aveTPhi/F")
treeI.Branch("aveTCos", aveTCos1,"aveTCos/F")

treeIIa = TTree("TIIa","double_cluster") # first cluster of the double-cluster, marked as 2a
evtID2a = array('l',[0]) # eventID2 is same
sumQ2a = array('f',[0])
aveQ2a = array('f',[0]) # average charge of cluster
maxQ2a = array('f',[0]) # max charge of cluster
aveTime2a = array('f',[0]) # average PMT hitTime 
modeTime2a = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
qwTime2a = array('f',[0])
earlyTime2a = array('f',[0]) # min(t)
lateTime2a = array('f',[0]) # max(t)
time_maxQ2a = array('f',[0])
evtx2a = array('f',[0])
evty2a = array('f',[0])
evtz2a = array('f',[0])
avePhi2a = array('f',[0]) # average of PMT phi-position
aveCos2a = array('f',[0]) # average of PMT cosTheta-position
aveQPhi2a = array('f',[0]) # Q-average of PMT phi-position
aveQCos2a = array('f',[0]) # Q-average of PMT cosTheta-position
aveTPhi2a = array('f',[0]) # T-average of PMT phi-position
aveTCos2a = array('f',[0]) # T-average of PMT cosTheta-position

runID2a = array('l',[0])
subrunID2a = array('l',[0])

treeIIa.Branch("runID",runID2a,"runID/L")
treeIIa.Branch("subrunID",subrunID2a,"subrunID/L")
treeIIa.Branch("eventID",evtID2a,"eventID/L")
treeIIa.Branch("sumQ",sumQ2a,"sumQ/F") #ULong64_
treeIIa.Branch("aveQ",aveQ2a,"aveQ/F") #ULong64_
treeIIa.Branch("maxQ", maxQ2a,"maxQ/F")
treeIIa.Branch("aveTime", aveTime2a,"aveTime/F")
treeIIa.Branch("modeTime", modeTime2a,"modeTime/F")
treeIIa.Branch("qwTime", qwTime2a, "qwTime/F")
treeIIa.Branch("earlyTime", earlyTime2a,"earlyTime/F")
treeIIa.Branch("lateTime", lateTime2a,"lateTime/F")
treeIIa.Branch("time_maxQ", time_maxQ2a, "time_maxQ/F")
treeIIa.Branch("avePhi", avePhi2a,"avePhi/F")
treeIIa.Branch("aveCos", aveCos2a,"aveCos/F")
treeIIa.Branch("aveQPhi", aveQPhi2a,"aveQPhi/F")
treeIIa.Branch("aveQCos", aveQCos2a,"aveQCos/F")
treeIIa.Branch("aveTPhi", aveTPhi2a,"aveTPhi/F")
treeIIa.Branch("aveTCos", aveTCos2a,"aveTCos/F")

treeIIb = TTree("TIIb","double_cluster") # first cluster of the double-cluster, marked as 2b
evtID2b = array('l',[0]) # eventID2 is same
sumQ2b = array('f',[0])
aveQ2b = array('f',[0]) # average charge of cluster
maxQ2b = array('f',[0]) # max charge of cluster
aveTime2b = array('f',[0]) # average PMT hitTime 
modeTime2b = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
qwTime2b = array('f',[0])
earlyTime2b = array('f',[0]) # min(t)
lateTime2b = array('f',[0]) # max(t)
time_maxQ2b = array('f',[0])
evtx2b = array('f',[0])
evty2b = array('f',[0])
evtz2b = array('f',[0])
avePhi2b = array('f',[0]) # average of PMT phi-position
aveCos2b = array('f',[0]) # average of PMT cosTheta-position
aveQPhi2b = array('f',[0]) # Q-average of PMT phi-position
aveQCos2b = array('f',[0]) # Q-average of PMT cosTheta-position
aveTPhi2b = array('f',[0]) # T-average of PMT phi-position
aveTCos2b = array('f',[0]) # T-average of PMT cosTheta-position
runID2b = array('l',[0])
subrunID2b = array('l',[0])

treeIIb.Branch("runID",runID2b,"runID/L")
treeIIb.Branch("subrunID",subrunID2b,"subrunID/L")
treeIIb.Branch("eventID",evtID2b,"eventID/L")
#treeIIb.Branch("evtx2b",evtx2b,"evt2b/F")
#treeIIb.Branch("evty2b",evty2b,"evt2b/F")
#treeIIb.Branch("evtz2b",evtz2b,"evt2b/F")
treeIIb.Branch("sumQ",sumQ2b,"sumQ/F") #ULong64_
treeIIb.Branch("aveQ",aveQ2b,"aveQ/F") #ULong64_
treeIIb.Branch("maxQ", maxQ2b,"maxQ/F")
treeIIb.Branch("aveTime", aveTime2b,"aveTime/F")
treeIIb.Branch("modeTime", modeTime2b,"modeTime/F")
treeIIb.Branch("qwTime", qwTime2b, "qwTime/F")
treeIIb.Branch("earlyTime", earlyTime2b,"earlyTime/F")
treeIIb.Branch("lateTime", lateTime2b,"lateTime/F")
treeIIb.Branch("time_maxQ", time_maxQ2b, "time_maxQ/F")
treeIIb.Branch("avePhi", avePhi2b,"avePhi/F")
treeIIb.Branch("aveCos", aveCos2b,"aveCos/F")
treeIIb.Branch("aveQPhi", aveQPhi2b,"aveQPhi/F")
treeIIb.Branch("aveQCos", aveQCos2b,"aveQCos/F")
treeIIb.Branch("aveTPhi", aveTPhi2b,"aveTPhi/F")
treeIIb.Branch("aveTCos", aveTCos2b,"aveTCos/F")

for fname in fileList:
   print "processing ", fname
   f = TFile(fname,"read")
   tree1 = f.Get("T1")
   tree2 = f.Get("T2")
   tree3 = f.Get("T3")
  
   tree4 = f.Get("T4")
   tree5 = f.Get("T5")
   tree6 = f.Get("T6")

   runID = array('l',[0])
   subrunID = array('l',[0])
   evtID = array('l',[0])
   evtx = array('f',[0])
   evty = array('f',[0])
   evtz = array('f',[0])
   qpeVal = array('f',[0]) #unsigned double 
   fpromptVal = array('f',[0])

   Nmax = 200

   #### pmt by pmt, single-cluster
   nPMTs1 = array('i',[0])
   pmtPhi1 = array('f',Nmax*[0])
   pmtCosTheta1 = array('f',Nmax*[0])
   pmttime1 = array('f',Nmax*[0])
   charge1 = array('f',Nmax*[0])
   
   #### pmt by pmt, double-cluster
   nPMTs2 = array('i',[0])
   pmtPhi2 = array('f',Nmax*[0])
   pmtCosTheta2 = array('f',Nmax*[0])
   pmttime2 = array('f',Nmax*[0])
   charge2 = array('f',Nmax*[0])
   
   #### pmt by pmt, triple-cluster
   nPMTs3 = array('i',[0])
   pmtPhi3 = array('f',Nmax*[0])
   pmtCosTheta3 = array('f',Nmax*[0])
   pmttime3 = array('f',Nmax*[0])
   charge3 = array('f',Nmax*[0])

   #### pmt by pmt, cluster 1
   nPMTs4 = array('i',[0])
   pmtPhi4 = array('f',Nmax*[0])
   pmtCosTheta4 = array('f',Nmax*[0])
   pmttime4 = array('f',Nmax*[0])
   charge4 = array('f',Nmax*[0])

   #### pmt by pmt, cluster 2
   nPMTs5 = array('i',[0])
   pmtPhi5 = array('f',Nmax*[0])
   pmtCosTheta5 = array('f',Nmax*[0])
   pmttime5 = array('f',Nmax*[0])
   charge5 = array('f',Nmax*[0])

   #### pmt by pmt, cluster 3
   nPMTs6 = array('i',[0])
   pmtPhi6 = array('f',Nmax*[0])
   pmtCosTheta6 = array('f',Nmax*[0])
   pmttime6 = array('f',Nmax*[0])
   charge6 = array('f',Nmax*[0])

   tree1.SetBranchAddress("runID",runID) #ULong64_
   tree1.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree1.SetBranchAddress("qpe", qpeVal)
   tree1.SetBranchAddress("eventID",evtID) #ULong64_
   tree1.SetBranchAddress("qpe", qpeVal)
   tree1.SetBranchAddress("fprompt", fpromptVal)
   tree1.SetBranchAddress("evtx", evtx)
   tree1.SetBranchAddress("evty", evty)
   tree1.SetBranchAddress("evtz", evtz)
   tree1.SetBranchAddress("nPMTs",nPMTs1)
   tree1.SetBranchAddress("pmtPhi", pmtPhi1)
   tree1.SetBranchAddress("pmtCosTheta", pmtCosTheta1)
   tree1.SetBranchAddress("pmttime", pmttime1)
   tree1.SetBranchAddress("charge", charge1)

   tree2.SetBranchAddress("runID",runID) #ULong64_
   tree2.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree2.SetBranchAddress("eventID",evtID) #ULong64_
   tree2.SetBranchAddress("qpe", qpeVal)
   tree2.SetBranchAddress("fprompt", fpromptVal)
   tree2.SetBranchAddress("evtx", evtx)
   tree2.SetBranchAddress("evty", evty)
   tree2.SetBranchAddress("evtz", evtz)
   tree2.SetBranchAddress("nPMTs",nPMTs2)
   tree2.SetBranchAddress("pmtPhi", pmtPhi2)
   tree2.SetBranchAddress("pmtCosTheta", pmtCosTheta2)
   tree2.SetBranchAddress("pmttime", pmttime2)
   tree2.SetBranchAddress("charge", charge2)

   tree3.SetBranchAddress("runID",runID) #ULong64_
   tree3.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree3.SetBranchAddress("eventID",evtID) #ULong64_
   tree3.SetBranchAddress("qpe", qpeVal)
   tree3.SetBranchAddress("fprompt", fpromptVal)
   tree3.SetBranchAddress("evtx", evtx)
   tree3.SetBranchAddress("evty", evty)
   tree3.SetBranchAddress("evtz", evtz)
   tree3.SetBranchAddress("nPMTs",nPMTs3)
   tree3.SetBranchAddress("pmtPhi", pmtPhi3)
   tree3.SetBranchAddress("pmtCosTheta", pmtCosTheta3)
   tree3.SetBranchAddress("pmttime", pmttime3)
   tree3.SetBranchAddress("charge", charge3)

   tree4.SetBranchAddress("runID",runID) #ULong64_
   tree4.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree4.SetBranchAddress("eventID",evtID) #ULong64_
   tree4.SetBranchAddress("qpe", qpeVal)
   tree4.SetBranchAddress("fprompt", fpromptVal)
   tree4.SetBranchAddress("evtx", evtx)
   tree4.SetBranchAddress("evty", evty)
   tree4.SetBranchAddress("evtz", evtz)
   tree4.SetBranchAddress("nPMTs",nPMTs4)
   tree4.SetBranchAddress("pmtPhi", pmtPhi4)
   tree4.SetBranchAddress("pmtCosTheta", pmtCosTheta4)
   tree4.SetBranchAddress("pmttime", pmttime4)
   tree4.SetBranchAddress("charge", charge4)

   tree5.SetBranchAddress("runID",runID) #ULong64_
   tree5.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree5.SetBranchAddress("eventID",evtID) #ULong64_
   tree5.SetBranchAddress("qpe", qpeVal)
   tree5.SetBranchAddress("fprompt", fpromptVal)
   tree5.SetBranchAddress("evtx", evtx)
   tree5.SetBranchAddress("evty", evty)
   tree5.SetBranchAddress("evtz", evtz)
   tree5.SetBranchAddress("nPMTs",nPMTs5)
   tree5.SetBranchAddress("pmtPhi", pmtPhi5)
   tree5.SetBranchAddress("pmtCosTheta", pmtCosTheta5)
   tree5.SetBranchAddress("pmttime", pmttime5)
   tree5.SetBranchAddress("charge", charge5)

   tree6.SetBranchAddress("runID",runID) #ULong64_
   tree6.SetBranchAddress("subrunID",subrunID) #ULong64_
   tree6.SetBranchAddress("eventID",evtID) #ULong64_
   tree6.SetBranchAddress("qpe", qpeVal)
   tree6.SetBranchAddress("fprompt", fpromptVal)
   tree6.SetBranchAddress("evtx", evtx)
   tree6.SetBranchAddress("evty", evty)
   tree6.SetBranchAddress("evtz", evtz)
   tree6.SetBranchAddress("nPMTs",nPMTs6)
   tree6.SetBranchAddress("pmtPhi", pmtPhi6)
   tree6.SetBranchAddress("pmtCosTheta", pmtCosTheta6)
   tree6.SetBranchAddress("pmttime", pmttime6)
   tree6.SetBranchAddress("charge", charge6)

   N1 = tree1.GetEntries()
   N2 = tree2.GetEntries()
   N3 = tree3.GetEntries()
  
   eventID_singleCluster = []
   eventID_doubleCluster = []
   eventID_tripleCluster = []

   for i in xrange(N1):
       tree1.GetEntry(i)
       #print "eventID ", evtID[0], nPMTs1[0]
       eventID_singleCluster.append(evtID[0])

   for i in xrange(N2):
       tree2.GetEntry(i)
       #print "eventID ", evtID[0], nPMTs2[0]
       eventID_doubleCluster.append(evtID[0])

   for i in xrange(N3):
       tree3.GetEntry(i)
       #print "eventID ", evtID[0], "pmts ",nPMTs3[0]
       eventID_tripleCluster.append(evtID[0])

   N4 = tree4.GetEntries()
   N5 = tree5.GetEntries()
   N6 = tree6.GetEntries()

   ### CLUSTER 1
   ### search for cluster1, get single-cluster or the first of the double-cluster
   for i in xrange(N4):
       tree4.GetEntry(i)
       ### remove 0 values
       ### check for single-cluster event!!
       eventID = evtID[0]
       if eventID in eventID_singleCluster: # CLUSTER 1
           print "process single cluster event", eventID, "npmt = ", nPMTs4[0]
           runID1[0] = runID[0]
           subrunID1[0] = subrunID[0]
           charge = []
           time = []
           sumTime = 0
           sumCharge = 0
           sumPhi = 0
           sumCosTheta = 0
           sumPhi_weightQ = 0
           sumCosTheta_weightQ = 0
           sumPhi_weightT = 0
           sumCosTheta_weightT = 0
           sumTime_weightQ = 0

           for kk in range(nPMTs4[0]):
                 #### calculate total hit-time and charge
                 charge.append(charge4[kk])
                 time.append(pmttime4[kk])

                 sumPhi      = sumPhi + pmtPhi4[kk]
                 sumCosTheta = sumCosTheta + pmtCosTheta4[kk]

                 sumPhi_weightQ = sumPhi_weightQ + pmtPhi4[kk]*charge4[kk]
                 sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta4[kk]*charge4[kk]

                 sumPhi_weightT = sumPhi_weightT + pmtPhi4[kk]*pmttime4[kk]
                 sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta4[kk]*pmttime4[kk]

                 hPhiCosTheta_singleCluster.Fill(pmtPhi4[kk], pmtCosTheta4[kk])
                 hPhiCosTheta_Qweighted_singleCluster.Fill(pmtPhi4[kk], pmtCosTheta4[kk], charge4[kk])
                 hPhiCosTheta_Tweighted_singleCluster.Fill(pmtPhi4[kk], pmtCosTheta4[kk], pmttime4[kk])

                 sumTime_weightQ = sumTime_weightQ + charge4[kk]*pmttime4[kk]

                 # 1D
                 hPhi_singleCluster.Fill(pmtPhi4[kk]); hCosTheta_singleCluster.Fill(pmtCosTheta4[kk])
                 hPhi_Qweighted_singleCluster.Fill(pmtPhi4[kk],charge4[kk]); hCosTheta_Qweighted_singleCluster.Fill(pmtCosTheta4[kk],charge4[kk])
                 hPhi_Tweighted_singleCluster.Fill(pmtPhi4[kk], pmttime4[kk]); hCosTheta_Tweighted_singleCluster.Fill(pmtCosTheta4[kk], pmttime4[kk])

           sumCharge = float(sum(charge))
           sumTime =   float(sum(time))
           hAvePhiCosTheta_singleCluster.Fill(sumPhi/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]))
           hAvePhiCosTheta_Qweighted_singleCluster.Fill(sumPhi/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]), sumCharge)
           hAvePhiCosTheta_Tweighted_singleCluster.Fill(sumCosTheta/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]), sumTime)

           hAvePhiCosTheta_Qfraction_singleCluster.Fill(sumPhi_weightQ/sumCharge,sumCosTheta_weightQ/sumCharge)
           hAvePhiCosTheta_Tfraction_singleCluster.Fill(sumPhi_weightT/sumTime,sumCosTheta_weightT/sumTime)

           hAvePhi_singleCluster.Fill(sumPhi/float(nPMTs4[0]))
           hAveCosTheta_singleCluster.Fill(sumCosTheta/float(nPMTs4[0]))
           hAvePhi_Qweighted_singleCluster.Fill(sumPhi/float(nPMTs4[0]),sum(charge))
           hAveCosTheta_Qweighted_singleCluster.Fill(sumCosTheta/float(nPMTs4[0]),sum(charge))
           hAvePhi_Tweighted_singleCluster.Fill(sumPhi/float(nPMTs4[0]),sum(time))
           hAveCosTheta_Tweighted_singleCluster.Fill(sumCosTheta/float(nPMTs4[0]),sum(time))

           sumQ1[0] = sumCharge
           evtID1[0] = eventID
           aveQ1[0] = sum(charge)/float(len(charge)) 
           maxQ1[0] = max(charge) 
           aveTime1[0] = sum(time)/float(len(time)) 
           modeTime1[0] = max(set(time), key = time.count)
           qwTime1[0] = sumTime_weightQ/sum(charge)
           earlyTime1[0] = min(time) 
           lateTime1[0] = max(time) 
           avePhi1[0] = sumPhi/float(nPMTs4[0]) 
           aveCos1[0] = sumCosTheta/float(nPMTs4[0])
           aveQPhi1[0] = sumPhi_weightQ/sumCharge
           aveQCos1[0] = sumCosTheta_weightQ/sumCharge
           aveTPhi1[0] = sumPhi_weightT/sumTime
           aveTCos1[0] = sumCosTheta_weightT/sumTime
           time_maxQ1[0] = time[charge.index(maxQ1[0])]
           treeI.Fill()

       ### check for double-cluster event in the first cluster tree!!
       if eventID  in eventID_doubleCluster: ## CLUSTER 1  == checkID:
           print "process double-cluster event", eventID, "the first cluster; still in tree1, npmt = ", nPMTs4[0]
           eventID = evtID[0]
           runID2a[0] = runID[0]
           subrunID2a[0] = subrunID[0]
           charge = []
           time = []
           sumTime = 0
           sumCharge = 0
           sumPhi = 0
           sumCosTheta = 0
           sumPhi_weightQ = 0
           sumCosTheta_weightQ = 0
           sumPhi_weightT = 0
           sumCosTheta_weightT = 0
           sumTime_weightQ = 0

           ### loop the first cluster, it could be single cluster, or the first cluster of the double-c, triple-c
           for kk in range(nPMTs4[0]):
                 #### calculate total hit-time and charge
                 charge.append(charge4[kk])
                 time.append(pmttime4[kk])
                 sumPhi      = sumPhi + pmtPhi4[kk]
                 sumCosTheta = sumCosTheta + pmtCosTheta4[kk]

                 sumPhi_weightQ = sumPhi_weightQ + pmtPhi4[kk]*charge4[kk]
                 sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta4[kk]*charge4[kk]

                 sumPhi_weightT = sumPhi_weightT + pmtPhi4[kk]*pmttime4[kk]
                 sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta4[kk]*pmttime4[kk]

                 hPhiCosTheta_doubleA.Fill(pmtPhi4[kk], pmtCosTheta4[kk])
                 hPhiCosTheta_Qweighted_doubleA.Fill(pmtPhi4[kk], pmtCosTheta4[kk], charge4[kk])
                 hPhiCosTheta_Tweighted_doubleA.Fill(pmtPhi4[kk], pmtCosTheta4[kk], pmttime4[kk])
                
                 sumTime_weightQ = sumTime_weightQ + charge4[kk]*pmttime4[kk]

                 # 1D
                 hPhi_doubleA.Fill(pmtPhi4[kk]); hCosTheta_singleCluster.Fill(pmtCosTheta4[kk])
                 hPhi_Qweighted_doubleA.Fill(pmtPhi4[kk],charge4[kk]); hCosTheta_Qweighted_singleCluster.Fill(pmtCosTheta4[kk],charge4[kk])
                 hPhi_Tweighted_doubleA.Fill(pmtPhi4[kk], pmttime4[kk]); hCosTheta_Tweighted_singleCluster.Fill(pmtCosTheta4[kk], pmttime4[kk])

           #print charge
           #print time
           sumCharge = float(sum(charge))
           sumTime =   float(sum(time))
           hAvePhiCosTheta_doubleA.Fill(sumPhi/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]))
           hAvePhiCosTheta_Qweighted_doubleA.Fill(sumPhi/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]), sumCharge)
           hAvePhiCosTheta_Tweighted_doubleA.Fill(sumCosTheta/float(nPMTs4[0]),sumCosTheta/float(nPMTs4[0]), sumTime)

           hAvePhiCosTheta_Qfraction_doubleA.Fill(sumPhi_weightQ/sumCharge,sumCosTheta_weightQ/sumCharge)
           hAvePhiCosTheta_Tfraction_doubleA.Fill(sumPhi_weightT/sumTime,sumCosTheta_weightT/sumTime)

           hAvePhi_doubleA.Fill(sumPhi/float(nPMTs4[0]))
           hAveCosTheta_doubleA.Fill(sumCosTheta/float(nPMTs4[0]))
           hAvePhi_Qweighted_doubleA.Fill(sumPhi/float(nPMTs4[0]),sum(charge))
           hAveCosTheta_Qweighted_doubleA.Fill(sumCosTheta/float(nPMTs4[0]),sum(charge))
           hAvePhi_Tweighted_doubleA.Fill(sumPhi/float(nPMTs4[0]),sum(time))
           hAveCosTheta_Tweighted_doubleA.Fill(sumCosTheta/float(nPMTs4[0]),sum(time))

           sumQ2a[0] = sumCharge
           evtID2a[0] = eventID
           aveQ2a[0] = sum(charge)/float(len(charge)) 
           maxQ2a[0] = max(charge) 
           aveTime2a[0] = sum(time)/float(len(time)) 
           modeTime2a[0] = max(set(time), key = time.count)
           qwTime2a[0] = sumTime_weightQ/sum(charge)
           earlyTime2a[0] = min(time) 
           lateTime2a[0] = max(time) 
           avePhi2a[0] = sumPhi/float(nPMTs4[0]) 
           aveCos2a[0] = sumCosTheta/float(nPMTs4[0])
           aveQPhi2a[0] = sumPhi_weightQ/sumCharge
           aveQCos2a[0] = sumCosTheta_weightQ/sumCharge
           aveTPhi2a[0] = sumPhi_weightT/sumTime
           aveTCos2a[0] = sumCosTheta_weightT/sumTime
           time_maxQ2a[0] = time[charge.index(maxQ2a[0])]
           treeIIa.Fill()

   #### CLUSTER 2
   #### search for the 2nd cluster for the 2nd of double-cluster     
   for i in xrange(N5):
       tree5.GetEntry(i)
       ### remove 0 values
       eventID = evtID[0]
       runID2b[0] = runID[0]
       subrunID2b[0] = subrunID[0]
       ### check for double-cluster event in the second cluster tree!!
       if eventID  in eventID_doubleCluster: ## CLUSTER 2  == checkID:
           print "process double-cluster event", eventID, "the second cluster; it's in tree2; so cluster 1->2, npmt = ", nPMTs5[0]
           charge = []
           time = []
           sumTime = 0
           sumCharge = 0
           sumPhi = 0
           sumCosTheta = 0
           sumPhi_weightQ = 0
           sumCosTheta_weightQ = 0
           sumPhi_weightT = 0
           sumCosTheta_weightT = 0
           sumTime_weightQ = 0
           ### loop the first cluster, it could be single cluster, or the first cluster of the double-c, triple-c
           for kk in range(nPMTs5[0]):
                 #### calculate total hit-time and charge
                 charge.append(charge5[kk])
                 time.append(pmttime5[kk])
                 sumPhi      = sumPhi + pmtPhi5[kk]
                 sumCosTheta = sumCosTheta + pmtCosTheta5[kk]
   
                 sumPhi_weightQ = sumPhi_weightQ + pmtPhi5[kk]*charge5[kk]
                 sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta5[kk]*charge5[kk]
   
                 sumPhi_weightT = sumPhi_weightT + pmtPhi5[kk]*pmttime5[kk]
                 sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta5[kk]*pmttime5[kk]
   
                 hPhiCosTheta_doubleB.Fill(pmtPhi5[kk], pmtCosTheta5[kk])
                 hPhiCosTheta_Qweighted_doubleB.Fill(pmtPhi5[kk], pmtCosTheta5[kk], charge5[kk])
                 hPhiCosTheta_Tweighted_doubleB.Fill(pmtPhi5[kk], pmtCosTheta5[kk], pmttime5[kk])
  
                 sumTime_weightQ = sumTime_weightQ + charge5[kk]*pmttime5[kk]

                 # 1D
                 hPhi_doubleB.Fill(pmtPhi5[kk]); hCosTheta_singleCluster.Fill(pmtCosTheta5[kk])
                 hPhi_Qweighted_doubleB.Fill(pmtPhi5[kk],charge5[kk]); hCosTheta_Qweighted_singleCluster.Fill(pmtCosTheta5[kk],charge5[kk])
                 hPhi_Tweighted_doubleB.Fill(pmtPhi5[kk], pmttime5[kk]); hCosTheta_Tweighted_singleCluster.Fill(pmtCosTheta5[kk], pmttime5[kk])
  
           #print charge
           #print time
           sumCharge = float(sum(charge))
           sumTime =   float(sum(time))
           hAvePhiCosTheta_doubleB.Fill(sumPhi/float(nPMTs5[0]),sumCosTheta/float(nPMTs5[0]))
           hAvePhiCosTheta_Qweighted_doubleB.Fill(sumPhi/float(nPMTs5[0]),sumCosTheta/float(nPMTs5[0]), sumCharge)
           hAvePhiCosTheta_Tweighted_doubleB.Fill(sumCosTheta/float(nPMTs5[0]),sumCosTheta/float(nPMTs5[0]), sumTime)
   
           hAvePhiCosTheta_Qfraction_doubleB.Fill(sumPhi_weightQ/sumCharge,sumCosTheta_weightQ/sumCharge)
           hAvePhiCosTheta_Tfraction_doubleB.Fill(sumPhi_weightT/sumTime,sumCosTheta_weightT/sumTime)
   
           hAvePhi_doubleB.Fill(sumPhi/float(nPMTs5[0]))
           hAveCosTheta_doubleB.Fill(sumCosTheta/float(nPMTs5[0]))
           hAvePhi_Qweighted_doubleB.Fill(sumPhi/float(nPMTs5[0]),sum(charge))
           hAveCosTheta_Qweighted_doubleB.Fill(sumCosTheta/float(nPMTs5[0]),sum(charge))
           hAvePhi_Tweighted_doubleB.Fill(sumPhi/float(nPMTs5[0]),sum(time))
           hAveCosTheta_Tweighted_doubleB.Fill(sumCosTheta/float(nPMTs5[0]),sum(time))
  
           sumQ2b[0] = sumCharge
           evtID2b[0] = eventID
           aveQ2b[0] = sum(charge)/float(len(charge)) 
           maxQ2b[0] = max(charge) 
           aveTime2b[0] = sum(time)/float(len(time)) 
           modeTime2b[0] = max(set(time), key = time.count)
           qwTime2b[0] = sumTime_weightQ/sum(charge)
           earlyTime2b[0] = min(time) 
           lateTime2b[0] = max(time) 
           avePhi2b[0] = sumPhi/float(nPMTs5[0]) 
           aveCos2b[0] = sumCosTheta/float(nPMTs5[0])
           aveQPhi2b[0] = sumPhi_weightQ/sumCharge
           aveQCos2b[0] = sumCosTheta_weightQ/sumCharge
           aveTPhi2b[0] = sumPhi_weightT/sumTime
           aveTCos2b[0] = sumCosTheta_weightT/sumTime
           time_maxQ2b[0] = time[charge.index(maxQ2b[0])]
           treeIIb.Fill()
   
   #### CLUSTER 3 
   #### search for the 2nd cluster for the 2nd of double-cluster     
   for i in xrange(N6):
       tree6.GetEntry(i)
       ### remove 0 values
       eventID = evtID[0]
       runID2b[0] = runID[0]
       subrunID2b[0] = subrunID[0]
       ### check for double-cluster event in the second cluster tree!!
       if eventID  in eventID_doubleCluster: ## CLUSTER 3   == checkID:
           print "process double-cluster event", eventID, "the second cluster; it's in tree3; so cluster 1->3, npmt = ", nPMTs6[0]
           charge = []
           time = []
           sumTime = 0
           sumCharge = 0
           sumPhi = 0
           sumCosTheta = 0
           sumPhi_weightQ = 0
           sumCosTheta_weightQ = 0
           sumPhi_weightT = 0
           sumCosTheta_weightT = 0
           sumTime_weightQ = 0
           ### loop the first cluster, it could be single cluster, or the first cluster of the double-c, triple-c
           for kk in range(nPMTs6[0]):
                 #### calculate total hit-time and charge
                 charge.append(charge6[kk])
                 time.append(pmttime6[kk])
                 sumPhi      = sumPhi + pmtPhi6[kk]
                 sumCosTheta = sumCosTheta + pmtCosTheta6[kk]
  
                 sumPhi_weightQ = sumPhi_weightQ + pmtPhi6[kk]*charge6[kk]
                 sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta6[kk]*charge6[kk]
  
                 sumPhi_weightT = sumPhi_weightT + pmtPhi6[kk]*pmttime6[kk]
                 sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta6[kk]*pmttime6[kk]
  
                 hPhiCosTheta_doubleB.Fill(pmtPhi6[kk], pmtCosTheta6[kk])
                 hPhiCosTheta_Qweighted_doubleB.Fill(pmtPhi6[kk], pmtCosTheta6[kk], charge6[kk])
                 hPhiCosTheta_Tweighted_doubleB.Fill(pmtPhi6[kk], pmtCosTheta6[kk], pmttime6[kk])
 
                 sumTime_weightQ = sumTime_weightQ + charge6[kk]*pmttime6[kk]

                 # 1D
                 hPhi_doubleB.Fill(pmtPhi6[kk]); hCosTheta_singleCluster.Fill(pmtCosTheta6[kk])
                 hPhi_Qweighted_doubleB.Fill(pmtPhi6[kk],charge6[kk]); hCosTheta_Qweighted_singleCluster.Fill(pmtCosTheta6[kk],charge6[kk])
                 hPhi_Tweighted_doubleB.Fill(pmtPhi6[kk], pmttime6[kk]); hCosTheta_Tweighted_singleCluster.Fill(pmtCosTheta6[kk], pmttime6[kk])
 
           #print charge
           #print time
           sumCharge = float(sum(charge))
           sumTime =   float(sum(time))
           hAvePhiCosTheta_doubleB.Fill(sumPhi/float(nPMTs6[0]),sumCosTheta/float(nPMTs6[0]))
           hAvePhiCosTheta_Qweighted_doubleB.Fill(sumPhi/float(nPMTs6[0]),sumCosTheta/float(nPMTs6[0]), sumCharge)
           hAvePhiCosTheta_Tweighted_doubleB.Fill(sumCosTheta/float(nPMTs6[0]),sumCosTheta/float(nPMTs6[0]), sumTime)
  
           hAvePhiCosTheta_Qfraction_doubleB.Fill(sumPhi_weightQ/sumCharge,sumCosTheta_weightQ/sumCharge)
           hAvePhiCosTheta_Tfraction_doubleB.Fill(sumPhi_weightT/sumTime,sumCosTheta_weightT/sumTime)
  
           hAvePhi_doubleB.Fill(sumPhi/float(nPMTs6[0]))
           hAveCosTheta_doubleB.Fill(sumCosTheta/float(nPMTs6[0]))
           hAvePhi_Qweighted_doubleB.Fill(sumPhi/float(nPMTs6[0]),sum(charge))
           hAveCosTheta_Qweighted_doubleB.Fill(sumCosTheta/float(nPMTs6[0]),sum(charge))
           hAvePhi_Tweighted_doubleB.Fill(sumPhi/float(nPMTs6[0]),sum(time))
           hAveCosTheta_Tweighted_doubleB.Fill(sumCosTheta/float(nPMTs6[0]),sum(time))
  
           evtID2b[0] = eventID
           sumQ2b[0] = sumCharge
           aveQ2b[0] = sumCharge/float(len(charge)) 
           maxQ2b[0] = max(charge) 
           aveTime2b[0] = sum(time)/float(len(time)) 
           modeTime2b[0] = max(set(time), key = time.count)
           qwTime2b[0] = sumTime_weightQ/sum(charge)
           earlyTime2b[0] = min(time) 
           lateTime2b[0] = max(time) 
           avePhi2b[0] = sumPhi/float(nPMTs6[0]) 
           aveCos2b[0] = sumCosTheta/float(nPMTs6[0])
           aveQPhi2b[0] = sumPhi_weightQ/sumCharge
           aveQCos2b[0] = sumCosTheta_weightQ/sumCharge
           aveTPhi2b[0] = sumPhi_weightT/sumTime
           aveTCos2b[0] = sumCosTheta_weightT/sumTime
           time_maxQ2b[0] = time[charge.index(maxQ2b[0])]
           treeIIb.Fill()
          
   fnew.cd()
   hAvePhiCosTheta_singleCluster.Write("hAvePhiCosTheta_single_%s"%fname[12:17])
   hAvePhiCosTheta_Qfraction_singleCluster.Write("hAvePhiCosTheta_Qfraction_single_%s"%fname[12:17])
   hAvePhiCosTheta_Tfraction_singleCluster.Write("hAvePhiCosTheta_Tfraction_single_%s"%fname[12:17])
 
   hAvePhiCosTheta_doubleA.Write("hAvePhiCosTheta_doubleA_%s"%fname[12:17])
   hAvePhiCosTheta_Qfraction_doubleA.Write("hAvePhiCosTheta_Qfraction_doubleA_%s"%fname[12:17])
   hAvePhiCosTheta_Tfraction_doubleA.Write("hAvePhiCosTheta_Tfraction_doubleA_%s"%fname[12:17])
 
   hAvePhiCosTheta_doubleB.Write("hAvePhiCosTheta_doubleB_%s"%fname[12:17])
   hAvePhiCosTheta_Qfraction_doubleB.Write("hAvePhiCosTheta_Qfraction_doubleB_%s"%fname[12:17])
   hAvePhiCosTheta_Tfraction_doubleB.Write("hAvePhiCosTheta_Tfraction_doubleB_%s"%fname[12:17])

   fnew.cd()
   treeI.Write()
   treeIIa.Write()
   treeIIb.Write()

   #fnew.Close()
