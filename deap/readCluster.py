import ROOT
import sys, os, getopt
from ROOT import *
from rat import *
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array

f = TFile("My_deap_cal_030708_0000.root")

tree1 = f.Get("T1")
tree2 = f.Get("T2")
tree3 = f.Get("T3")

evtID = array('l',[0])
evtx = array('f',[0])
evty = array('f',[0])
evtz = array('f',[0])
qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])

Nmax = 200

#### pmt by pmt, cluster 1
nPMTs1 = array('i',[0])
pmtPhi1 = array('f',Nmax*[0])
pmtCosTheta1 = array('f',Nmax*[0])
pmttime1 = array('f',Nmax*[0])
charge1 = array('f',Nmax*[0])

#### pmt by pmt, cluster 2
nPMTs2 = array('i',[0])
pmtPhi2 = array('f',Nmax*[0])
pmtCosTheta2 = array('f',Nmax*[0])
pmttime2 = array('f',Nmax*[0])
charge2 = array('f',Nmax*[0])

#### pmt by pmt, cluster 3
nPMTs3 = array('i',[0])
pmtPhi3 = array('f',Nmax*[0])
pmtCosTheta3 = array('f',Nmax*[0])
pmttime3 = array('f',Nmax*[0])
charge3 = array('f',Nmax*[0])

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

hAvePhiCosTheta_cluster1 = TH2F("hAvePhiCosTheta_cluster1",";#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Qweighted_cluster1 = TH2F("hAvePhiCosTheta_Qweighted_cluster1", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Tweighted_cluster1 = TH2F("hAvePhiCosTheta_Tweighted_cluster1", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)

hAvePhiCosTheta_cluster2 = TH2F("hAvePhiCosTheta_cluster2",";#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Qweighted_cluster2 = TH2F("hAvePhiCosTheta_Qweighted_cluster2", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Tweighted_cluster2 = TH2F("hAvePhiCosTheta_Tweighted_cluster2", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)

hAvePhiCosTheta_cluster3 = TH2F("hAvePhiCosTheta_cluster3",";#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Qweighted_cluster3 = TH2F("hAvePhiCosTheta_Qweighted_cluster3", "charge;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)
hAvePhiCosTheta_Tweighted_cluster3 = TH2F("hAvePhiCosTheta_Tweighted_cluster3", "hitTime;#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1)

N1 = tree1.GetEntries()
N2 = tree2.GetEntries()
N3 = tree3.GetEntries()


print "process cluster 1"
for i in xrange(N1):
    tree1.GetEntry(i)
    ### remove 0 values
    sumTime = 0
    sumCharge = 0
    sumPhi = 0
    sumCosTheta = 0
    sumPhi_weightQ = 0
    sumCosTheta_weightQ = 0
    sumPhi_weightT = 0
    sumCosTheta_weightT = 0

    print "eventID ", evtID[0], nPMTs1[0] 
    for kk in range(nPMTs1[0]):
        #### calculate total hit-time and charge
        sumTime     = sumTime + pmttime1[kk]
        sumCharge   = sumCharge + charge1[kk]
        sumPhi      = sumPhi + pmtPhi1[kk]
        sumCosTheta = sumCosTheta + pmtCosTheta1[kk]
          
        sumPhi_weightQ = sumPhi_weightQ + pmtPhi1[kk]*charge1[kk]
        sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta1[kk]*charge1[kk]
    
        sumPhi_weightT = sumPhi_weightT + pmtPhi1[kk]*pmttime1[kk]
        sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta1[kk]*pmttime1[kk]

    hAvePhiCosTheta_cluster1.Fill(sumPhi/float(nPMTs1[0]),sumCosTheta/float(nPMTs1[0]))
    hAvePhiCosTheta_Qweighted_cluster1.Fill(sumPhi_weightQ/(sumCharge*float(nPMTs1[0])),sumCosTheta/(sumCharge*float(nPMTs1[0])))
    hAvePhiCosTheta_Tweighted_cluster1.Fill(sumPhi_weightT/(sumTime*float(nPMTs1[0])),sumCosTheta/(sumTime*float(nPMTs1[0])))

print "process cluster 2"
for i in xrange(N2):
    tree2.GetEntry(i)
    ### remove 0 values
    sumTime = 0
    sumCharge = 0
    sumPhi = 0
    sumCosTheta = 0
    sumPhi_weightQ = 0
    sumCosTheta_weightQ = 0
    sumPhi_weightT = 0
    sumCosTheta_weightT = 0

    print "eventID ", evtID[0], nPMTs2[0]
    for kk in range(nPMTs2[0]):
        #### calculate total hit-time and charge
        sumTime     = sumTime + pmttime2[kk]
        sumCharge   = sumCharge + charge2[kk]
        sumPhi      = sumPhi + pmtPhi2[kk]
        sumCosTheta = sumCosTheta + pmtCosTheta2[kk]

        sumPhi_weightQ = sumPhi_weightQ + pmtPhi2[kk]*charge2[kk]
        sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta2[kk]*charge2[kk]

        sumPhi_weightT = sumPhi_weightT + pmtPhi2[kk]*pmttime2[kk]
        sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta2[kk]*pmttime2[kk]

    hAvePhiCosTheta_cluster2.Fill(sumPhi/float(nPMTs2[0]),sumCosTheta/float(nPMTs2[0]))
    hAvePhiCosTheta_Qweighted_cluster2.Fill(sumPhi_weightQ/(sumCharge*float(nPMTs2[0])),sumCosTheta/(sumCharge*float(nPMTs2[0])))
    hAvePhiCosTheta_Tweighted_cluster2.Fill(sumPhi_weightT/(sumTime*float(nPMTs2[0])),sumCosTheta/(sumTime*float(nPMTs2[0])))

print "process cluster 3"
for i in xrange(N3):
    tree3.GetEntry(i)
    ### remove 0 values
    sumTime = 0
    sumCharge = 0
    sumPhi = 0
    sumCosTheta = 0
    sumPhi_weightQ = 0
    sumCosTheta_weightQ = 0
    sumPhi_weightT = 0
    sumCosTheta_weightT = 0

    print "eventID ", evtID[0], "pmts ",nPMTs3[0]
    for kk in range(nPMTs3[0]):
        #### calculate total hit-time and charge
        sumTime     = sumTime + pmttime3[kk]
        sumCharge   = sumCharge + charge3[kk]
        sumPhi      = sumPhi + pmtPhi3[kk]
        sumCosTheta = sumCosTheta + pmtCosTheta3[kk]

        sumPhi_weightQ = sumPhi_weightQ + pmtPhi3[kk]*charge3[kk]
        sumCosTheta_weightQ = sumCosTheta_weightQ + pmtCosTheta3[kk]*charge3[kk]

        sumPhi_weightT = sumPhi_weightT + pmtPhi3[kk]*pmttime3[kk]
        sumCosTheta_weightT = sumCosTheta_weightT + pmtCosTheta3[kk]*pmttime3[kk]

    hAvePhiCosTheta_cluster3.Fill(sumPhi/float(nPMTs3[0]),sumCosTheta/float(nPMTs3[0]))
    hAvePhiCosTheta_Qweighted_cluster3.Fill(sumPhi_weightQ/(sumCharge*float(nPMTs3[0])),sumCosTheta/(sumCharge*float(nPMTs3[0])))
    hAvePhiCosTheta_Tweighted_cluster3.Fill(sumPhi_weightT/(sumTime*float(nPMTs3[0])),sumCosTheta/(sumTime*float(nPMTs3[0])))


fnew = TFile("saveClusters.root","recreate")
fnew.cd()
hAvePhiCosTheta_cluster1.Write()
hAvePhiCosTheta_Qweighted_cluster1.Write()
hAvePhiCosTheta_Tweighted_cluster1.Write()

hAvePhiCosTheta_cluster2.Write()
hAvePhiCosTheta_Qweighted_cluster2.Write()
hAvePhiCosTheta_Tweighted_cluster2.Write()

hAvePhiCosTheta_cluster3.Write()
hAvePhiCosTheta_Qweighted_cluster3.Write()
hAvePhiCosTheta_Tweighted_cluster3.Write()




