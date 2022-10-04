import ROOT
import sys, os, getopt
from ROOT import *
from rat import *
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
from numpy import sqrt
import csv
from array import array
import operator
### mostly care about 2a to 2b for the double-cluster

r = 850.0 # AV radius
c_light = 299.792*0.05  # alpha particle is 5% of speed of light
fList = ["saveOptMC_Merged_OptHalfQ_data.root"]

#saveOptMC_OptHalfQ_MCusher_Po210_AV_5000evts_thickness0p00001.root"]

#saveOptMC_OptHalfQ_MC_Run30681_Po210_vacuum_avBulk_thin0p00001_5000evts.root"]

fname = fList[0]

fnew = TFile("results_"+fname,"recreate")

### for single cluster events
H2D_singleC_ReconHitPoint1_ave = TH2D("H2D_singleC_ReconHitPoint1_ave","hitpoint1 of singleCluster, use averaged cos,phi",100, 0, 6, 100, -1, 1);
H2D_singleC_ReconHitPoint1_maxQ = TH2D("H2D_singleC_ReconHitPoint1_maxQ","hitpoint1 of singleCluster, use maxQ cos,phi",100, 0, 6, 100, -1, 1);
H2D_singleC_ReconHitPoint1_aveQ = TH2D("H2D_singleC_ReconHitPoint1_aveQ","hitpoint1 of singleCluster, use Q-weighted cos,phi",100, 0, 6, 100, -1, 1);
### for double cluster events
H2D_doubleC_ReconHitPoint1_ave = TH2D("H2D_doubleC_ReconHitPoint1_ave","hitpoint1 of doubleCluster, use averaged cos,phi",100, 0, 6, 100, -1, 1);
H2D_doubleC_ReconHitPoint1_maxQ = TH2D("H2D_doubleC_ReconHitPoint1_maxQ","hitpoint1 of doubleCluster, use maxQ cos,phi",100, 0, 6, 100, -1, 1);
H2D_doubleC_ReconHitPoint1_aveQ = TH2D("H2D_doubleC_ReconHitPoint1_aveQ","hitpoint1 of doubleCluster, use Q-weighted cos,phi",100, 0, 6, 100, -1, 1);

HToFdiff = TH1F("HToFdiff","time difference of the two clusters;ns;",100, 0, 10)
HTaveDiff = TH1F("HTaveDiff","time difference of the two clusters;ns;",80, -200, 200)
HTmodeDiff = TH1F("HTmodeDiff","mode time difference of the two cluster;ns;",80, -200, 200)
HTqwDiff = TH1F("HTqwDiff","qw time difference of the two cluster;ns;",80, -200, 200)
HT21 = TH1F("HT21","2b late - 2a early between the two cluster;ns;",80, -200, 200)
HT21_1 = TH1F("HT21_1","2b early - 2a late between the two cluster;ns;",80, -200, 200)

HTearlyDiff = TH1F("HTearlyDiff","early time difference of the two clusters;ns;",80, -200, 200)
### dist vs early time
H2D_distVsT1_early = TH2D("H2D_distVsT1_early",";ns;mm", 300, 0, 3000, 100,0,1700)
H2D_distVsT2_early = TH2D("H2D_distVsT2_early",";ns;mm", 300, 0, 3000, 100,0,1700)
H2D_distVsTdiff_early = TH2D("H2D_distVsTdiff_early",";ns;mm", 300, 0, 3000, 100,0,1700)

H2D_distVsT1_mode  = TH2D("H2D_distVsT1_mode",";ns;mm",300, 0, 3000, 100,0,1700)
H2D_distVsT2_mode  = TH2D("H2D_distVsT2_mode",";ns;mm",300, 0, 3000, 100,0,1700)
H2D_distVsTdiff_mode  = TH2D("H2D_distVsTdiff_mode",";ns;mm", 80, -200, 200, 100, 0, 1700)

H2D_distVsT1_qw  = TH2D("H2D_distVsT1_qw",";ns;mm", 300, 0, 3000, 100, 0, 1700)
H2D_distVsT2_qw  = TH2D("H2D_distVsT2_qw",";ns;mm", 300, 0, 3000, 100, 0, 1700)
H2D_distVsTdiff_qw  = TH2D("H2D_distVsTdiff_qw",";ns;mm", 100, 0, 200, 100, 0, 1700)

H2D_distQVsTdiff_ave  = TH2D("H2D_distQVsTdiff_ave",";ns;mm", 80, -200, 200, 100, 0, 1700)
H2D_distQVsTdiff_mode  = TH2D("H2D_distQVsTdiff_mode",";ns;mm", 80, -200, 200, 100, 0, 1700)
H2D_distQVsTdiff_qw  = TH2D("H2D_distQVsTdiff_qw",";ns;mm", 80, -200, 200, 100, 0, 1700)

H2D_distMaxQvsTdiff = TH2D("H2D_distMaxQvsTdiff",";ns;mm", 80, -200, 200, 100, 0, 1700)

eventID1 = []
Q_ave1 = []
phi_ave1 = []
ct_ave1 = []
phi_maxQ1 = []
ct_maxQ1 = []
phi_aveQ1 = []
ct_aveQ1 = []
pmtt_ave1 = []
pmtt_mode1 = []
pmtt_qw1 = []
pmtt_early1 = []
pmtt_late1 = []
s_maxQ1 = []
s_aveQ1 = []

eventID2a = []
Q_sum2a = []
Q_ave2a = []
phi_ave2a = []
ct_ave2a = []
phi_maxQ2a = []
ct_maxQ2a = []
phi_aveQ2a = []
ct_aveQ2a = []
pmtt_ave2a = []
pmtt_mode2a = []
pmtt_qw2a = []
pmtt_early2a = []
pmtt_late2a = []
s_maxQ2a = []
s_aveQ2a = []
s_maxQPhi2a = []
s_maxQCos2a = []
s_tmaxq2a = []

eventID2b = []
Q_sum2b = []
Q_ave2b = []
phi_ave2b = []
ct_ave2b = []
phi_maxQ2b = []
ct_maxQ2b = []
phi_aveQ2b = []
ct_aveQ2b = []
pmtt_ave2b = []
pmtt_mode2b = []
pmtt_qw2b = []
pmtt_early2b = []
pmtt_late2b = []
s_maxQ2b = []
s_aveQ2b = []
s_maxQPhi2b = []
s_maxQCos2b = []
s_tmaxq2b = []

pmttDiff_ave = []
pmttDiff_early = []
pmttDiff_late = []

pmtt_maxQ_2a = []
pmtt_maxQ_2b = []

val_runID1 = []
val_subrunID1 = []

val_runID2a = []
val_subrunID2a = []

val_runID2b = []
val_subrunID2b = []

for fname in fList:
  f = TFile(fname,"read")
  treeI = f.Get("TI")
  evtID1 = array('l',[0])
  aveQ1 = array('f',[0]) # average charge of cluster
  maxQ1 = array('f',[0]) # max charge of cluster
  aveTime1 = array('f',[0]) # average PMT hitTime 
  modeTime1 = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
  qwTime1 = array('f',[0])
  earlyTime1 = array('f',[0]) # min(t)
  lateTime1 = array('f',[0]) # max(t)
  avePhi1 = array('f',[0]) # average of PMT phi-position
  aveCos1 = array('f',[0]) # average of PMT cosTheta-position
  maxQPhi1 = array('f',[0])
  maxQCos1 = array('f',[0])
  aveQPhi1 = array('f',[0]) # Q-average of PMT phi-position
  aveQCos1 = array('f',[0]) # Q-average of PMT cosTheta-position
  aveTPhi1 = array('f',[0]) # T-average of PMT phi-position
  aveTCos1 = array('f',[0]) # T-average of PMT cosTheta-position
  time_maxQ1 = array('f',[0])
  runID1 = array('f',[0]) 
  subrunID1 = array('f',[0])

  treeI.SetBranchAddress("runID",runID1)
  treeI.SetBranchAddress("subrunID",subrunID1)
  treeI.SetBranchAddress("eventID",evtID1)
  treeI.SetBranchAddress("aveQ",aveQ1)
  treeI.SetBranchAddress("maxQ", maxQ1)
  treeI.SetBranchAddress("aveTime", aveTime1)
  treeI.SetBranchAddress("modeTime", modeTime1)
  treeI.SetBranchAddress("qwTime", qwTime1)
  treeI.SetBranchAddress("earlyTime", earlyTime1)
  treeI.SetBranchAddress("lateTime", lateTime1)
  treeI.SetBranchAddress("time_maxQ", time_maxQ1)
  treeI.SetBranchAddress("avePhi", avePhi1)
  treeI.SetBranchAddress("aveCos", aveCos1)
  treeI.SetBranchAddress("maxQPhi", maxQPhi1)
  treeI.SetBranchAddress("maxQCos", maxQCos1)
  treeI.SetBranchAddress("aveQPhi", aveQPhi1)
  treeI.SetBranchAddress("aveQCos", aveQCos1)
  treeI.SetBranchAddress("aveTPhi", aveTPhi1)
  treeI.SetBranchAddress("aveTCos", aveTCos1)

  treeIIa = f.Get("TIIa")
  evtID2a = array('l',[0])
  sumQ2a = array('f',[0])
  aveQ2a = array('f',[0]) # average charge of cluster
  maxQ2a = array('f',[0]) # max charge of cluster
  aveTime2a = array('f',[0]) # average PMT hitTime 
  modeTime2a = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
  qwTime2a = array('f',[0])
  earlyTime2a = array('f',[0]) # min(t)
  lateTime2a = array('f',[0]) # max(t)
  avePhi2a = array('f',[0]) # average of PMT phi-position
  aveCos2a = array('f',[0]) # average of PMT cosTheta-position
  maxQPhi2a = array('f',[0])
  maxQCos2a = array('f',[0])
  aveQPhi2a = array('f',[0]) # Q-average of PMT phi-position
  aveQCos2a = array('f',[0]) # Q-average of PMT cosTheta-position
  aveTPhi2a = array('f',[0]) # T-average of PMT phi-position
  aveTCos2a = array('f',[0]) # T-average of PMT cosTheta-position
  tmaxQ2a = array('f',[0])
  runID2a = array('f',[0])
  subrunID2a = array('f',[0])

  treeIIa.SetBranchAddress("runID",runID2a)
  treeIIa.SetBranchAddress("subrunID",subrunID2a)
  treeIIa.SetBranchAddress("eventID",evtID2a)
  treeIIa.SetBranchAddress("aveQ",aveQ2a)
  treeIIa.SetBranchAddress("sumQ",sumQ2a)
  treeIIa.SetBranchAddress("maxQ", maxQ2a)
  treeIIa.SetBranchAddress("aveTime", aveTime2a)
  treeIIa.SetBranchAddress("modeTime", modeTime2a)
  treeIIa.SetBranchAddress("qwTime", qwTime2a)
  treeIIa.SetBranchAddress("earlyTime", earlyTime2a)
  treeIIa.SetBranchAddress("lateTime", lateTime2a)
  treeIIa.SetBranchAddress("avePhi", avePhi2a)
  treeIIa.SetBranchAddress("aveCos", aveCos2a)
  treeIIa.SetBranchAddress("maxQPhi",maxQPhi2a)
  treeIIa.SetBranchAddress("maxQCos",maxQCos2a)
  treeIIa.SetBranchAddress("aveQPhi", aveQPhi2a)
  treeIIa.SetBranchAddress("aveQCos", aveQCos2a)
  treeIIa.SetBranchAddress("aveTPhi", aveTPhi2a)
  treeIIa.SetBranchAddress("aveTCos", aveTCos2a)
  treeIIa.SetBranchAddress("aveTCos", aveTCos2a)
  treeIIa.SetBranchAddress("tmaxQ", tmaxQ2a)
  treeIIa.SetBranchAddress("maxQPhi", maxQPhi2a)
  treeIIa.SetBranchAddress("maxQCos", maxQCos2a)

  treeIIb = f.Get("TIIb")
  qave2b = array('l',[0])
  evtID2b = array('l',[0])
  sumQ2b = array('f',[0])
  aveQ2b = array('f',[0]) # average charge of cluster
  maxQ2b = array('f',[0]) # max charge of cluster
  aveTime2b = array('f',[0]) # average PMT hitTime 
  modeTime2b = array('f',[0]) # PMT mode time, max(set(t), key = t.count)
  qwTime2b = array('f',[0])
  earlyTime2b = array('f',[0]) # min(t)
  lateTime2b = array('f',[0]) # max(t)
  avePhi2b = array('f',[0]) # average of PMT phi-position
  aveCos2b = array('f',[0]) # average of PMT cosTheta-position
  maxQPhi2b = array('f',[0])
  maxQCos2b = array('f',[0])
  aveQPhi2b = array('f',[0]) # Q-average of PMT phi-position
  aveQCos2b = array('f',[0]) # Q-average of PMT cosTheta-position
  aveTPhi2b = array('f',[0]) # T-average of PMT phi-position
  aveTCos2b = array('f',[0]) # T-average of PMT cosTheta-position
  runID2b = array('f',[0])
  subrunID2b = array('f',[0])
  tmaxQ2b = array('f',[0])

  treeIIb.SetBranchAddress("runID",runID2b)
  treeIIb.SetBranchAddress("subrunID",subrunID2b)
  treeIIb.SetBranchAddress("eventID",evtID2b)
  treeIIb.SetBranchAddress("aveQ",aveQ2b)
  treeIIb.SetBranchAddress("sumQ",sumQ2b)
  treeIIb.SetBranchAddress("maxQ", maxQ2b)
  treeIIb.SetBranchAddress("aveTime", aveTime2b)
  treeIIb.SetBranchAddress("modeTime", modeTime2b)
  treeIIb.SetBranchAddress("qwTime", qwTime2b)
  treeIIb.SetBranchAddress("earlyTime", earlyTime2b)
  treeIIb.SetBranchAddress("lateTime", lateTime2b)
  treeIIb.SetBranchAddress("avePhi", avePhi2b)
  treeIIb.SetBranchAddress("aveCos", aveCos2b)
  treeIIb.SetBranchAddress("maxQPhi", maxQPhi2b)
  treeIIb.SetBranchAddress("maxQCos", maxQCos2b)
  treeIIb.SetBranchAddress("aveQPhi", aveQPhi2b)
  treeIIb.SetBranchAddress("aveQCos", aveQCos2b)
  treeIIb.SetBranchAddress("aveTPhi", aveTPhi2b)
  treeIIb.SetBranchAddress("aveTCos", aveTCos2b)
  treeIIb.SetBranchAddress("tmaxQ", tmaxQ2b)
  treeIIb.SetBranchAddress("maxQPhi", maxQPhi2b)
  treeIIb.SetBranchAddress("maxQCos", maxQCos2b)

  N1 = treeI.GetEntries()
  N2a = treeIIa.GetEntries()
  N2b = treeIIb.GetEntries()
  print N2a
  for i in range(N1): 
     treeI.GetEntry(i)
     val_runID1.append(runID1[0])
     val_subrunID1.append(subrunID1[0])
     Q_ave1.append(aveQ1[0])
     phi_ave1.append(avePhi1[0])
     ct_ave1.append(aveCos1[0])
     phi_aveQ1.append(aveQPhi1[0])
     ct_aveQ1.append(aveQCos1[0])
     pmtt_ave1.append(aveTime1[0])
     pmtt_early1.append(earlyTime1[0])
     pmtt_late1.append(lateTime1[0])
     pmtt_mode1.append(modeTime1[0])
     pmtt_qw1.append(qwTime1[0])
     s_maxQ1.append(maxQ1[0])
     eventID1.append(evtID1[0])

     H2D_singleC_ReconHitPoint1_ave.Fill(avePhi1[0], aveCos1[0])
     H2D_singleC_ReconHitPoint1_maxQ.Fill(maxQPhi1[0], maxQCos1[0])
     H2D_singleC_ReconHitPoint1_aveQ.Fill(aveQPhi1[0], aveQCos1[0])

  for i in range(N2a): 
     treeIIa.GetEntry(i)
     eventID2a.append(evtID2a[0])
     val_runID2a.append(runID2a[0])
     val_subrunID2a.append(subrunID2a[0])
     Q_sum2a.append(sumQ2a[0])
     Q_ave2a.append(aveQ2a[0])
     phi_ave2a.append(avePhi2a[0])
     ct_ave2a.append(aveCos2a[0])
     phi_aveQ2a.append(aveQPhi2a[0])
     ct_aveQ2a.append(aveQCos2a[0])
     pmtt_ave2a.append(aveTime2a[0])
     pmtt_early2a.append(earlyTime2a[0])
     pmtt_late2a.append(lateTime2a[0])
     pmtt_mode2a.append(modeTime2a[0])
     pmtt_qw2a.append(qwTime2a[0])
     s_maxQ2a.append(maxQ2a[0])

     s_maxQPhi2a.append(maxQPhi2a[0])
     s_maxQCos2a.append(maxQCos2a[0])
     s_tmaxq2a.append(tmaxQ2a[0])

     #HAveCosTheta_Phi_doubleA.Fill(avePhi2a[0],aveCos2a[0])
     #HAveQCosTheta_Phi_doubleA.Fill(aveQPhi2a[0],aveQCos2a[0])

  ### save PMT information
  ###     0      1       2      3   4     5    6      7       8           9      10     11        12      13        14       15
  ### pmt t_ave, t_mode, t_qw, (phi,ct), (phiQ,ctQ)   runID   subrunID    aveQ   sumQ   t_early   t_late  phi_maxQ  ct_maxQ  t_maxQ
  id_pmt2a = dict(zip(eventID2a, zip(pmtt_ave2a, pmtt_mode2a, pmtt_qw2a, phi_ave2a, ct_ave2a, phi_aveQ2a, ct_aveQ2a, val_runID2a, val_subrunID2a, Q_ave2a, Q_sum2a,  pmtt_early2a, pmtt_late2a, s_maxQPhi2a, s_maxQCos2a, s_tmaxq2a)))

  for i in range(N2b):
     treeIIb.GetEntry(i)
     eventID2b.append(evtID2b[0])
     val_runID2b.append(runID2b[0])
     val_subrunID2b.append(subrunID2b[0])
     Q_sum2b.append(sumQ2b[0])
     Q_ave2b.append(aveQ2b[0])
     phi_ave2b.append(avePhi2b[0])
     ct_ave2b.append(aveCos2b[0])
     phi_aveQ2b.append(aveQPhi2b[0])
     ct_aveQ2b.append(aveQCos2b[0])
     pmtt_ave2b.append(aveTime2b[0])
     pmtt_early2b.append(earlyTime2b[0])
     pmtt_late2b.append(lateTime2b[0])
     pmtt_mode2b.append(modeTime2b[0])
     pmtt_qw2b.append(qwTime2b[0])
     s_maxQ2b.append(maxQ2b[0])

     s_maxQPhi2b.append(maxQPhi2b[0])
     s_maxQCos2b.append(maxQCos2b[0])
     s_tmaxq2b.append(tmaxQ2a[0])

     #HAveCosTheta_Phi_doubleB.Fill(avePhi2b[0],aveCos2b[0])
     #HAveQCosTheta_Phi_doubleB.Fill(aveQPhi2b[0],aveQCos2b[0])

  id_pmt2b = dict(zip(eventID2b, zip(pmtt_ave2b, pmtt_mode2b, pmtt_qw2b, phi_ave2b, ct_ave2b, phi_aveQ2b, ct_aveQ2b, val_runID2b, val_subrunID2b, Q_ave2b, Q_sum2b,  pmtt_early2b, pmtt_late2b, s_maxQPhi2b, s_maxQCos2b, s_tmaxq2b)))

  sorted_pmt2a = dict(sorted(id_pmt2a.items(), key=operator.itemgetter(0)))
  sorted_pmt2b = dict(sorted(id_pmt2b.items(), key=operator.itemgetter(0)))

  #print "pmt2a", sorted_pmt2a
  #print "pmt2b", sorted_pmt2b
  #print id_pmt2a

  #print sorted_pmt2a

  #print id_pmt2a
  if len(eventID2a) != len(eventID2b):
      print "!!!! something wrong!"

  ###     0      1       2      3   4     5    6     7      8          9      10      11      12      13       14          15 
  ### pmt t_ave, t_mode, t_qw, (phi,ct), (phiQ,ctQ)  runID  subrunID   aveQ   sumQ    tearly  tlate   (phi_qmax cos_qmax)  t_maxq
  sort_aveTime2a = []
  sort_modeTime2a = []
  sort_earlyTime2a = []
  sort_lateTime2a = []
  sort_tqw2a = [] # q-weighted time
  sort_pos2a = []
  sort_posQ2a = []
  sort_eventID2a = []
  sort_runID2a = []
  sort_subrunID2a = []
  sort_qAve2a = []
  sort_qSum2a = []
  sort_pos_qmax2a = []
  sort_tmaxq_2a = []
  sort_avePhiCos_2a = []
  sort_maxQPhiCos_2a = []
  sort_aveQPhiCos_2a = []

  i = 0
  for key in sorted_pmt2a:
      #if key != sorted_pmt2b.keys()[i]:
      #   print 'things wrong!'
      #   continue
      sort_aveTime2a.append(id_pmt2a[key][0])
      sort_modeTime2a.append(id_pmt2a[key][1])
      sort_earlyTime2a.append(id_pmt2a[key][11])
      sort_lateTime2a.append(id_pmt2a[key][12])
      sort_tqw2a.append(id_pmt2a[key][2])
      sort_tmaxq_2a.append(id_pmt2a[key][15])
      sort_eventID2a.append(key)
      ## ave positions
      phi = id_pmt2a[key][3] 
      ct = id_pmt2a[key][4]
      sort_avePhiCos_2a.append((phi,ct))
      xx = TVector3(r*cos(phi)*sqrt(1-ct*ct),r*sin(phi)*sqrt(1-ct*ct),r*ct)
      sort_pos2a.append(xx)
      
      ## q-weighted positions
      phi1 = id_pmt2a[key][5]
      ct1 = id_pmt2a[key][6]
      sort_aveQPhiCos_2a.append((phi1, ct1))
      xx1 = TVector3(r*cos(phi1)*sqrt(1-ct1*ct1),r*sin(phi1)*sqrt(1-ct1*ct1),r*ct1)
      sort_posQ2a.append(xx1)
      
      ## maxQ positions
      phi_qmax = id_pmt2a[key][13]
      ct_qmax = id_pmt2a[key][14]
      sort_maxQPhiCos_2a.append((phi_qmax, ct_qmax))
      xx_qmax = TVector3(r*cos(phi_qmax)*sqrt(1-ct_qmax*ct_qmax),r*sin(phi_qmax)*sqrt(1-ct_qmax*ct_qmax),r*ct_qmax)
      sort_pos_qmax2a.append(xx_qmax)

      sort_runID2a.append(id_pmt2a[key][7])
      sort_subrunID2a.append(id_pmt2a[key][8])
      sort_qAve2a.append(id_pmt2a[key][9])
      sort_qSum2a.append(id_pmt2a[key][10])
      i = i+1

  ###     0      1       2      3   4     5    6     7      8          9      10     11       12      13       14        15
  ### pmt t_ave, t_mode, t_qw, (phi,ct), (phiQ,ctQ)  runID  subrunID   aveQ   sumQ   tearly   tlate   phi_qmax cos_qmax  t_maxq
  sort_aveTime2b = []
  sort_modeTime2b = []
  sort_earlyTime2b = []
  sort_lateTime2b = []
  sort_tqw2b = [] # q-weighted time
  sort_pos2b = []
  sort_posQ2b = []
  sort_eventID2b = []
  sort_runID2b = []
  sort_subrunID2b = []
  sort_qAve2b = []
  sort_qSum2b = []
  sort_avePhiCos_2b = []
  sort_maxQPhiCos_2b = []
  sort_aveQPhiCos_2b = []

  sort_pos_qmax2b = []
  sort_tmaxq_2b = []
  i = 0
  for key in sorted_pmt2b:
      if key != sorted_pmt2b.keys()[i]:
         print 'unmathced, things wrong!'
         continue
      sort_aveTime2b.append(id_pmt2b[key][0])
      sort_modeTime2b.append(id_pmt2b[key][1])
      sort_earlyTime2b.append(id_pmt2b[key][11])
      sort_lateTime2b.append(id_pmt2b[key][12])
      sort_tqw2b.append(id_pmt2b[key][2])

      sort_tmaxq_2b.append(id_pmt2b[key][15])

      sort_eventID2b.append(key)
      phi = id_pmt2b[key][3] 
      ct = id_pmt2b[key][4]
      sort_avePhiCos_2b.append((phi,ct))
      # print phi, ct
      xx = TVector3(r*cos(phi)*sqrt(1-ct*ct),r*sin(phi)*sqrt(1-ct*ct),r*ct)
      sort_pos2b.append(xx)

      phi1 = id_pmt2b[key][5]
      ct1 = id_pmt2b[key][6]
      sort_aveQPhiCos_2b.append((phi1,ct1))
      xx1 = TVector3(r*cos(phi1)*sqrt(1-ct1*ct1),r*sin(phi1)*sqrt(1-ct1*ct1),r*ct1)
      sort_posQ2b.append(xx1)

      phi_qmax = id_pmt2b[key][13]
      ct_qmax = id_pmt2b[key][14]
      sort_maxQPhiCos_2b.append((phi_qmax, ct_qmax))
      xx_qmax = TVector3(r*cos(phi_qmax)*sqrt(1-ct_qmax*ct_qmax),r*sin(phi_qmax)*sqrt(1-ct_qmax*ct_qmax),r*ct_qmax)
      sort_pos_qmax2b.append(xx_qmax)

      sort_runID2b.append(id_pmt2b[key][7])
      sort_subrunID2b.append(id_pmt2b[key][8])
      sort_qAve2b.append(id_pmt2b[key][9])
      sort_qSum2b.append(id_pmt2b[key][10])

      i = i +1

  ### match and draw plots
  count_ave = 0; count_mode = 0; count_qw = 0
  if len(sort_eventID2a)==len(sort_eventID2b):
      print "data length mismatch, everything is fine"

  size00 = 0
  if len(sort_eventID2b)>len(sort_eventID2a):
     size00 = len(sort_eventID2a)
  else:
     print "!!! length mis-match"
     size00 = len(sort_eventID2b)

  for i in range( size00 ):
      delta_t_mode = sort_modeTime2b[i] - sort_modeTime2a[i]
      delta_t_ave = sort_aveTime2b[i] - sort_aveTime2a[i]
      delta_t_qw = abs(sort_tqw2b[i] - sort_tqw2a[i])

      delta_t_early = sort_earlyTime2b[i] - sort_earlyTime2a[i]
      delta_t_21 = sort_lateTime2b[i] - sort_earlyTime2a[i] # late in 2b - early in 2a 
      delta_t_21_1 = sort_lateTime2a[i] - sort_earlyTime2b[i] # late in 2b - early in 2a 

      HTaveDiff.Fill(delta_t_ave)
      HTmodeDiff.Fill(delta_t_mode)
      HTqwDiff.Fill(delta_t_qw)
      HT21.Fill(delta_t_21)
      HT21_1.Fill(delta_t_21_1) 

      xa = sort_pos2a[i]
      xb = sort_pos2b[i]
      dist = (xb-xa).Mag()

      tof = dist/c_light
      HToFdiff.Fill(tof)
      
      xaQ = sort_posQ2a[i]
      xbQ = sort_posQ2b[i]
      distQ = (xbQ-xaQ).Mag()

      X_maxq2a = sort_pos_qmax2a[i]
      X_maxq2b = sort_pos_qmax2b[i]
      dist_maxq = (X_maxq2a - X_maxq2b).Mag()
      delta_t_maxq = sort_tmaxq_2a[i] - sort_tmaxq_2b[i]

      phi_ave = sort_avePhiCos_2a[i][0]
      ct_ave = sort_avePhiCos_2a[i][1]

      phi_Qave = sort_aveQPhiCos_2a[i][0]
      ct_Qave = sort_aveQPhiCos_2a[i][1]

      phi_qmax = sort_maxQPhiCos_2a[i][0] 
      ct_qmax = sort_maxQPhiCos_2a[i][1]

      if delta_t_maxq<0:
          phi_ave = sort_avePhiCos_2b[i][0]
          ct_ave = sort_avePhiCos_2b[i][1]

          phi_Qave = sort_aveQPhiCos_2b[i][0]
          ct_Qave = sort_aveQPhiCos_2b[i][1]

          phi_qmax = sort_maxQPhiCos_2b[i][0]
          ct_qmax = sort_maxQPhiCos_2b[i][1]

      H2D_doubleC_ReconHitPoint1_ave.Fill(phi_ave, ct_ave)
      H2D_doubleC_ReconHitPoint1_aveQ.Fill(phi_Qave, ct_Qave)
      H2D_doubleC_ReconHitPoint1_maxQ.Fill(phi_qmax, ct_qmax)

      #if(delta_t_ave<0):
      #    print sort_eventID2a[i], "tave<0" #sort_eventID2b[i]
      #    print "run",val_runID2a[i], val_subrunID2a[i]
      #    count_ave += 1

      #if(delta_t_mode>10 and delta_t_mode<100):
      print "delta t_mode, runID, runID, eventID, chargeAve2a, chargeAve2b, chargeSum"
      print "tmode= ", delta_t_mode, ">0, run", int(sort_runID2a[i]), "event", sort_eventID2a[i], sort_qAve2a[i], sort_qAve2b[i], sort_qSum2a[i]
      #H_checkQ1a.Fill(sort_qAve2a[i])
      #H_checkQ1b.Fill(sort_qAve2b[i])
      count_mode += 1
      H2D_distVsT1_qw.Fill(sort_tqw2a[i], dist)
      H2D_distVsT2_qw.Fill(sort_tqw2b[i], dist)
      H2D_distVsTdiff_qw.Fill(delta_t_qw, dist)

      if(delta_t_mode<0):

        ### save the 2nd cluster point as the 1st point
        print "tmode= ", delta_t_mode, "<0, run", int(sort_runID2a[i]), "event", sort_eventID2a[i], sort_qAve2a[i], sort_qAve2b[i], sort_qSum2b[i]
        #H_checkQ2a.Fill(sort_qAve2a[i])
        #H_checkQ2b.Fill(sort_qAve2b[i])
        count_mode += 1

      #if(delta_t_qw<0):
      #    print sort_eventID2a[i], "tqw<0"
      #    print val_runID2a[i], val_subrunID2a[i]
      #    count_qw += 1

      #Hdist.Fill(dist)
      #HdistQ.Fill(distQ)
      #Hdist_maxq.Fill(dist_maxq)

      HToFdiff.Fill(tof)
      #H2D_distVsT1_ave.Fill(sort_aveTime2a[i],dist)
      #H2D_distVsT2_ave.Fill(sort_aveTime2b[i],dist)
      #H2D_distVsTdiff_ave.Fill(delta_t_ave,dist)

      #H2D_distVsT1_early.Fill(sort_earlyTime2a[i],dist)
      #H2D_distVsT2_early.Fill(sort_earlyTime2b[i],dist)
      #H2D_distVsTdiff_early.Fill(delta_t_early,dist)

      #H2D_distVsT1_mode.Fill(sort_modeTime2a[i],dist)
      #H2D_distVsT2_mode.Fill(sort_modeTime2b[i],dist)
      H2D_distVsTdiff_mode.Fill(delta_t_mode,dist)

      H2D_distQVsTdiff_ave.Fill(delta_t_ave, distQ)
      H2D_distQVsTdiff_mode.Fill(delta_t_mode, distQ)
      H2D_distQVsTdiff_qw.Fill(delta_t_qw, distQ)

      #H2D_distQVsT21.Fill(delta_t_21, dist)
      #H2D_distVsT21.Fill(delta_t_21, distQ)

      H2D_distMaxQvsTdiff.Fill(delta_t_maxq, dist_maxq) 

  print "single-cluster", N1
  print "double-cluster", N2a, N2b

print "count deltaT<0", count_ave, count_mode, count_qw #, count_early
#Hdist.Draw()
H2D_distVsTdiff_qw.Draw('colz')
raw_input()

fnew.cd()
H2D_singleC_ReconHitPoint1_ave.Write()
H2D_singleC_ReconHitPoint1_maxQ.Write()
H2D_singleC_ReconHitPoint1_aveQ.Write()

H2D_doubleC_ReconHitPoint1_ave.Write()
H2D_doubleC_ReconHitPoint1_maxQ.Write()
H2D_doubleC_ReconHitPoint1_aveQ.Write()

H2D_distVsT1_mode.Write()
H2D_distVsT2_mode.Write()
H2D_distVsTdiff_mode.Write()

H2D_distVsT1_qw.Write()
H2D_distVsT2_qw.Write()
H2D_distVsTdiff_qw.Write()

H2D_distQVsTdiff_ave.Write()
H2D_distQVsTdiff_mode.Write()
H2D_distQVsTdiff_qw.Write()

H2D_distMaxQvsTdiff.Write()
