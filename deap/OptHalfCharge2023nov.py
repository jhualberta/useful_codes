#!/usr/bin/python2
# double-cluster analysis
# @author Jie Hu 2022 Sep
# updated Nov 2023
# to run: python OptHalfCharge2023nov.py /project/6004969/data/v5.14.0/cal/run030627/deap_cal_030627_0000.root
#### NOTE: Modified on Apr 2023, for data
# for an event, tag as single- or double-cluster event by checking the 6-highest charge PMT locations.
# do the tagging by using 4 types of pmt charges! pmt_qpe, pmt_qnscb, pmt_qprompt, pmt_qaap (removed alberta after-pulse charge).
# save the 6-highest charge PMT info (id, charge, time).
# for a double-cluster event, separate the saved PMTs into the 1st and 2nd clusters.
# we also fill hit time residual with different charges for the event.
# This comes with RAT, so source rat.env first!

import os
import couchdb
import sys
import glob
from rat import *
#from rat import PMTInfoUtil
from ROOT import *
import numpy as np
from numpy import *
import math
from datetime import datetime, time as datetime_time, timedelta
from array import array
import pickle
global qc

PI = np.pi
radius = 850
QmaxCut = 0 #5 #100 ## the qmax PMT must > 100 pC
QminCut = 0 #5 #20 ## the minimum q of selected PMT must > 20 pC
Nmax = 255 ## save how many PMTs for an event
Nsubpeak = 50 ## default npulses*nsubpeak, as 50

pmtInfo = RAT.PMTInfoUtil.GetPMTInfoUtil()

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
fileName = os.path.basename(file0) ## remove the directory path

fin = TFile(file0); ## with full directory path

fout = TFile("opt2023Nov_"+fileName,"RECREATE");
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])

Dir=sys.argv[1]
#AllFiles=glob.glob(Dir)
#AllFiles = sorted(AllFiles)

### save as different cluster events
tree1 = TTree("T1","all events")

### each tree saves event variables as well as PMT clusters
### the event variables are same for each tree; but each tree has different PMT clusters

### these are for each event, common for all cluster
runID = array('l',[0])
subrunID = array('l',[0])
evtID = array('l',[0])
mbx = array('f',[0])
mby = array('f',[0])
mbz = array('f',[0])
tf2x = array('f',[0])
tf2y = array('f',[0])
tf2z = array('f',[0])
eventTimeVal = array('f',[0])
tf2_ch_tVal = array('f',[0])
ncluster = array('i',[0])
qpeVal = array('f',[0]) #unsigned double 
fmaxpeVal = array('f',[0])
fpromptVal = array('f',[0])
nSCBayesVal = array('f',[0])
rprompt60BayesVal = array('f',[0])
nhitVal = array('i',[0])
#### these values are reserved for master cuts, ChSource STR cuts
numEarlyPulsesVal = array('f',[0])
subeventNval = array('f',[0])
eventTimeVal = array('f',[0])
pulseindexfirstgarVal = array('f',[0]) ## EV.GetPulseIndexFirstGAr()
cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qPE<0.1
### scintlikelihood calculation, qPE V1740 variables
scintlikeVal = array('f',[0])
scintlike16Val = array('f',[0])
scintlikeqVal = array('f',[0])
llneutronflashVal = array('f',[0])
neckVetoVal = array('f',[0])
#qPEnoSat_10000Val = array('f',[0])
#qPEnoSat_5000Val = array('f',[0])

#### pmt by pmt, cluster 1
nPMTs1 = array('i',[0])
listPMTid1 = array('I', Nmax*[0])
nPulse1 = array('i', [0])
nSubpeaks1 = array('i',[0])
### all pmts in an event
pmtPhi1 = array('f',Nmax*[0])
pmtCosTheta1 = array('f',Nmax*[0])
#pmttime1 = array('f',Nmax*[0]) # earliest SubpeakTime
#pmttime_uni1 = array('f',Nmax*[0]) #uniform 
#pmttime_unicor1 = array('f',Nmax*[0]) # corrected
subpeaktime1 = array('f',Nmax*Nsubpeak*[0])
pmtQpe1  = array('f', Nmax*[0]) # pmtqpe
pmtQnscb1 = array('f', Nmax*[0])
pmtQprompt1 = array('f', Nmax*[0])
pmtQaap1 = array('f', Nmax*[0]) 

## default = 0, if tagged by pmt.qpe=1; pmt.qnscb=2; pmt.qprompt=3; pmt.qaap=4;
tag_singleCluster = array('i',[0])
tag_doubleCluster = array('i',[0])

### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### searched 6PMTs by different charge methods 

NsortPMTs = 6 # select 6 highest charge PMTs
## --- qpe method --- (1+2) * (costheta, phi, Q)
## A1: 1st of double; B1: 2nd of double
pmtPhiQpe_single1 = array('f',NsortPMTs*[0])
pmtCosThetaQpe_single1 = array('f',NsortPMTs*[0])
pmtQpe_single1  = array('f', NsortPMTs*[0]) # pmtqpe
pmtPhiQpe_doubleA1 = array('f',NsortPMTs*[0])
pmtCosThetaQpe_doubleA1 = array('f',NsortPMTs*[0])
pmtQpe_doubleA1 = array('f',NsortPMTs*[0])
pmtPhiQpe_doubleB1 = array('f',NsortPMTs*[0])
pmtCosThetaQpe_doubleB1 = array('f',NsortPMTs*[0])
pmtQpe_doubleB1 = array('f',NsortPMTs*[0])
## --- qnscb method ---
pmtPhiQnscb_single1 = array('f',NsortPMTs*[0])
pmtCosThetaQnscb_single1 = array('f',NsortPMTs*[0])
pmtQnscb_single1  = array('f', NsortPMTs*[0]) # pmtqnscb
pmtPhiQnscb_doubleA1 = array('f',NsortPMTs*[0])
pmtCosThetaQnscb_doubleA1 = array('f',NsortPMTs*[0])
pmtQnscb_doubleA1 = array('f',NsortPMTs*[0])
pmtPhiQnscb_doubleB1 = array('f',NsortPMTs*[0])
pmtCosThetaQnscb_doubleB1 = array('f',NsortPMTs*[0])
pmtQnscb_doubleB1 = array('f',NsortPMTs*[0])
## --- qprompt method ---
pmtPhiQprompt_single1 = array('f',NsortPMTs*[0])
pmtCosThetaQprompt_single1 = array('f',NsortPMTs*[0])
pmtQprompt_single1  = array('f', NsortPMTs*[0]) # pmtqpei
pmtPhiQprompt_doubleA1 = array('f',NsortPMTs*[0])
pmtCosThetaQprompt_doubleA1 = array('f',NsortPMTs*[0])
pmtQprompt_doubleA1 = array('f',NsortPMTs*[0])
pmtPhiQprompt_doubleB1 = array('f',NsortPMTs*[0])
pmtCosThetaQprompt_doubleB1 = array('f',NsortPMTs*[0])
pmtQprompt_doubleB1 = array('f',NsortPMTs*[0])
## --- qaap method ---
pmtPhiQaap_single1 = array('f',NsortPMTs*[0])
pmtCosThetaQaap_single1 = array('f',NsortPMTs*[0])
pmtQaap_single1  = array('f', NsortPMTs*[0]) # pmtqpe
pmtPhiQaap_doubleA1 = array('f',NsortPMTs*[0])
pmtCosThetaQaap_doubleA1 = array('f',NsortPMTs*[0])
pmtQaap_doubleA1 = array('f',NsortPMTs*[0])
pmtPhiQaap_doubleB1 = array('f',NsortPMTs*[0])
pmtCosThetaQaap_doubleB1 = array('f',NsortPMTs*[0])
pmtQaap_doubleB1 = array('f',NsortPMTs*[0])

tree1.Branch("runID",runID,"runID/L")
tree1.Branch("subrunID",subrunID,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("numEarlyPulses", numEarlyPulsesVal,"numEarlyPulses/F") 
tree1.Branch("subeventN",  subeventNval, " subeventN/F") 
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("nSCBayes", nSCBayesVal, "nSCBayes/F") # GetNSCBayes()
tree1.Branch("rprompt60Bayes", rprompt60BayesVal, "rprompt60Bayes/F")
tree1.Branch("fmaxpe", fmaxpeVal, "fmaxpe/F")
tree1.Branch("mbx", mbx, "mbx/F")
tree1.Branch("mby", mby, "mby/F")
tree1.Branch("mbz", mbz, "mbz/F")
tree1.Branch("tf2x", tf2x, "tf2x/F")
tree1.Branch("tf2y", tf2y, "tf2y/F")
tree1.Branch("tf2z", tf2z, "tf2z/F")
tree1.Branch("eventTime",eventTimeVal, "eventTime/F")
tree1.Branch("ncluster",ncluster,"ncluster/I")
### new added variables: nhits, ..., qPEnoSat_5000
tree1.Branch("nhit", nhitVal, "nhit/I")
tree1.Branch("pulseindexfirstgar", pulseindexfirstgarVal, "pulseindexfirstgar/F")
tree1.Branch("tf2_ch_time", tf2_ch_tVal, "tf2_ch_time/F") ## timefit2 ch time
tree1.Branch("cft2r", cft2r,"cft2r/F")
tree1.Branch("cfb3r", cfb3r,"cfb3r/F")
tree1.Branch("neckVeto", neckVetoVal, "neckVeto/F")
tree1.Branch("scintlike", scintlikeVal, "scintlike/F")
tree1.Branch("scintlike16", scintlike16Val, "scintlike16/F")
tree1.Branch("scintlikeq", scintlikeqVal, "scintlikeq/F")
tree1.Branch("llneutronflash", llneutronflashVal, "llneutronflash/F")
tree1.Branch("nPMTs",nPMTs1, "nPMTs/I")
tree1.Branch("listPMTid",listPMTid1, "listPMTid[nPMTs]/I")
tree1.Branch("pmtPhi", pmtPhi1, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtCosTheta1, "pmtCosTheta[nPMTs]/F")
tree1.Branch("pmttime", pmttime1, "pmttime[nPMTs]/F")
tree1.Branch("nPulse", nPulse1, "nPulse/I")
tree1.Branch("nSubpeaks", nSubpeaks1, "nSubpeaks/I") # nPMT*nPulse*nSubpeak
tree1.Branch("subpeaktime", subpeaktime1, "subpeaktime[nSubpeaks]/F")
tree1.Branch("pmtQpe", pmtQpe1, "pmtQpe[nPMTs]/F")
tree1.Branch("pmtQnscb", pmtQnscb1, "pmtQnscb[nPMTs]/F")
tree1.Branch("pmtQprompt", pmtQprompt1, "pmtQprompt[nPMTs]/F")
tree1.Branch("pmtQaap", pmtQaap, "pmtQaap[nPMTs]/F")
## pmt.qpe: tag == 1; pmt.qnscb: tag == 2; pmt.qprompt: tag == 3; pmt.qaap: tag == 4
tree1.Branch("tag_singleCluster",tag_singleCluster,"tag_singleCluster/I")
tree1.Branch("tag_doubleCluster",tag_doubleCluster,"tag_doubleCluster/I")

## --- qpe method ---
tree1.Branch("pmtPhiQpe_single", pmtPhiQpe_single1,"pmtPhiQpe_single[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQpe_single", pmtCosThetaQpe_single1, "pmtCosThetaQpe_single[NsortPMTs]/F")
tree1.Branch("pmtQpe_single", pmtQpe_single1, "pmtQpe_single[NsortPMTs]/F")
tree1.Branch("pmtPhiQpe_doubleA", pmtPhiQpe_doubleA1, "pmtPhiQpe_doubleA[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQpe_doubleA", pmtCosThetaQpe_doubleA1,"pmtCosThetaQpe_doubleA[NsortPMTs]/F")
tree1.Branch("pmtQpe_doubleA", pmtQpe_doubleA1, "pmtQpe_doubleA[NsortPMTs]/F")
tree1.Branch("pmtPhiQpe_doubleB", pmtPhiQpe_doubleB1, "pmtPhiQpe_doubleB[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQpe_doubleB", pmtCosThetaQpe_doubleB1, "pmtCosThetaQpe_doubleB[NsortPMTs]/F")
tree1.Branch("pmtQpe_doubleB", pmtQpe_doubleB1, "pmtQpe_doubleB1[NsortPMTs]/F")
## --- qnscb method ---
tree1.Branch("pmtPhiQnscb_single", pmtPhiQnscb_single1,"pmtPhiQnscb_single[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQnscb_single", pmtCosThetaQnscb_single1, "pmtCosThetaQnscb_single[NsortPMTs]/F")
tree1.Branch("pmtQnscb_single", pmtQnscb_single1, "pmtQnscb_single[NsortPMTs]/F")
tree1.Branch("pmtPhiQnscb_doubleA", pmtPhiQnscb_doubleA1, "pmtPhiQnscb_doubleA[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQnscb_doubleA", pmtCosThetaQnscb_doubleA1,"pmtCosThetaQnscb_doubleA[NsortPMTs]/F")
tree1.Branch("pmtQnscb_doubleA", pmtQnscb_doubleA1, "pmtQnscb_doubleA[NsortPMTs]/F")
tree1.Branch("pmtPhiQnscb_doubleB", pmtPhiQnscb_doubleB1, "pmtPhiQnscb_doubleB[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQnscb_doubleB", pmtCosThetaQnscb_doubleB1, "pmtCosThetaQnscb_doubleB[NsortPMTs]/F")
tree1.Branch("pmtQnscb_doubleB", pmtQnscb_doubleB1, "pmtQnscb_doubleB1[NsortPMTs]/F")
## --- qprompt method ---
tree1.Branch("pmtPhiQprompt_single", pmtPhiQprompt_single1,"pmtPhiQprompt_single[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQprompt_single", pmtCosThetaQprompt_single1, "pmtCosThetaQprompt_single[NsortPMTs]/F")
tree1.Branch("pmtQprompt_single", pmtQprompt_single1, "pmtQprompt_single[NsortPMTs]/F")
tree1.Branch("pmtPhiQprompt_doubleA", pmtPhiQprompt_doubleA1, "pmtPhiQprompt_doubleA[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQprompt_doubleA", pmtCosThetaQprompt_doubleA1,"pmtCosThetaQprompt_doubleA[NsortPMTs]/F")
tree1.Branch("pmtQprompt_doubleA", pmtQprompt_doubleA1, "pmtQprompt_doubleA[NsortPMTs]/F")
tree1.Branch("pmtPhiQprompt_doubleB", pmtPhiQprompt_doubleB1, "pmtPhiQprompt_doubleB[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQprompt_doubleB", pmtCosThetaQprompt_doubleB1, "pmtCosThetaQprompt_doubleB[NsortPMTs]/F")
tree1.Branch("pmtQprompt_doubleB", pmtQprompt_doubleB1, "pmtQprompt_doubleB1[NsortPMTs]/F")
## --- qaap method --- 
tree1.Branch("pmtPhiQaap_single", pmtPhiQaap_single1,"pmtPhiQaap_single[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQaap_single", pmtCosThetaQaap_single1, "pmtCosThetaQaap_single[NsortPMTs]/F")
tree1.Branch("pmtQaap_single", pmtQaap_single1, "pmtQaap_single[NsortPMTs]/F")
tree1.Branch("pmtPhiQaap_doubleA", pmtPhiQaap_doubleA1, "pmtPhiQaap_doubleA[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQaap_doubleA", pmtCosThetaQaap_doubleA1,"pmtCosThetaQaap_doubleA[NsortPMTs]/F")
tree1.Branch("pmtQaap_doubleA", pmtQaap_doubleA1, "pmtQaap_doubleA[NsortPMTs]/F")
tree1.Branch("pmtPhiQaap_doubleB", pmtPhiQaap_doubleB1, "pmtPhiQaap_doubleB[NsortPMTs]/F")
tree1.Branch("pmtCosThetaQaap_doubleB", pmtCosThetaQaap_doubleB1, "pmtCosThetaQaap_doubleB[NsortPMTs]/F")
tree1.Branch("pmtQaap_doubleB", pmtQaap_doubleB1, "pmtQaap_doubleB1[NsortPMTs]/F")

#######################################################

#### PMT results ####
## data/DEAP-3600/pmt_positions.py

## For MC, use pre-saved PMT data in pickle, avoid connecting to couchdb
##server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
##db = server["deapdb"]
##results = db.view('WebView/PMTPos', include_docs=True)

Phi = {} #np.zeros(255)
CosTheta = {}#np.zeros(255)
index = [] #np.zeros(255,dtype = int)

#  x = math.sin(math.radians(row.doc["locationTheta"]))*math.cos(math.radians(row.doc["locationPhi"]))
#  y = math.sin(math.radians(row.doc["locationTheta"]))*math.sin(math.radians(row.doc["locationPhi"]))
#  z = math.cos(math.radians(row.doc["locationTheta"]))

## buil PMT look-up tables
fpmtdata = open("savePMTinfo.pickle","r")
my_dict = pickle.load(fpmtdata)
for row in my_dict:
    key = row.keys()[0]
    Phi[key] = row[key][0]
    CosTheta[key] = row[key][1]
    index.append(key)

#for row in results:
#	if row.doc["run_range"][0] == 0:
#		 #print(int(row.key), row.doc["locationPhi"])
#		 Phi[int(row.key)] = row.doc["locationPhi"]*PI/180
#		 CosTheta[int(row.key)] = np.cos(row.doc["locationTheta"]*PI/180)
#		 index[int(row.key)] = int(row.key)

#PMTlocation = []  ## Cartesian coordination
#for x in zip(CosTheta, Phi):
#    PMTlocation.append(TVector3(np.sqrt(1.0-x[0]**2)*np.cos(x[1]),np.sqrt(1.0-x[0]**2)*np.sin(x[1]),x[0]))
#
#zz = [] ## save PMT id if cos()<16.3 degree
#for x in PMTlocation:
#    zz.append([y[1] for y in zip(PMTlocation,index) if x.Dot(y[0])>np.cos(16.3*PI/180)])
#print(zz)
n1 = []
n2 = []
n3 = []

#def addPMTs(id,cluster,q_histo,t_histo):
#	global qc
#	for k in zz[id]:
#           if charges[k] > 0 and (k not in cluster):
#              cluster.append(k)
#              qc = qc + charges[k]
#              t_histo.Fill(Phi[k],CosTheta[k],pmtTime[k])
#              q_histo.Fill(Phi[k],CosTheta[k],charges[k])
#              if charges[k] > 10: # find PMT q>10 pe,
#                 addPMTs(k,cluster,q_histo,t_histo)

def search_cluster( listPMTid ): ## sort 6 highest-charge PMTs
    return tag  ##





countTrig = 0
countTrig_afterCuts = 0
print "\n********************************** File " , fin, " loaded... **********************************\n"
# data = File.Get("T")
data = fin.Get("T_satCorr");#TTree
nentries = data.GetEntries();
update = int(nentries/100);
for event in range(nentries):
    #if (event+1)%update == 0:
    #    print(event, "analyzed...")
    data.GetEntry(event);DS = data.ds;#(required)
    if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
        TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);
        try: CalTrigTime = CAL.GetEventTime()
        except: continue
        countTrig += 1
        # Low Level cuts. NOTE: qPE>60 cuts for Cherenkov threshold. No Fmaxpe cuts.
        if ( CAL.GetCuts().GetCutWord()&0x31f8 ):  continue; # Low level cut, NOTE: use 0x199 if looking at data between May-Aug 2016 processed in 2016
        if ( TS.dtmTrigSrc&0x82                ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
        if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
        if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
        if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
        if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
        if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
        #if ( EV.GetFmaxpe() > 0.4              ):  continue; #
        #if ( CAL.GetFprompt() > 0.55           ):  continue; # was >0.58
        if ( CAL.GetQPE() <= 60               ):  continue; # Skip events that have no charge information (primarily pre-scaled events)
        H2_CosTheta_Phi.Reset();
        T2_CosTheta_Phi.Reset();

        eventID = TS.GetEventID();
        evtID[0] = eventID
        runID[0] = DS.GetRunID();
        subrunID[0] = DS.GetSubrunID();

        numEarlyPulsesVal[0] = CAL.numEarlyPulses;
        subeventNval[0] = CAL.GetSubeventCount();
        nhitVal[0] = int(ord(CAL.GetLateNhit()))

        ##print "!!!!! initial count", countTrig_afterCuts
        eventTimeVal[0] = CalTrigTime # daniel 6th cut: eventTime > 2250 and eventTime < 2700
        fprompt = CAL.GetFprompt(); fpromptVal[0] = fprompt
        qpe = CAL.GetQPE(); qpeVal[0] = qpe
        fmaxpeVal[0] = EV.GetFmaxpe();

        ### get event position (MBLikelihood), qpe, fpromt, etc
        posx = -9999; posy = -9999; posz = -9999;
        tf2x = -9999; tf2y = -9999; tf2z = -9999;
        #posPhi = -9999; posCosTheta = -9999;
        time_tf2_ch = -9999
        ### Event Reconstruction MB likelihood valid
        try:
            posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
            #posPhi = EV.mblikelihood.GetPosition().Phi(); posCosTheta = EV.mblikelihood.GetPosition().CosTheta();
        except:
            posx = -9999; posy = -9999; posz = -9999;
        
        mbx[0] = posx; mby[0] = posy; mbz[0] = posz;
        ### Event Reconstruction of TF2 valid
        try:
            tf2x = EV.timefit2.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
        except:
            tf2x = -9999; posy = -9999; posz = -9999;
        try:
            time_tf2_ch = EV.timefit2.ch_t
        except:
            time_tf2_ch = -9999

        tf2_ch_tVal[0] = time_tf2_ch

        ## Fred Schuckman cuts
        qTopRing = EV.GetChargeTopRing(); qSecondRing = EV.GetChargeSecondRing();
        qBotRing = EV.GetChargeBottomRing(); qSecondBotRing = EV.GetChargeSecondBottomRing(); qThirdBotRing = EV.GetChargeThirdBottomRing();
        ### for liquid level, just save values, for further cuts!!!
        pulseindexfirstgarVal[0] =EV.pulseindexfirstgar
        cft2r[0] = (qTopRing + qSecondRing)/qpe
        cfb3r[0] = (qBotRing + qSecondBotRing + qThirdBotRing)/qpe

        nSCBayesVal[0] = CAL.GetNSCBayes()
        rprompt60BayesVal[0] = EV.GetRprompt60Bayes()

        scintlikeVal[0] = EV.GetScintillationLikelihood()
        scintlike16Val[0] = EV.GetScintillationLikelihood16ns()
        scintlikeqVal[0] = EV.GetChargeGOF()
        llneutronflashVal[0] = EV.GetLLNeutronFlash()
        #qPEnoSat_10000Val[0] = CAL.GetQPEnoSat_10000()
        #qPEnoSat_5000Val[0] = CAL.GetQPEnoSat_5000()
        neckVetoVal[0] = CAL.GetNeckVetoPMTCount()

        ### check pmts, instead of np.zeros, use {} dictionary
        charges_qpe = {}; charges_qnscb = {}; charges_qprompt = {}; charges_qaap = {};
        pmtTime = {}; ### time of each PMT
        savePmtPhi = {}; savePmtCosTheta = {}; 
        pmtTime1 = {}; pmtTime2 = {}; 
        pmtInfo.SetRunAndSubrun(runID[0], subrunID[0])
        saveAllSubpeakTime = {} ## for each pmtID, save all pulses and subpeak time
        nHitPMTs = CAL.GetPMTCount()
        nPMTs1[0] = nHitPMTs 
        for ipmt in range (nHitPMTs):### loop PMTs
            PMT = CAL.GetPMT(ipmt)
            PMTid = PMT.id
            ## NOTE: select good pmts, Only for data, DONOT use for MC!!
            if pmtInfo.IsPMTGoodForAnalysis(PMT) == False: continue ### !!! added 15 Nov 2022
            spe = pmtInfo.GetSPECharge(PMT) ### !!! added 15 Nov 2022
            if not(spe>0): continue ### !!! added 15 Nov 2022
            ### look at 4 charges!!
            pmt_qPE = PMT.qPE
            pmt_qNSCB   = PMT.nscBayes
            pmt_qPrompt = PMT.qPrompt #PMT.qPrompt/spe
            pmt_qAAP = pmt_qPE - PMT.GetQAP() ## to remove the Alberta AfterPulsing charge
            ## ensure all the 4 charges exist
            if not(pmt_qPE>0): continue
            if not(pmt_qPrompt>0): continue
            if not(pmt_qNSCB>0): continue
            if not(pmt_qAAP>0): continue 

            #specharge = PMTInfoUtil().GetSPECharge(PMTid) # James
            ###NOTE: default algorithm used PMT.qPrompt, without /spe.
            charges_qpe[PMTid] = pmt_qPE 
            charges_qnscb[PMTid] = pmt_qNSCB
            charges_qpromt[PMTid] = pmt_qPrompt
            charges_qaap[PMTid] = pmt_qAAP
            listPMTid[ipmt] = PMTid

            pmtTime[PMTid] = 100000
            pmtTime1[PMTid] = 100000
            pmtTime2[PMTid] = 100000

            ### three types of subpeak time
            saveSubpeakTime = []
            saveSubpeakTime1 = []
            saveSubpeakTime2 = []
            nPulses = PMT.GetPulseCount()
            for j in range (nPulses): ## loop subpeak for a PMT
                PULSE = PMT.GetPulse(j)
                ### check for nPulses!!
                for ipeak in range(PULSE.GetSubpeakCount()): ## loop subpeaks for a PMT
                    #print "!!!! PMT id", PMTid, "pulse:", j, "subpeak:", ipeak
                    time = PULSE.GetSubpeakTime(ipeak) ##
                    time1 = PULSE.GetSubpeakTimeUniform(ipeak)
                    time2 = PULSE.GetSubpeakTimeUniformWithTOFCorrection(ipeak)
                    if time < pmtTime[PMTid]:
                        saveSubpeakTime.append(time)
                    if time1 < pmtTime1[PMTid]:
                        saveSubpeakTime1.append(time1)
                    if time2 < pmtTime2[PMTid]:        
                        saveSubpeakTime2.append(time2)

            saveAllSubpeakTime[PMTid] = saveSubpeakTime ## here only check PULSE.GetSubpeakTime
            ### always using earliest subpeak time as the pmttime !! added 10 Nov 2022
            if len(saveSubpeakTime)>0:
               pmtTime[PMTid] = min(saveSubpeakTime)
            else: pmtTime[PMTid] = 0
            if len(saveSubpeakTime1)>0:
               pmtTime1[PMTid] = min(saveSubpeakTime1)
            else: pmtTime1[PMTid] = 0
            if len(saveSubpeakTime2)>0:
               pmtTime2[PMTid] = min(saveSubpeakTime2)
            else: pmtTime2[PMTid] = 0

            H2_CosTheta_Phi.Fill(Phi[PMTid], CosTheta[PMTid], q)
            T2_CosTheta_Phi.Fill(Phi[PMTid], CosTheta[PMTid], pmtTime[PMTid])

        ### look for cluster 1
        cluster1 = []
        qcluster1 = 0
        ## Sort the PMT ids based on their charge in an decreasing order; i.e first item is the brightest PMT
        sortPMT = sorted(charges.items(), key=lambda x:x[1], reverse=True)
        nCheckPMTs = len(sortPMT)
        ### qsort = [(pmtID0, q0), (pmtID1, q1), ...]

        if nCheckPMTs<1:
            print "no PMT selected, pass this event !!!"
            count_ncluster = 0
            continue

        #qsort_pmtid, qSorted = [], []
        #for pmtdata in sortPMT:
        #     ## print "!!! pmt data", pmtdata[0], pmtdata[1]
        #     qsort_pmtid.append(pmtdata[0])
        #     qSorted.append(pmtdata[1])

        ### select the 6 largest Q pmts (from 0 to 5), and q5>20
        nCheck = 6
        if nCheckPMTs<nCheck:
            nCheck = nCheckPMTs
        ## ensure q5>20, so not necessary 6 PMTs!
        q_selected, id_selected  = [], []
        for i in range(nCheck):
            qq = sortPMT[i][1]
            idx = sortPMT[i][0]
            if qq>QminCut:
                q_selected.append(sortPMT[i][1])
                id_selected.append(sortPMT[i][0])

        nSelected = len(q_selected) ## actually selected PMTs, not necessarily 6

        if nSelected<1:
            #print "too low PMT charges, no obvious hit points, pass this event !!!"
            count_ncluster = 0
            continue

        if nSelected<2:
            #print "just a single high-charge PMT, count as single-cluster, done this event!!!"
            count_ncluster = 1
            continue

        #print "sort list of pmt q", selectQ_pmtList
        qmax = q_selected[0] ## the 1st is the largest
        qmin = q_selected[nSelected-1]

        if qmax<QmaxCut:
            #print "too low PMT charges, no obvious hit points, pass this event !!!"
            count_ncluster = 0
            continue

        select_pmtPos = [] ## from 0 to nSelected
        for i in range(nSelected):
            idx = id_selected[i]
            phi = Phi[idx] 
            cost = CosTheta[idx]
            th = arccos(cost)
            pmtPos = TVector3(radius*cos(phi)*sin(th),radius*sin(phi)*sin(th), radius*cost)
            select_pmtPos.append(pmtPos)

        pmtPos_Qmax = select_pmtPos[0]
        #print "======= checking event", evtID[0], "qmax", selectQ_pmtList

        flag_single = True ## by default as single-cluster
        flag_double = False
        flag_triple = False
        count_nearPMT = 0
        cluster_phi = []
        cluster_cosTheta = []
        ### search PMT id by looking up in charges[]
        searchPMTid_cluster1 = [] # could be single-cluster or the 1st of the double
        searchPMTid_cluster2 = [] # the 2nd of the double
        searchPMTid_cluster3 = [] # the 3rd of the triple
        searchPMTid_cluster1.append(id_selected[0]) # save the qmax first; this is the seed of the 1st cluster
        seed_cluster2_idx = -10000 # if there is a PMT far away from qmax(i.e q0) PMT, is this PMT the seed of the 2nd cluster?

        ## First check the PMTs around q0 with dist, angles between 0,1; 0,2; ..., each q PMT > 20 pC
        for i in range(1,nSelected): # nSelected>1 is ensured previously, then compare q0 to q1, ..., q5
            dist = (select_pmtPos[i] - pmtPos_Qmax).Mag()
            angle = pmtPos_Qmax.Angle(select_pmtPos[i])
            if angle<(16.3*pi/180) and dist<250:
                count_nearPMT += 1
                searchPMTid_cluster1.append(id_selected[i]) # Note: selectQ_pmtList from 1, 2, ... instead of 0
                continue
            if angle>(16.3*pi/180) and dist>250:
                # print "eventID", evtID[0],"double cluster! dist = ", dist
                flag_double = True
                ## save PMT charges in 2nd cluster
                searchPMTid_cluster2.append(id_selected[i]) ## Note: !!!selectQ_pmtList from 1, 2, ... instead of 0
                seed_cluster2_idx = i # this index is for the selected PMT list, from 1 to 5, not the PMT id
                #print "event ", evtID[0], "there is a seed2 for double cluster, it is the q%d PMT" %seed_cluster2_idx, "with charge %f pC" %selectQ_pmtList[seed_cluster2_idx]
                break;## break the loop around the seed q0, done with single-cluster and then check the loop around the new seed

        #print "checking!!!!!!!"
        #print "selectQ_pmtList", selectQ_pmtList
        #print "pmt0_angles", pmt0_angles

        ## check the 2nd cluster from the seed
        searchPMTid_cluster3_candidate = [] # save the ids from 1 to 5 in select_pmtPos
        if flag_double == True:
             pos_seed_cluster2 = select_pmtPos[seed_cluster2_idx]
             for i in range(seed_cluster2_idx+1, nSelected):# automatically ensures seed_cluster2_idx <= nSelected
                   ### first check whether the PMTs after this 2nd seed are around the q0
                   angle0 = pmtPos_Qmax.Angle(select_pmtPos[i])
                   dist0 = (pmtPos_Qmax - select_pmtPos[i]).Mag()
                   ## is it around q0 or belonging to single cluster?
                   if angle0<(16.3*pi/180) and dist0<250:
                       count_nearPMT += 1
                       searchPMTid_cluster1.append(id_selected[i])
                   ## then check whether it is around the 2nd seed
                   angle = pos_seed_cluster2.Angle(select_pmtPos[i]) ## angle between TVector3
                   dist = (select_pmtPos[i] - pos_seed_cluster2).Mag()
                   if angle<(16.3*pi/180) and dist<250:
                        searchPMTid_cluster2.append(id_selected[i])
                   if angle>(16.3*pi/180) and dist>250: ## is there a possible triple-cluster?
                        if angle0>(16.3*pi/180) and dist0>250: ## it is also not surround cluster1
                            flag_triple = True # if there is a PMT far away from the 2nd seed, then could be triple/multiple, needn't check the rest
                            searchPMTid_cluster3_candidate.append(i)
                            
        if flag_triple == True:
             # print "a possible triple??"
             seed_cluster3_idx = searchPMTid_cluster3_candidate[0]
             pos_seed_cluster3 = select_pmtPos[seed_cluster3_idx] #position of the seed of cluster 3
             searchPMTid_cluster3.append(id_selected[i]) ## save the seed of cluster 3
             for i in range(seed_cluster3_idx + 1, nSelected):
                   angle3 = pos_seed_cluster3.Angle(select_pmtPos[i]) 
                   dist3 = (pos_seed_cluster3 - select_pmtPos[i]).Mag()
                   if angle3<(16.3*pi/180) and dist3<250:
                            searchPMTid_cluster3.append(id_selected[i])

        ## check the summed charges of the clusters
        check_cluster1_charges = []
        check_cluster2_charges = []
        check_cluster3_charges = []
        ## for double-clusters, it requires the sum of the PMTs> 100 pC
        #print "searched charges in cluster2", searchPMTid_cluster2
        if len(searchPMTid_cluster2) != 0:
             qsum_double = 0
             ## qsum_double = charges[searchPMTid_cluster2[0]] #just look at the first seed
             for i in searchPMTid_cluster2:## selected pmtIds
                  # print "charges in cluster2", i, charges[i]
                  qsum_double =  qsum_double + charges[i] ## look up the charge in the original PMT loop
             if qsum_double<QmaxCut:
                  #print "pmt charge sum %f<QmaxCut="%qsum_double, QmaxCut," pC, no double cluster any more!"
                  flag_double = False
        else:
             flag_double = False

        ## then check triple-cluster
        if len(searchPMTid_cluster3) != 0:
             qsum_triple = 0
             for i in searchPMTid_cluster3:
                  # print "charges in cluster2", i, charges[i]
                  qsum_triple =  qsum_triple + charges[i]
             if qsum_triple<QmaxCut:
                  #print "pmt charge sum %f<QmaxCut="%qsum_triple, QmaxCut,"pC, no triple cluster any more!"
                  flag_triple = False
             else:
                  #print "catch a triple/multiple cluster! very rare, cluster1",searchPMTid_cluster1, "cluster2", searchPMTid_cluster2, "cluster3", searchPMTid_cluster3
                  #print "charge of 3rd is", qsum_triple
                  flag_triple = True
                  flag_double = False
                  #for ind in searchPMTid_cluster3:
                  #     check_cluster3_charges.append(charges[ind])
        else:
             flag_triple = False

        if flag_double == True or flag_triple == True:
             flag_single = False
             # print "eventID", evtID[0], "single cluster!!"


        for idx in searchPMTid_cluster1:
             check_cluster1_charges.append(charges[idx])

        for idx in searchPMTid_cluster2:
             check_cluster2_charges.append(charges[idx])

        for idx in searchPMTid_cluster3:
             check_cluster2_charges.append(charges[idx])
             #check_cluster3_charges.append(charges[idx])

        ### Fill single-Cluster
        if flag_single == True:
             # print "evt", evtID[0], "is a single-cluster, with id", searchPMTid_cluster1, "charge:", check_cluster1_charges
             count_ncluster = 1 ## count for single cluster
             ncluster[0] = count_ncluster
             #id1 = charges.index(qmax)
             #ClusterPos_CosTheta_Phi1.Fill(Phi[id1],CosTheta[id1]) ## PMT pos, no weighting
             nPMTs1[0] = len(searchPMTid_cluster1)
             kk = 0
             # print "single", searchPMTid_cluster1
             for kk in range(nPMTs1[0]):
                    pmtid = searchPMTid_cluster1[kk]
                    listPMTid1[kk] = pmtid
                    pmtPhi1[kk] = Phi[pmtid]
                    pmtCosTheta1[kk] = CosTheta[pmtid]
                    pmttime1[kk] = pmtTime[pmtid]
                    pmttime_uni1[kk] = pmtTime1[pmtid] 
                    pmttime_unicor1[kk] = pmtTime2[pmtid]
                    pmtQpe1[kk] = charges[pmtid] #pmt.qPrompt
                    chargeNSCB1[kk] = charges_nscb[pmtid]
                    chargeQPE1[kk] = charges_qpe[pmtid]

                    for time in saveAllSubpeakTime[pmtid]:
                        HsubpkTime_cluster1.Fill (time)
                        HsubpkTime_qpromptWt_cluster1.Fill( time, pmtQpe1[kk] )
                        HsubpkTime_qnscbWt_cluster1.Fill( time, chargeNSCB1[kk] )
                        HsubpkTime_qpeWt_cluster1.Fill( time, chargeQPE1[kk] )
                        HtimeRes_cluster1.Fill( time-CalTrigTime )
                        HtimeRes_qpromptWt_cluster1.Fill( time-CalTrigTime, pmtQpe1[kk])
                        HtimeRes_qnscbWt_cluster1.Fill( time-CalTrigTime, chargeNSCB1[kk])
                        HtimeRes_qpeWt_cluster1.Fill( time-CalTrigTime, chargeQPE1[kk])

                    # tagCluster1[kk] = 1
                    # print "pmtid", pmtid, kk, pmtPhi1[kk], pmtCosTheta1[kk], pmttime1[kk], pmtQpe1[kk], pmttime_unicor1[kk]
             ## print "!!!!! tree1.Fill count", countTrig_afterCuts       
             tree1.Fill()
             Multiplicity.Fill(ncluster[0])

        if flag_double == True:
             # print "evt", evtID[0], "is a double-cluster, with 1st id", searchPMTid_cluster1, "charge:", check_cluster1_charges, "with 2nd id", searchPMTid_cluster2, "charge:", check_cluster2_charges
             ncluster[0] = 2
             #id1 = charges.index(qmax)
             #ClusterPos_CosTheta_Phi2.Fill(Phi[id1],CosTheta[id1]) ## PMT pos, no weighting
             nPMTs2[0] = len(searchPMTid_cluster1) ## the 1st cluster is around q0
             kk = 0 ## save the 1st cluster of double-cluster
             #print "double", searchPMTid_cluster1
             for kk in range(nPMTs2[0]):
                    pmtid = searchPMTid_cluster1[kk]
                    listPMTid2[kk] = pmtid
                    pmtPhi2[kk] = Phi[pmtid]
                    pmtCosTheta2[kk] = CosTheta[pmtid]
                    pmttime2[kk] = pmtTime[pmtid]
                    pmttime_uni2[kk] = pmtTime1[pmtid]
                    pmttime_unicor2[kk] = pmtTime2[pmtid]
                    charge2[kk] = charges[pmtid]
                    chargeNSCB2[kk] = charges_nscb[pmtid]
                    chargeQPE2[kk] = charges_qpe[pmtid]

                    for time in saveAllSubpeakTime[pmtid]:
                        HsubpkTime_cluster2a.Fill (time)
                        HsubpkTime_qpromptWt_cluster2a.Fill( time, charge2[kk] )
                        HsubpkTime_qnscbWt_cluster2a.Fill( time, chargeNSCB2[kk] )
                        HsubpkTime_qpeWt_cluster2a.Fill( time, chargeQPE2[kk] )
                        HtimeRes_cluster2a.Fill( time-CalTrigTime )
                        HtimeRes_qpromptWt_cluster2a.Fill( time-CalTrigTime, charge2[kk])
                        HtimeRes_qnscbWt_cluster2a.Fill( time-CalTrigTime, chargeNSCB2[kk])
                        HtimeRes_qpeWt_cluster2a.Fill( time-CalTrigTime, chargeQPE2[kk])

                    # tagCluster2[kk] = 2
                    # print "pmtid", pmtid, pmtPhi2[kk], pmtCosTheta2[kk], pmttime2[kk]
             ## print "!!!!! tree2.Fill count", countTrig_afterCuts
             tree2.Fill()

             nPMTs3[0] = len(searchPMTid_cluster2) ## the 2nd cluster is around qi
             kk = 0 ## save the 2nd cluster of double-cluster
             for kk in range(nPMTs3[0]):
                    pmtid = searchPMTid_cluster2[kk]
                    listPMTid3[kk] = pmtid
                    pmtPhi3[kk] = Phi[pmtid]
                    pmtCosTheta3[kk] = CosTheta[pmtid]
                    pmttime3[kk] = pmtTime[pmtid]
                    pmttime_uni3[kk] = pmtTime1[pmtid]
                    pmttime_unicor3[kk] = pmtTime2[pmtid]
                    charge3[kk] = charges[pmtid]
                    chargeNSCB3[kk] = charges_nscb[pmtid]
                    chargeQPE3[kk] = charges_qpe[pmtid]

                    for time in saveAllSubpeakTime[pmtid]: 
                        HsubpkTime_cluster2b.Fill (time)
                        HsubpkTime_qpromptWt_cluster2b.Fill( time, charge3[kk] )
                        HsubpkTime_qnscbWt_cluster2b.Fill( time, chargeNSCB3[kk] )
                        HsubpkTime_qpeWt_cluster2b.Fill( time, chargeQPE3[kk] )
                        HtimeRes_cluster2b.Fill( time-CalTrigTime )
                        HtimeRes_qpromptWt_cluster2b.Fill( time-CalTrigTime, charge3[kk] )
                        HtimeRes_qnscbWt_cluster2b.Fill( time-CalTrigTime, chargeNSCB3[kk] )
                        HtimeRes_qpeWt_cluster2b.Fill( time-CalTrigTime, chargeQPE3[kk] )
                    # tagCluster3[kk] = 2
             ## print "!!!!! tree3.Fill count", countTrig_afterCuts       
             tree3.Fill()
             Multiplicity.Fill(ncluster[0])

        if flag_triple == True:
             #print "!!!Ok, get a triple-cluster, with 1st id", searchPMTid_cluster1, "charge:", check_cluster1_charges, "with 2nd id", searchPMTid_cluster2, "charge:", check_cluster2_charges, "with 3rd id", searchPMTid_cluster3, "charge:", check_cluster3_charges
             ncluster[0] = 3
             Multiplicity.Fill(ncluster[0])

## end event loop       
fout.cd()
tree1.Write()
tree2.Write()
tree3.Write()
Multiplicity.Write("nclusters")
## save subpeak times!!
HsubpkTime_cluster1.Write()
HsubpkTime_cluster2a.Write()
HsubpkTime_cluster2b.Write()

HsubpkTime_qpromptWt_cluster1.Write()
HsubpkTime_qpromptWt_cluster2a.Write()
HsubpkTime_qpromptWt_cluster2b.Write()

HsubpkTime_qnscbWt_cluster1.Write()
HsubpkTime_qnscbWt_cluster2a.Write()
HsubpkTime_qnscbWt_cluster2b.Write()

HsubpkTime_qpeWt_cluster1.Write()
HsubpkTime_qpeWt_cluster2a .Write()
HsubpkTime_qpeWt_cluster2b.Write()

HtimeRes_cluster1.Write()
HtimeRes_cluster2a.Write()
HtimeRes_cluster2b.Write()

HtimeRes_qpromptWt_cluster1.Write()
HtimeRes_qpromptWt_cluster2a.Write()
HtimeRes_qpromptWt_cluster2b.Write()

HtimeRes_qnscbWt_cluster1.Write()
HtimeRes_qnscbWt_cluster2a.Write()
HtimeRes_qnscbWt_cluster2b.Write()

HtimeRes_qpeWt_cluster1.Write()
HtimeRes_qpeWt_cluster2a.Write()
HtimeRes_qpeWt_cluster2b.Write()

#_%d"%event
##########################################################################################################################################
##########################################################################################################################################
fin.Close()

print "***************************"
print "total trigger event count", countTrig
print "triggered event after cuts", countTrig_afterCuts

fout.Close()
