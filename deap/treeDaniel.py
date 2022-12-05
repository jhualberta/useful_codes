## for analyzing vacuum data to check Th-232 source, using Daniel's cuts
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

global qc

PI = np.pi
radius = 850

#runs_name = ["U232Source_SLB005_2020_CalC_L0", "U232Source_SLB005_2020_CalC_RT465_L0", "U232Source_SLB005_2020_CalC_RT480_L0", "U232Source_SLB005_2020_CalF_L0", "U232Source_SLB005_2020_CalFpos6_RT465_L0", "U232Source_SLB005_2020_CalFpos6_RT480_L0", "U232Source_SLB005_2020_CalFpos7_RT465_L0", "U232Source_SLB005_2020_CalFpos7_RT480_L0", "U232Source_SLB005_2021_CalFpos6_RT480_L0"]
#
#RunList = runs
#
#Run            = sys.argv[1]

#### Loading PMT results ####
## data/DEAP-3600/pmt_positions.py
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
Dir=sys.argv[1]
AllFiles=glob.glob(Dir)
AllFiles = sorted(AllFiles)

results = db.view('WebView/PMTPos', include_docs=True)
Phi = np.zeros(255)
CosTheta = np.zeros(255)
pmtX = np.zeros(255)
pmtY = np.zeros(255)
pmtZ = np.zeros(255)

index = np.zeros(255,dtype = int)

#  x = math.sin(math.radians(row.doc["locationTheta"]))*math.cos(math.radians(row.doc["locationPhi"]))
#  y = math.sin(math.radians(row.doc["locationTheta"]))*math.sin(math.radians(row.doc["locationPhi"]))
#  z = math.cos(math.radians(row.doc["locationTheta"]))


## buil PMT look-up tables
for row in results:
        if row.doc["run_range"][0] == 0:
                 ind = int(row.key)
                 index[ind] = ind
                 phi = row.doc["locationPhi"]*PI/180
                 cost = np.cos(row.doc["locationTheta"]*PI/180)
                 Phi[ind] = phi 
                 CosTheta[ind] = cost 
                 pmtX[ind] = radius*cos(phi)*sin(arccos(cost))
                 pmtY[ind] = radius*sin(phi)*sin(arccos(cost))
                 pmtZ[ind] = radius*cost

#PMTlocation = []  ## Cartesian coordination
#for x in zip(CosTheta, Phi):
#        PMTlocation.append(TVector3(np.sqrt(1.0-x[0]**2)*np.cos(x[1]),np.sqrt(1.0-x[0]**2)*np.sin(x[1]),x[0]))

file0 = str(sys.argv[1])
fileName = os.path.basename(file0)
###########################################################################################################################
#Files = AllFiles[subrun];

File = TFile(file0)
print ("\n********************************** File " , File , " loaded... **********************************\n")
data = File.Get("T_satCorr")
nentries = data.GetEntries();           
update = int(nentries/10);
fout = TFile("Daniel_level4Cuts_"+fileName, "RECREATE")
###########################################################################################################################
#foutfile="/project/6004969/dpapi/cherenkov/skim/{}/{}".format(runs_name[int(Run)-1], value.split("/")[-1])
#fout = TFile(foutfile,"RECREATE")
#dstreeclone = data.CloneTree(0);
#dstreeclone.SetDirectory(fout);
###########################################################################################################################

##### Tree information #####
### save as different cluster events
tree1 = TTree("T","processed")

### these are for each event, values of variables
runID = array('l',[0])
subrunID = array('l',[0])
evtID = array('l',[0])
evtx = array('f',[0])
evty = array('f',[0])
evtz = array('f',[0])
qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])
fmaxpeVal = array('f',[0])
nSCBayesVal = array('f',[0])
rprompt60BayesVal = array('f',[0])

dtmTrigSrcVal = array('f',[0])
subeventNVal = array('f',[0])
eventTimeVal = array('f',[0])
calCutsVal = array('f',[0])
deltatVal  = array('f',[0])
numEarlyPulsesVal = array('f',[0])

Nmax = 255 ## save how many PMTs for an event
NpeaksMax = 2550

#### pmt by pmt information
nPMTs = array('i',[0])
nSubpeaks = array('i',[0])
pmtphi = array('f',Nmax*[0])
pmtcosTheta = array('f',Nmax*[0])

pmtposx = array('f',Nmax*[0])
pmtposy = array('f',Nmax*[0])
pmtposz = array('f',Nmax*[0])

pmttime = array('f',NpeaksMax*[0]) # SubpeakTime
pmttime1 = array('f',NpeaksMax*[0]) # SubpeakTime
pmttime2 = array('f',NpeaksMax*[0]) # SubpeakTime

pmtq = array('f',Nmax*[0])

tree1.Branch("runID",runID,"runID/L")
tree1.Branch("subrunID",subrunID,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("nSCBayes", nSCBayesVal, "nSCBayes/F") # GetNSCBayes()
tree1.Branch("rprompt60Bayes", rprompt60BayesVal, "rprompt60Bayes/F")
tree1.Branch("fmaxpe", fmaxpeVal, "fmaxpe/F")

tree1.Branch("subeventN",subeventNVal, "subeventN/F")
tree1.Branch("eventTime",eventTimeVal, "eventTime/F")
tree1.Branch("dtmTrigSrc", dtmTrigSrcVal, "dtmTrigSrc/F")
tree1.Branch("calCuts",calCutsVal, "calCuts/F")
tree1.Branch("deltat",deltatVal, "deltat/F")
tree1.Branch("numEarlyPulses",numEarlyPulsesVal, "numEarlyPulses/F")

#tree1.Branch("evtx", evtx, "evtx/F")
#tree1.Branch("evty", evty, "evty/F")
#tree1.Branch("evtz", evtz, "evtz/F")
tree1.Branch("nPMTs",nPMTs, "nPMTs/I")
tree1.Branch("nSubpeaks",nSubpeaks, "nSubpeaks/I")

tree1.Branch("pmtq", pmtq, "pmtq[nPMTs]/F")
tree1.Branch("pmtPhi", pmtphi, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtcosTheta, "pmtCosTheta[nPMTs]/F")
#tree1.Branch("pmtPosX", pmtposx, "pmtPosX[nPMTs]/F")
#tree1.Branch("pmtPosY", pmtposy, "pmtPosY[nPMTs]/F")
#tree1.Branch("pmtPosZ", pmtposz, "pmtPosZ[nPMTs]/F")

tree1.Branch("pmttime", pmttime, "pmttime[nSubpeaks]/F")
#tree1.Branch("pmttime1", pmttime1, "pmttime1[nSubpeaks]/F")
#tree1.Branch("pmttime2", pmttime2, "pmttime2[nSubpeaks]/F")

#tree1.Branch("pmttime_uni", pmttime_uni1, "pmttime_uni[nPMTs]/F")
#tree1.Branch("pmttime_unicor", pmttime_unicor1, "pmttime_unicor[nPMTs]/F")

for event in range (nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."

    data.GetEntry(event); DS = data.ds;#(required)
    if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
        TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);
        try:CalTrigTime = CAL.GetEventTime()
        except:continue
        
        ### have to put this to reduce data !! Got to level 4 cuts
        if ( TS.dtmTrigSrc&0x82 ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
        if ( CAL.GetCuts().GetCutWord()&0x31f8 ): continue;
        if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
        if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
        #if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
        #if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
        #if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
        #if ( CAL.GetQPE() <= 60               ):  continue; # Skip events that have no charge information (primarily pre-scaled events)

        eventID = TS.GetEventID();
        if eventID != 363842: continue
        evtID[0] = eventID
        runID[0] = DS.GetRunID();
        subrunID[0] = DS.GetSubrunID();
        nSCBayesVal[0] = CAL.GetNSCBayes()
        fpromptVal[0] = CAL.GetFprompt()
        rprompt60BayesVal[0] = EV.GetRprompt60Bayes()

        dtmTrigSrcVal[0] = TS.dtmTrigSrc # daniel 1st cut: dtmTrigSrc&0x82 == 0
        calCutsVal[0] = CAL.GetCuts().GetCutWord() # daniel 2nd cut: calcut&0x31f8 == 0
        deltatVal[0] = EV.deltat # daniel 3rd cut: deltat > 20000
        numEarlyPulsesVal[0] = CAL.numEarlyPulses # daniel 4th cut: numEarlyPulses <= 3
        subeventNVal[0] = CAL.GetSubeventCount() # daniel 5th cut: subeventN==1 
        eventTimeVal[0] = CAL.GetEventTime() # daniel 6th cut: eventTime > 2250and eventTime < 2700
        qpeVal[0] = CAL.GetQPE() # daniel 7th cut: qPE>60
        fmaxpeVal[0] = EV.GetFmaxpe()
        #posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
        #evtx[0] = posx; evty[0] = posy; evtz[0] = posz;

        charges = []; pmtid = []; 
        pmtTime = []; ### time of each PMT
        pmtTime1 = []; pmtTime2 = []; 
        pmtPhi = []; pmtCosTheta = [];
        pmtPosX = []; pmtPosY = []; pmtPosZ = [];

        ## print "!!!! pmtCount", CAL.GetPMTCount()
        countPMTzeroQ = 0
        for i in range (CAL.GetPMTCount()): ### loop PMTs
             PMT = CAL.GetPMT(i)
             PMTid = PMT.id
             q = PMT.qPrompt #haojie
             ### just keep PMTs with charge>0
             if q == 0 or q<0:
                 countPMTzeroQ += 1
                 continue
             # specharge = PMTInfoUtil().GetSPECharge(PMTid) # James
             charges.append(q)
             pmtPhi.append(Phi[PMTid])
             pmtCosTheta.append(CosTheta[PMTid])
             #pmtPosX.append(pmtX[PMTid])
             #pmtPosY.append(pmtY[PMTid])
             #pmtPosZ.append(pmtZ[PMTid])
             #pmtid.append(PMTid)
             #print "event", event, q, Phi[PMTid], PMTid
             timeThresh = 100000 # 100 us
             ## print "!!!! pulseCount", PMT.GetPulseCount()
             for j in range (PMT.GetPulseCount()):
                  PULSE = PMT.GetPulse(j)
                  for ipeak in range(PULSE.GetSubpeakCount()):
                       time = PULSE.GetSubpeakTime(ipeak) ##
                       #time1 = PULSE.GetSubpeakTimeUniform(ipeak)
                       #time2 = PULSE.GetSubpeakTimeUniformWithTOFCorrection(ipeak)
                       if 0<time and time<timeThresh:
                               pmtTime.append(time)
                       #if time1 < timeThresh and time1 >0:
                       #        pmtTime1.append(time1)
                       #if time2 < timeThresh and time2 >0:
                       #        pmtTime2.append(time2)

             ## put conditions here; dimensions should be the same!!
             #charges_filter = filter(lambda x: x != 0, charges)
             #pmtPhi_filter = filter(lambda x: x != 0, pmtPhi)

        ## end of PMT loop
        if len(charges) != len(pmtPhi):
            print "pmt data dimension error, quit the event loop ..."
            continue
        #id_filter = []
        #for qq in charges_filter: ## non-zero PMT id
        #      idpmt = charges.index(qq)
        #      id_filter.append(idpmt)

        ## print "event", event
        ## print charges
        ## print pmtTime 
        ## print pmtCosTheta

        nPMTs[0] = len(charges)
        nSubpeaks[0] = len(pmtTime)
        ##print "!!! nPMTs", nPMTs[0], "zero Q pmts:", countPMTzeroQ, "nSubpeaks", nSubpeaks[0]
        ## must swtich list to array to save tree
        for i in range(nPMTs[0]):
              #indx = id_filter[i]
              pmtq[i] = charges[i]
              pmtphi[i] = pmtPhi[i]
              pmtcosTheta[i] = pmtCosTheta[i]
              #pmtposx[i] = pmtPosX[i]
              #pmtposy[i] = pmtPosY[i]
              #pmtposz[i] = pmtPosZ[i]
        for j in range(nSubpeaks[0]):
              pmttime[i] = pmtTime[i]
              #pmttime1[i] = pmtTime1[i]
              #pmttime2[i] = pmtTime2[i] 

        tree1.Fill()
#################################################################################################################################

#dstreeclone.Write()
#fout.Close()
#File.Close()

fout.cd()
tree1.Write()
fout.Close()
