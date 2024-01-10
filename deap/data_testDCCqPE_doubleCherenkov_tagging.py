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
from array import array
global qc
PI = np.pi
lastGArPMTID = 64     # The highest PMT ID with a light-guide area above the LAr fill level at z = 551 mm

## using pmt.qPE!!
# Default cut values
fportionCut = 0.30
fdistanceCut = 600.
fpmtFmaxpeCut = 0.06
nHiQpmt = 6
rAV = 851.

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
#check_eventID = str(sys.argv[2])

fileName = os.path.basename(file0)
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
Dir=sys.argv[1]

fout = TFile("CheckDcc_qPE_"+fileName,"RECREATE");
##NOTE: this is for data!
#pmtSPE = [9.401531, 9.478522, 9.533119, 9.598516, 9.51875, 9.603341, 9.511438, 9.463538, 9.478821, 9.378164, 9.459953, 9.446279, 9.588317, 9.465664, 9.639305, 9.462516, 9.474754, 9.442683, 9.480303, 9.496368, 9.602288, 9.431365, 9.54215, 9.519255, 9.537716, 9.390802, 9.560265, 9.593805, 9.538998, 9.510449, 9.570398, 9.487828, 9.549857, 9.637281, 9.604923, 9.557969, 9.57493, 9.645522, 9.489523, 9.434468, 9.637426, 9.577977, 9.456603, 9.460023, 9.5574, 9.631948, 9.572963, 9.75191, 9.609713, 9.455954, 9.579628, 9.659296, 9.540004, 9.493335, 9.549732, 9.559233, 9.523397, 9.516285, 9.605323, 9.494011, 9.551182, 9.589623, 9.518515, 9.405429, 9.564447, 9.575242, 9.662695, 9.448769, 9.579314, 9.587907, 9.492312, 9.54867, 9.557077, 9.543069, 9.638461, 9.59615, 9.642604, 9.512714, 9.612689, 9.577344, 9.521589, 9.427127, 9.589479, 9.54513, 9.599389, 9.552323, 9.549178, 9.561992, 9.699215, 9.605592, 9.466375, 9.541047, 9.42873, 9.425537, 9.38876, 9.405875, 9.56666, 9.667546, 9.495486, 9.40529, 9.505435, 9.374628, 9.473295, 9.627826, 9.534452, 9.575734, 9.510781, 9.610623, 9.548651, 9.516392, 9.4374, 9.532756, 9.543339, 9.673977, 9.531166, 9.377868, 9.469941, 9.608672, 9.551737, 9.532205, 9.448813, 9.41417, 9.413856, 9.423391, 9.538473, 9.559955, 9.596799, 9.320759, 7.98968, 9.619632, 9.465214, 9.501684, 9.330815, 9.449166, 9.592426, 9.643352, 9.675405, 9.441998, 9.369229, 9.509341, 9.417464, 9.510898, 9.705049, 9.413175, 9.588417, 9.663846, 9.504849, 9.457065, 9.626133, 9.35901, 9.543874, 9.991369, 9.56696, 9.674172, 9.504127, 9.564596, 9.010438, 9.424538, 9.680238, 9.506499, 9.587976, 9.579459, 9.57697, 9.589566, 9.458425, 9.424939, 9.545506, 9.706876, 9.580699, 9.421514, 9.621499, 9.322694, 9.405054, 9.422732, 9.752479, 9.612691, 9.613838, 9.383172, 9.686356, 9.648059, 9.417195, 9.591515, 9.526299, 9.427006, 9.615754, 9.570264, 9.560497, 9.646378, 9.615693, 9.522192, 9.646722, 9.64217, 9.558745, 9.522792, 9.378283, 9.497169, 9.389177, 9.47464, 9.416678, 9.215664, 9.493145, 9.562009, 9.442619, 9.425348, 9.466465, 9.403482, 9.461758, 9.509473, 9.510313, 9.494834, 9.50248, 9.317768, 9.59908, 9.629463, 9.543371, 9.514894, 9.703927, 9.613575, 9.511826, 9.359103, 9.375991, 9.439625, 9.557396, 9.564541, 9.625321, 9.521011, 9.224092, 9.63166, 9.650788, 9.463997, 9.229194, 9.455678, 9.653828, 9.601283, 9.588856, 9.61808, 9.6089, 9.581942, 9.561008, 9.512361, 9.59211, 9.13298, 9.625385, 9.561963, 9.186693, 9.485404, 9.510303, 9.481409, 9.394967, 9.436327, 9.438083, 9.448631, 9.636225, 9.36327, 9.562619]

results = db.view('WebView/PMTPos', include_docs=True)
Phi = {} #np.zeros(255)
CosTheta = {}#np.zeros(255)
index = {} #np.zeros(255,dtype = int)

## buil PMT look-up tables
for row in results:
    if row.doc["run_range"][0] == 0:
         #print(int(row.key), row.doc["locationPhi"])
         Phi[int(row.key)] = row.doc["locationPhi"]*PI/180
         CosTheta[int(row.key)] = np.cos(row.doc["locationTheta"]*PI/180)
         index[int(row.key)] = int(row.key)

pmtPos = {}
for pmtId in CosTheta.keys():
     phi = Phi[pmtId]
     costheta = CosTheta[pmtId]
     pmtPos[pmtId] = TVector3(np.sqrt(1.0-costheta*costheta)*np.cos(phi), np.sqrt(1.0-costheta*costheta)*np.sin(phi), costheta)

#pmtPos = []  ## Cartesian coordination
#for x in zip(CosTheta, Phi):
#    pmtPos.append(TVector3(np.sqrt(1.0-x[0]*x[0])*np.cos(x[1]), np.sqrt(1.0-x[0]*x[0])*np.sin(x[1]), x[0]))

runIDval = array('l',[0])
subrunIDval = array('l',[0])
evtID = array('l',[0])
mbx = array('f',[0])
mby = array('f',[0])
mbz = array('f',[0])

tf2x = array('f',[0])
tf2y = array('f',[0])
tf2z = array('f',[0])

qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])
nSCBayesVal = array('f',[0])
rprompt60BayesVal = array('f',[0])
fmaxpeVal = array('f',[0])
nhitVal = array('i',[0])
pulseGar = array('f',[0])
cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qpe<0.1
neckVetoVal = array('f',[0])

Nmax = 255
#### pmt by pmt, save cluster patterns
nPMTs1 = array('i',[0])
pmtPhi = array('f', Nmax*[0])
pmtCosTheta = array('f', Nmax*[0])
### warning!! here pmtQPE is not pmtCharge. pmtCharge is a dictionary {}
### for sorting the DCC algorithm
pmtQPE = array('f', Nmax*[0]) #pmt.qpe
pmtQNSCB = array('f', Nmax*[0])
pmtQPrompt = array('f', Nmax*[0])
listPMTid = array('i', Nmax*[0])

pmtFmaxpeVal = array('f',[0])
pmtFmaxpePromptVal = array('f',[0])

pmtPhi6 = array('f',nHiQpmt*[0])
pmtCosTheta6 = array('f',nHiQpmt*[0])
#pmttime6 = array('f',nHiQpmt*[0]) # SubpeakTime
#timeRes6 = array('f',nHiQpmt*[0])
pmtPhi1st = array('f',[0])
pmtCosTheta1st = array('f',[0])
pmtCharge1st = array('f',[0])
pmtPhi2nd = array('f',[0])
pmtCosTheta2nd = array('f',[0])
pmtCharge2nd = array('f',[0])

pmtCharge0 = array('f',[0])
pmtCharge1 = array('f',[0]) 
pmtCharge2 = array('f',[0]) 
pmtCharge3 = array('f',[0]) 
pmtCharge4 = array('f',[0])
pmtCharge5 = array('f',[0]) 

pmtChargePrompt0 = array('f',[0])
pmtChargePrompt1 = array('f',[0])
pmtChargePrompt2 = array('f',[0])
pmtChargePrompt3 = array('f',[0])
pmtChargePrompt4 = array('f',[0])
pmtChargePrompt5 = array('f',[0])

dbCherenkov_portion1 = array('f',[0]) 
dbCherenkov_portion2 = array('f',[0])
dbCherenkov_portion3 = array('f',[0])
dbCherenkov_portion4 = array('f',[0])
dbCherenkov_portion5 = array('f',[0])

dbCherenkov_dist1 = array('f',[0])
dbCherenkov_dist2 = array('f',[0])
dbCherenkov_dist3 = array('f',[0])
dbCherenkov_dist4 = array('f',[0])
dbCherenkov_dist5 = array('f',[0])

checkTagVal = array('I', [0])

tree1 = TTree("T1","saveEvent")
tree1.Branch("runID",runIDval,"runID/L")
tree1.Branch("subrunID",subrunIDval,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("fmaxpe", fmaxpeVal, "fmaxpe/F")
tree1.Branch("nSCBayes", nSCBayesVal, "nSCBayes/F")
tree1.Branch("rprompt60Bayes", rprompt60BayesVal, "rprompt60Bayes/F")
tree1.Branch("mbx", mbx, "mbx/F")
tree1.Branch("mby", mby, "mby/F")
tree1.Branch("mbz", mbz, "mbz/F")
tree1.Branch("tf2x", tf2x, "tf2x/F")
tree1.Branch("tf2y", tf2y, "tf2y/F")
tree1.Branch("tf2z", tf2z, "tf2z/F")
tree1.Branch("nPMTs1", nPMTs1, "nPMTs1/I")
tree1.Branch("pmtPhi", pmtPhi, "pmtPhi[nPMTs1]/F")
tree1.Branch("pmtCosTheta", pmtCosTheta, "pmtCosTheta[nPMTs1]/F")
tree1.Branch("pmtCharge", pmtQPE, "pmtQPE[nPMTs1]/F")
tree1.Branch("pmtChargeNSCB", pmtQNSCB, "pmtQNSCB[nPMTs1]/F")
tree1.Branch("pmtChargePrompt", pmtQPrompt, "pmtQPrompt[nPMTs1]/F")
tree1.Branch("listPMTid", listPMTid, "listPMTid[nPMTs1]/I")
tree1.Branch("neckVeto", neckVetoVal, "neckVeto/F")
tree1.Branch("pulseGar", pulseGar, "pulseGar/F")
tree1.Branch("cft2r", cft2r,"cft2r/F")
tree1.Branch("cfb3r", cfb3r,"cfb3r/F")
tree1.Branch("pmtFmaxpe", pmtFmaxpeVal, "pmtFmaxpe/F")
tree1.Branch("pmtFmaxpePrompt", pmtFmaxpePromptVal, "pmtFmaxpePrompt/F")

tree1.Branch("pmtPhi6", pmtPhi6, "pmtPhi6[6]/F") 
tree1.Branch("pmtCosTheta6", pmtCosTheta6, "pmtCosTheta6[6]/F")
tree1.Branch("pmtPhi1st", pmtPhi1st, "pmtPhi1st/F")
tree1.Branch("pmtCosTheta1st", pmtCosTheta1st, "pmtCosTheta1st/F")
tree1.Branch("pmtCharge1st", pmtCharge1st, "pmtCharge1st/F")
tree1.Branch("pmtPhi2nd", pmtPhi2nd, "pmtPhi2nd/F")
tree1.Branch("pmtCosTheta2nd", pmtCosTheta2nd, "pmtCosTheta2nd/F")
tree1.Branch("pmtCharge2nd", pmtCharge2nd, "pmtCharge2nd/F")
tree1.Branch("pmtCharge0", pmtCharge0, "pmtCharge0/F")
tree1.Branch("pmtCharge1", pmtCharge1, "pmtCharge1/F")
tree1.Branch("pmtCharge2", pmtCharge2, "pmtCharge2/F")
tree1.Branch("pmtCharge3", pmtCharge3, "pmtCharge3/F")
tree1.Branch("pmtCharge4", pmtCharge4, "pmtCharge4/F")
tree1.Branch("pmtCharge5", pmtCharge5, "pmtCharge5/F")
tree1.Branch("pmtChargePrompt0", pmtChargePrompt0, "pmtChargePrompt0/F")
tree1.Branch("pmtChargePrompt1", pmtChargePrompt1, "pmtChargePrompt1/F")
tree1.Branch("pmtChargePrompt2", pmtChargePrompt2, "pmtChargePrompt2/F")
tree1.Branch("pmtChargePrompt3", pmtChargePrompt3, "pmtChargePrompt3/F")
tree1.Branch("pmtChargePrompt4", pmtChargePrompt4, "pmtChargePrompt4/F")
tree1.Branch("pmtChargePrompt5", pmtChargePrompt5, "pmtChargePrompt5/F")
tree1.Branch("dbCherenkov_portion1", dbCherenkov_portion1, "dbCherenkov_portion1/F") 
tree1.Branch("dbCherenkov_portion2", dbCherenkov_portion2, "dbCherenkov_portion2/F")
tree1.Branch("dbCherenkov_portion3", dbCherenkov_portion3, "dbCherenkov_portion3/F")
tree1.Branch("dbCherenkov_portion4", dbCherenkov_portion4, "dbCherenkov_portion4/F")
tree1.Branch("dbCherenkov_portion5", dbCherenkov_portion5, "dbCherenkov_portion5/F")
tree1.Branch("dbCherenkov_dist1", dbCherenkov_dist1, "dbCherenkov_dist1/F")
tree1.Branch("dbCherenkov_dist2", dbCherenkov_dist2, "dbCherenkov_dist2/F")
tree1.Branch("dbCherenkov_dist3", dbCherenkov_dist3, "dbCherenkov_dist3/F")
tree1.Branch("dbCherenkov_dist4", dbCherenkov_dist4, "dbCherenkov_dist4/F")
tree1.Branch("dbCherenkov_dist5", dbCherenkov_dist5, "dbCherenkov_dist5/F")
tree1.Branch("checkTag", checkTagVal, "checkTag/I")

File = TFile(file0);
# data = File.Get("T")
data = File.Get("T_satCorr");#TTree
nentries = data.GetEntries();
update = int(nentries/100);
#pmtInfo = RAT.PMTInfoUtil.GetPMTInfoUtil()

countTrig = 0
for event in range(nentries):
        #if (event+1)%update == 0:
        #       print(event+1, "analyzed...")
        data.GetEntry(event);    DS = data.ds;#(required)
        if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
            TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);
            ### NOTE: ONLY for MC
            ## MC = DS.GetMC()
            try: CalTrigTime = CAL.GetEventTime()
            except: continue
            # Low Level cuts
            if ( CAL.GetCuts().GetCutWord()&0x31f8 ):  continue; # Low level cut, NOTE: use 0x199 if looking at data between May-Aug 2016 processed in 2016
            if ( TS.dtmTrigSrc&0x82                ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
            if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
            if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
            if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
            if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
            if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform

            countTrig += 1

            qpe = CAL.GetQPE()
            if (qpe <= 0):
                continue 
            fprompt = CAL.GetFprompt()

            nscb = CAL.GetNSCBayes()

            ### NOTE: this is for ar39 acceptance study!!!
            #if nscb >200 or nscb<90:
            #    continue

            posx = -9999999; posy = -9999999; posz = -9999999;
            ## tf2 positions
            posx1 = -9999999; posy1 = -9999999; posz1 = -9999999;
            try:
               posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
            except:
               posx = -9999999; posy = -9999999; posz = -9999999;
            try:
                posx1 = EV.timefit2.GetPosition().X();    posy1 = EV.timefit2.GetPosition().Y(); posz1 = EV.timefit2.GetPosition().Z();
            except:
                posx1 = -9999999; posy1 = -9999999; posz1 = -9999999;

            mbx[0] = posx; mby[0] = posy; mbz[0] = posz;
            tf2x[0] = posx1; tf2y[0] = posy1; tf2z[0] = posz1;

            pmtCharge0[0] = -9999999.
            pmtCharge1[0] = -9999999.
            pmtCharge2[0] = -9999999.
            pmtCharge3[0] = -9999999.
            pmtCharge4[0] = -9999999.
            pmtCharge5[0] = -9999999.

            pmtChargePrompt0[0] = -9999999.
            pmtChargePrompt1[0] = -9999999.
            pmtChargePrompt2[0] = -9999999.
            pmtChargePrompt3[0] = -9999999.
            pmtChargePrompt4[0] = -9999999.
            pmtChargePrompt5[0] = -9999999.

            # the distance of the 0th PMT to itself is 0, so only save from the 1st one
            dbCherenkov_dist1[0] = -9999999.
            dbCherenkov_dist2[0] = -9999999.
            dbCherenkov_dist3[0] = -9999999.
            dbCherenkov_dist4[0] = -9999999.
            dbCherenkov_dist5[0] = -9999999.

            # the portion of the 0th PMT to itself is 100%, so only save from the 1st one
            dbCherenkov_portion1[0] = -9999999.
            dbCherenkov_portion2[0] = -9999999.
            dbCherenkov_portion3[0] = -9999999.
            dbCherenkov_portion4[0] = -9999999.
            dbCherenkov_portion5[0] = -9999999.

            ### Fill the tree
            evtID[0] = TS.GetEventID();
            runIDval[0] = DS.GetRunID();
            subrunIDval[0] = DS.GetSubrunID();
            nSCBayesVal[0] = CAL.GetNSCBayes()
            rprompt60BayesVal[0] = EV.GetRprompt60Bayes()
            qpeVal[0] = qpe
            fpromptVal[0] = fprompt
            fmaxpeVal[0] = EV.GetFmaxpe()
            nhitVal[0] = int(ord(CAL.GetLateNhit()))

            pmtCharge = {}
            pmtChargePrompt = {}
            #pmtPos = {}

            nPMT = CAL.GetPMTCount()
            nPMTs1[0] = nPMT
            # Sort PMTs
            ## NOTE: for spe calculation
            pmtInfo = RAT.PMTInfoUtil.GetPMTInfoUtil()
            pmtInfo.SetRunAndSubrun(runIDval[0], 0)

            for ipmt in range(nPMT):
                rpmt = CAL.GetPMT(ipmt)
                pmtid = rpmt.GetID()
                listPMTid[ipmt] = pmtid
                #spe = pmtSPE[pmtid]
                spe = pmtInfo.GetSPECharge(rpmt) ### !!! added 15 Nov 2022
                if spe <= 0: continue ### !!! added 15 Nov 2022

                #if pmtInfo.IsPMTGoodForAnalysis(rpmt) == False: continue # reserved for future
                #NOTE: use pmt.qPE
                pmt_qPE = rpmt.GetQPE() # might use pmt.GetPromptQ() as alternative
                pmt_qNSCB = rpmt.GetNSCBayes()
                pmt_qPrompt = rpmt.GetPromptQ()/spe
                pmtPhi[ipmt] = Phi[pmtid]
                pmtCosTheta[ipmt] = CosTheta[pmtid]
                ###NOTE: these are array saving to TTree
                pmtQPE[ipmt] = pmt_qPE
                pmtQNSCB[ipmt] = pmt_qNSCB
                pmtQPrompt[ipmt] = pmt_qPrompt
                pmtposTemp = TVector3()
                ###NOTE: these are dictionary for DCC sorting
                pmtCharge[pmtid] = pmt_qPE ### !!! use pmt.qPE for sorting
                pmtChargePrompt[pmtid] = pmt_qPrompt 
                pmtPos[pmtid].SetMag(rAV)

            ## Sort the PMT ids based on their charge in an decreasing order; i.e the 0th item is the brightest PMT
            pmt_charge_sorted = sorted(pmtCharge.items(), key=lambda x:x[1], reverse=True)
            pmt_chargePrompt_sorted = sorted(pmtChargePrompt.items(), key=lambda x:x[1], reverse=True)
            pmtidSorted, chargeSorted, pmtidPromptSorted, chargePromptSorted = [], [], [], []
            for item in pmt_charge_sorted:
                pmtidSorted.append(item[0]) # saves the pmt ids after sorting
                chargeSorted.append(item[1]) # saves the pmt charges after sorting

            for item in pmt_chargePrompt_sorted:
                pmtidPromptSorted.append(item[0]) # saves the pmt ids after sorting
                chargePromptSorted.append(item[1])

            if len(pmtidSorted)<2: ## less than 2 PMTs, abort this event
               continue 

            ## for pmt.qPrompt, not used
            #if len(pmtidPromptSorted)<2: ## less than 2 PMTs, abort this event
            #     continue

            ################################################################################################################################################
            # Step3 : take the charge ratio (portion) and the distance of the PMTs to the brightest PMT in the event
            maxQ = chargeSorted[0] # the charge of the brightest PMT, or Q0
            if (maxQ <= 0):
               continue 
            pmtid_maxQ = pmtidSorted[0] # PMT ID of the brightest PMT
            pmtPos_maxQ = TVector3(pmtPos[pmtid_maxQ][0], pmtPos[pmtid_maxQ][1], pmtPos[pmtid_maxQ][2]) # position of the brightest PMT, X0

            rangePMT = nHiQpmt # investigate 6 highest-charge PMTs by default
            if len(pmtidSorted)<nHiQpmt: # in case less than 6 (but must >= 2)
                rangePMT = len(pmtidSorted)

            ### sort by chargePrompt
            maxQprompt = chargePromptSorted[0] # the charge of the brightest PMT, or Q0
            #if (maxQprompt <= 0):
            #   continue
            #pmtid_maxQ = pmtidPromptSorted[0] # PMT ID of the brightest PMT
            #pmtPos_maxQ = TVector3(pmtPos[pmtid_maxQ][0], pmtPos[pmtid_maxQ][1], pmtPos[pmtid_maxQ][2]) # position of the brightest PMT, X0

            #rangePMT = nHiQpmt # investigate 6 highest-charge PMTs by default
            #if len(pmtidPromptSorted)<nHiQpmt: # in case less than 6 (but must >= 2)
            #    rangePMT = len(pmtidPromptSorted)

            dbCherenkov_pmtMaxpe = maxQ/qpe # the charge of the brightest pmt/event qPE
            pmtFmaxpeVal[0] = dbCherenkov_pmtMaxpe

            qPrompt = qpe*fprompt
            dbCherenkov_pmtMaxpePrompt = maxQprompt/qPrompt
            pmtFmaxpePromptVal[0] = dbCherenkov_pmtMaxpePrompt

            portionsTemp  = []
            distanceTemp  = []

            # Now loop the sorted PMTs and select 6 (nHiQpmt) of them.
            # calculate the PMTs' portions and distances to the max-charge/brightest PMT.
            for i in range(rangePMT):# i==0 is the brightest PMT it
                pmtid = pmtidSorted[i]
                temp_vector = TVector3(pmtPos[pmtid][0], pmtPos[pmtid][1], pmtPos[pmtid][2])
                pmtPhi6[i] = temp_vector.Phi()
                pmtCosTheta6[i] = temp_vector.CosTheta()

                ### NOTE: here use pmt.qPE!!
                portion = chargeSorted[i]/maxQ
                #portion = chargePromptSorted[i]/maxQprompt

                portionsTemp.append(portion)
                dist = (pmtPos_maxQ - temp_vector).Mag()
                distanceTemp.append(dist)
                if i == 0:
                   pmtCharge0[0] = chargeSorted[i] 
                   pmtChargePrompt0[0] = chargePromptSorted[i]
                if i == 1:
                    pmtCharge1[0] = chargeSorted[i]
                    pmtChargePrompt1[0] = chargePromptSorted[i]
                    dbCherenkov_dist1[0] = dist
                    dbCherenkov_portion1[0] = portion
                if i == 2:
                    pmtCharge2[0] = chargeSorted[i]
                    pmtChargePrompt2[0] = chargePromptSorted[i]
                    dbCherenkov_dist2[0] = dist
                    dbCherenkov_portion2[0] = portion
                if i == 3:
                    pmtCharge3[0] = chargeSorted[i]
                    pmtChargePrompt3[0] = chargePromptSorted[i]
                    dbCherenkov_dist3[0] = dist
                    dbCherenkov_portion3[0] = portion
                if i == 4:
                    pmtCharge4[0] = chargeSorted[i]
                    pmtChargePrompt4[0] = chargePromptSorted[i]
                    dbCherenkov_dist4[0] = dist
                    dbCherenkov_portion4[0] = portion
                if i == 5:
                    pmtCharge5[0] = chargeSorted[i]
                    pmtChargePrompt5[0] = chargePromptSorted[i]
                    dbCherenkov_dist5[0] = dist
                    dbCherenkov_portion5[0] = portion

            is_dbCherenkov_like = 0
            index_2nd = 0
            # Search for the double clusters. Note that the "cluster" here actually means a single PMT.
            # Except for the brightest PMT, loop the other 5 PMTs. Check whether there exists a PMT which has enough portion and distance to the brightest one.
            ## Change in Aug-2023, only apply fpmtFmaxpeCut>0.06 and distance>600 mm
            for i in range(1, rangePMT): # exclude the brightest PMT (i == 0)
                if dbCherenkov_pmtMaxpe < fpmtFmaxpeCut:
                    break
                else:
                    if distanceTemp[i]>fdistanceCut:
                        index_2nd = i
                        is_dbCherenkov_like = 1 # tagged, dbCherenkov_like = 1, done this event
                        break
                    else:
                        continue # go to check next PMT

            ### if index_2nd == 0, only 1 PMT is tagged!
            pmtPhi1st[0] = pmtPos_maxQ.Phi()
            pmtCosTheta1st[0] = pmtPos_maxQ.CosTheta()
            pmtCharge1st[0] = maxQ ##maxQprompt ## use pmt.qPrompt

            pmtID_2nd = pmtidSorted[index_2nd] 
            pmtPhi2nd[0] = pmtPos[pmtID_2nd].Phi()
            pmtCosTheta2nd[0] = pmtPos[pmtID_2nd].CosTheta()
            pmtCharge2nd[0] = chargeSorted[index_2nd] ##chargePromptSorted[index_2nd] ## use pmt.qPrompt
            #print "checking", index_2nd, maxQ, chargeSorted[0], chargeSorted[index_2nd]#pmtPhi2nd[0], pmtPhi6[index_2nd], " ", pmtCosTheta2nd[0], pmtCosTheta6[index_2nd]
 
            # fill Others
            qTopRing = EV.GetChargeTopRing(); qSecondRing = EV.GetChargeSecondRing();
            qBotRing = EV.GetChargeBottomRing(); qSecondBotRing = EV.GetChargeSecondBottomRing(); qThirdBotRing = EV.GetChargeThirdBottomRing();
            ### for liquid level, just save values, for further cuts!!!
            pulseGar[0] = EV.pulseindexfirstgar
            cft2r[0] = (qTopRing + qSecondRing)/qpe
            cfb3r[0] = (qBotRing + qSecondBotRing + qThirdBotRing)/qpe
            ### NOTE: neckVeto is different in MC
            #neckVetoVal[0] = MC.GetVetoPMTCount()
            neckVetoVal[0] = CAL.GetNeckVetoPMTCount()

            # Ensure the PMT with the maximum charge is in GAr or above the liquid level
            #if dbCherenkov_pmtID0 > lastGArPMTID:
            #    is_dbCherenkov_like = 0 # if not, not tagged as DCC

            checkTagVal[0] = is_dbCherenkov_like

            del portionsTemp
            del distanceTemp
            tree1.Fill()
            #if eventID == int(check_eventID):
            #    break

fout.cd()
tree1.Write()
