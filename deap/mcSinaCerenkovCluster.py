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
portionCut = 10
def Distance(V1, V2):
    '''the distance betwwen two PMTs on the surface of the sphere/detector'''
    dist = (V1-V2).Mag();
    return dist

def N_near_pmts(N_pmts, event_position):
    '''Find the N nearest PMT to the given PMT'''
    PMT_Distance = {}

    for i in range(255):
        PMTPosition = TVector3()
        PMTPosition.SetXYZ(pmtPos[i][0], pmtPos[i][1], pmtPos[i][2])
        PMTPosition.SetMag(851) ### Project PMT position to AV sphere!!
        PMT_Distance[i] = Distance(PMTPosition, event_position)

    Sorted_Dic = sorted(PMT_Distance.items(), key = lambda kv:(kv[1], kv[0]))# Sort the dictionary into pairs based on the distance
    ID = []
    for i in range(N_pmts):
        ID.append(Sorted_Dic[i][0])
    return Sorted_Dic[:N_pmts], ID


print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
#check_eventID = str(sys.argv[2])

fileName = os.path.basename(file0)
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
Dir=sys.argv[1]

fout = TFile("outSinaCerenNoPortion_"+fileName,"RECREATE");

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
evtx = array('f',[0])
evty = array('f',[0])
evtz = array('f',[0])
qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])
nSCBayesVal = array('f',[0])
rprompt60BayesVal = array('f',[0])
fmaxpeVal = array('f',[0])
checkTagVal = array('I', [0])
nhitVal = array('i',[0])
pulseindexfirstgarVal = array('f',[0])
cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qPE<0.1
neckVetoVal = array('f',[0])

Nmax = 255
#### pmt by pmt, save cluster patterns
nPMTs1 = array('i',[0])
pmtPhi1 = array('f',Nmax*[0])
pmtCosTheta1 = array('f',Nmax*[0])
pmttime1 = array('f',Nmax*[0]) # SubpeakTime
charge1 = array('f',Nmax*[0])
pmtPortion1 = array('f',Nmax*[0])
pmtFmaxpeVal = array('f',[0])

tree1 = TTree("T1","saveEvent")
tree1.Branch("runID",runIDval,"runID/L")
tree1.Branch("subrunID",subrunIDval,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("evtx", evtx, "evtx/F")
tree1.Branch("evty", evty, "evty/F")
tree1.Branch("evtz", evtz, "evtz/F")
tree1.Branch("nPMTs",nPMTs1, "nPMTs/I")
### PMT arrays
tree1.Branch("pmtPhi", pmtPhi1, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtCosTheta1, "pmtCosTheta[nPMTs]/F")
tree1.Branch("pmttime", pmttime1, "pmttime[nPMTs]/F")
tree1.Branch("charge", charge1, "charge[nPMTs]/F")
tree1.Branch("pmtPortion", pmtPortion1, "pmtPortion[nPMTs]/F")

tree1.Branch("pmtFmaxpe", pmtFmaxpeVal, "pmtFmaxpe/F")
tree1.Branch("nSCBayes", nSCBayesVal, "nSCBayes/F") # GetNSCBayes()
tree1.Branch("rprompt60Bayes", rprompt60BayesVal, "rprompt60Bayes/F")
tree1.Branch("fmaxpe", fmaxpeVal, "fmaxpe/F")
tree1.Branch("checkTag", checkTagVal, "checkTag/I")
tree1.Branch("nhit", nhitVal, "nhit/I")
tree1.Branch("pulseindexfirstgar", pulseindexfirstgarVal, "pulseindexfirstgar/F")
tree1.Branch("cft2r", cft2r,"cft2r/F")
tree1.Branch("cfb3r", cfb3r,"cfb3r/F")
tree1.Branch("neckVeto", neckVetoVal, "neckVeto/F")

File = TFile(file0);
data = File.Get("T")
#data = File.Get("T_satCorr");#TTree
nentries = data.GetEntries();
update = int(nentries/100);

eventIDToRemove = []
countTrig = 0
for event in range(nentries):
        #if (event+1)%update == 0:
        #       print(event+1, "analyzed...")
        data.GetEntry(event);    DS = data.ds;#(required)
        if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
            TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);
            ### NOTE: ONLY for MC
            MC = DS.GetMC()
            try: CalTrigTime = CAL.GetEventTime()
            except: continue
            countTrig += 1
            # Low Level cuts
            #if ( CAL.GetCuts().GetCutWord()&0x31f8 ):  continue; # Low level cut, NOTE: use 0x199 if looking at data between May-Aug 2016 processed in 2016
            #if ( TS.dtmTrigSrc&0x82                ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
            #if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
            #if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
            #if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
            #if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
            #if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
            ### if ( EV.GetFmaxpe() < 0.15             ):  continue; #
            #if ( EV.GetFmaxpe() > 0.4              ):  continue; #
            #if ( CAL.GetFprompt() > 0.55           ):  continue; # was >0.58
            #if ( CAL.GetQPE() <= 60               ):  continue; # Skip events that have no charge information (primarily pre-scaled events), was <200

            ################################################################################################################################################
            #Inside the for loop of your events:
            eventID = TS.GetEventID();
            qpe = CAL.GetQPE() 
            ## just check interested eventID
            ## if eventID != int(check_eventID): continue

            # Step 1: take the charge and the pmt ids
            pmt_charge = {}
            nPMT = CAL.GetPMTCount()
            nPMTs1[0] = nPMT
            pmtTime = {} #np.zeros(255);
            for iPMT in range(nPMT):
                pmt = CAL.GetPMT(iPMT)
                pmt_id = pmt.GetID()
                pmt_qPE = pmt.qPE
                #or qPrompt
                pmt_charge[pmt_id] = pmt_qPE
                saveSubpeakTime = []
                pmtTime[pmt_id] = 100000
                for j in range (pmt.GetPulseCount()): ## loop PMTs
                    PULSE = pmt.GetPulse(j)
                    for ipeak in range(PULSE.GetSubpeakCount()): ## loop subpeaks for a PMT
                        time = PULSE.GetSubpeakTime(ipeak) ##
                        if time < pmtTime[pmt_id]:
                            saveSubpeakTime.append(time)
                if len(saveSubpeakTime)>0:
                    pmtTime[pmt_id] = min(saveSubpeakTime)
                else: pmtTime[pmt_id] = 0

            ################################################################################################################################################
            # Step 2: Sort the PMT ids based on their charge in an decreasing order; i.e first item is the brightest PMT
            pmt_charge_sorted = sorted(pmt_charge.items(), key=lambda x:x[1], reverse=True)
            pmt_id, charge = [], []
            for pmtdata in pmt_charge_sorted:
                ## print "!!! pmt data", pmtdata[0], pmtdata[1]
                pmt_id.append(pmtdata[0])
                charge.append(pmtdata[1])
                ## print "!!!!! sorted PMT", pmt_charge_sorted
            ################################################################################################################################################
            # Step3 : take the charge ratio and the distance of the 5 next brightest PMTs, to THE brightest PMT in the event
            
            # NOTE: **pmtPos** is a list of offline PMT vector positions where the zeroth element corresponds to PMTID=0
            maxcharge_vector = TVector3(pmtPos[pmt_id[0]][0], pmtPos[pmt_id[0]][1], pmtPos[pmt_id[0]][2])
            maxcharge_vector.SetMag(851) ## Project PMT position to AV sphere!!!
            
            distances = []
            portions  = []
            rangePMT = 6
            if len(pmt_id)<6:
                rangePMT = len(pmt_id)

            for i in range(rangePMT):#i==0 is the brightest PMT itself.
                pmtid = pmt_id[i]
                temp_vector = TVector3(pmtPos[pmtid][0], pmtPos[pmtid][1], pmtPos[pmtid][2])
                temp_vector.SetMag(851)
                # print( pmt_id[i], int(charge[i]*100./charge[0]), int(charge[i]*100./qPE), round((maxcharge_vector-temp_vector).Mag(),2))
            
                # print( pmt_id[i], int(charge[i]*100./charge[0]), int(charge[i]*100./qPE), round(Distance(maxcharge_vector, temp_vector),2))
                distances.append(Distance(maxcharge_vector, temp_vector))
                portions.append(charge[i]*100./charge[0])

                pmtPhi1[i] = Phi[pmtid]
                pmtCosTheta1[i] = CosTheta[pmtid]
                pmttime1[i] = pmtTime[pmtid]
                charge1[i] = charge[i]

                pmtPortion1[i] = charge[i]*100./charge[0]
                # tagCluster1[kk] = 1
                # print "pmtid", pmtid, kk, pmtPhi1[kk], pmtCosTheta1[kk], pmttime1[kk], charge1[kk], pmttime_unicor1[kk]
                ## print "!!!!! tree1.Fill count", countTrig_afterCuts       

            ################################################################################################################################################
            # Step 4: apply the cut
            remove_event = False
            rangePMT = 6
            if len(pmt_id)<6:
                rangePMT = len(pmt_id)
            for jj in range(1, rangePMT):
                # if (portions[jj]>30 and distances[jj]>700) or (portions[jj]>60 and distances[jj]>900):
                # if (portions[jj]>35 and distances[jj]>700):# or (portions[jj]>25 and distances[jj]>900):
                # **qPE** is the qPE of the event
                ## if (portions[jj]>30 and distances[jj]>600) and (round(charge[0]*100./qPE)>7):
                if (round(charge[0]*100./qpe)<7):# or portions[jj]<portionCut:# or (round(charge[0]*100./qpe)<7):
                    continue 
                else:
                  if distances[jj]>600:
                      remove_event = True
                      break
            ################################################################################################################################################
            check = 0
            if remove_event == True:
                eventIDToRemove.append(eventID)
                check = 1

            ### Fill the tree
            runIDval[0] = DS.GetRunID(); 
            subrunIDval[0] = DS.GetSubrunID();
            evtID[0] = eventID

            ### Event Reconstruction MB likelihood valid
            try:
                posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
            except:
                posx = -9999; posy = -9999; posz = -9999;

            evtx[0] = posx; evty[0] = posy; evtz[0] = posz;
            
            nSCBayesVal[0] = CAL.GetNSCBayes()
            rprompt60BayesVal[0] = EV.GetRprompt60Bayes()
            qpe = CAL.GetQPE()
            qpeVal[0] = qpe
            fpromptVal[0] = CAL.GetFprompt()
            checkTagVal[0] = check
            fmaxpeVal[0] = EV.GetFmaxpe()
            nhitVal[0] = int(ord(CAL.GetLateNhit()))

            ## Fred Schuckman cuts
            qTopRing = EV.GetChargeTopRing(); qSecondRing = EV.GetChargeSecondRing();
            qBotRing = EV.GetChargeBottomRing(); qSecondBotRing = EV.GetChargeSecondBottomRing(); qThirdBotRing = EV.GetChargeThirdBottomRing();
            ### for liquid level, just save values, for further cuts!!!
            pulseindexfirstgarVal[0] =EV.pulseindexfirstgar
            cft2r[0] = (qTopRing + qSecondRing)/qpe
            cfb3r[0] = (qBotRing + qSecondBotRing + qThirdBotRing)/qpe
            #neckVetoVal[0] = CAL.GetNeckVetoPMTCount() 
            ### NOTE: neckVeto is different in MC
            neckVetoVal[0] = MC.GetVetoPMTCount()
            pmtFmaxpeVal[0] = charge[0]*100./qpe

            #print eventID, checkROI, nSCBayesVal[0], rprompt60BayesVal[0], qpeVal[0], fpromptVal[0]
            tree1.Fill()
            #if eventID == int(check_eventID):
            #    break

fout.cd()
tree1.Write()
#fout.Close()
#print eventIDToRemove
