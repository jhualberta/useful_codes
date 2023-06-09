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
fportionCut = 30
fdistCut = 600
nHiQpmt = 6 
switchOffIso = True

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

def calcuIsotropy( pmtDir ):
    if switchOffIso:
       result = (-9999, -9999)
       return result 

    thetaij = 0.0; p1 = 0.0; p2 = 0.0; p3 = 0.0; p4 = 0.0; tij = 0.0;
    result = (0, 0)
    #result.beta14 = 0;
    #result.thetaij = 0;
    if len(pmtDir) == 0: return result;

    for ipmt in range(len(pmtDir)):
        pmtDir1 = pmtDir[ipmt];
        for ipmt2 in range(ipmt + 1, len(pmtDir)):
            pmtDir2 = pmtDir[ipmt2]
            thetaij = pmtDir[ipmt].Angle(pmtDir[ipmt2])
            tij += thetaij
            cosTheta12 = pmtDir1.Dot( pmtDir2 )
            p1 += cosTheta12; ## px is Legendre Polynomial order x
            cosThetaSquared = cosTheta12 * cosTheta12
            p2 += ( 3.0 * cosThetaSquared - 1.0 ) / 2.0;
            p3 += ( 5.0 * cosThetaSquared - 3.0 ) * cosTheta12 / 2.0;
            p4 += ( 35.0 * cosThetaSquared * cosThetaSquared - 30.0 * cosThetaSquared + 3.0 ) / 8.0;

    p1 = 2.0 * p1 / ( len(pmtDir) * ( len(pmtDir) - 1 ) );
    p2 = 2.0 * p2 / ( len(pmtDir) * ( len(pmtDir) - 1 ) );
    p3 = 2.0 * p3 / ( len(pmtDir) * ( len(pmtDir) - 1 ) );
    p4 = 2.0 * p4 / ( len(pmtDir) * ( len(pmtDir) - 1 ) );
    tij = 2.0 * tij /( len(pmtDir) * ( len(pmtDir) - 1 ) );
    beta14_temp = p1 + 4 * p4
    result = (beta14_temp, tij)
    #result.thetaij = tij;
    return result

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
#check_eventID = str(sys.argv[2])

fileName = os.path.basename(file0)
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
Dir=sys.argv[1]

fout = TFile("outSinaBeta14noPortion_"+fileName,"RECREATE");

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
         #print int(row.key), Phi[int(row.key)], CosTheta[int(row.key)]

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
mcx = array('f',[0])
mcy = array('f',[0])
mcz = array('f',[0])

mcx_eEnergyRMSval = array('f',[0])
mcy_eEnergyRMSval = array('f',[0])
mcz_eEnergyRMSval = array('f',[0])

mcx_eEnergyCentroidVal = array('f',[0])
mcy_eEnergyCentroidVal = array('f',[0])
mcz_eEnergyCentroidVal = array('f',[0])

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
timeRes1 = array('f',Nmax*[0]) 
charge1 = array('f',Nmax*[0]) # save pmt.qPE
chargePrompt1 = array('f',Nmax*[0]) # save pmt.qPrompt
#pmtPortion1 = array('f',Nmax*[0])
pmtFmaxpeVal = array('f',[0])

pmtPhi6 = array('f',nHiQpmt*[0])
pmtCosTheta6 = array('f',nHiQpmt*[0])
pmttime6 = array('f',nHiQpmt*[0]) # SubpeakTime
timeRes6 = array('f',nHiQpmt*[0])
charge6 = array('f',nHiQpmt*[0])
chargePrompt6 = array('f',nHiQpmt*[0])
pmtPortion6 = array('f',(nHiQpmt-1)*[0]) ##NOTE: portion to the P0, so only 5 PMTs
pmtPortion2 = array('f',1*[0])
beta14Val = array('f',[0])
thetaijVal = array('f', [0])
beta14_6Val = array('f',[0]) ## just use 6 pmts
thetaij_6Val = array('f', [0])

beta14_truePosVal = array('f',[0])
thetaij_truePosVal = array('f', [0])
beta14_6_truePosVal = array('f',[0]) ## just use 6 pmts
thetaij_6_truePosVal = array('f', [0])

dist_cherVal = array('f', (nHiQpmt-1)*[0]) ##NOTE: distances to the P0, so only 5 PMTs
dist_cher2Val = array('f', [0])
deltaT_cherVal = array('f', [0])
velocity_cherVal = array('f',[0])

tree1 = TTree("T1","saveEvent")
tree1.Branch("runID",runIDval,"runID/L")
tree1.Branch("subrunID",subrunIDval,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("evtx", evtx, "evtx/F")
tree1.Branch("evty", evty, "evty/F")
tree1.Branch("evtz", evtz, "evtz/F")
tree1.Branch("mcx", mcx, "mcx/F")
tree1.Branch("mcy", mcy, "mcy/F")
tree1.Branch("mcz", mcz, "mcz/F")

tree1.Branch("mcx_eEnergyRMS", mcx_eEnergyRMSval, "mcx_eEnergyRMS")
tree1.Branch("mcy_eEnergyRMS", mcy_eEnergyRMSval, "mcy_eEnergyRMS")
tree1.Branch("mcz_eEnergyRMS", mcz_eEnergyRMSval, "mcz_eEnergyRMS")

tree1.Branch("mcx_eEnergyCentroid", mcx_eEnergyCentroidVal, "mcx_eEnergyCentroid")
tree1.Branch("mcy_eEnergyCentroid", mcy_eEnergyCentroidVal, "mcy_eEnergyCentroid")
tree1.Branch("mcz_eEnergyCentroid", mcz_eEnergyCentroidVal, "mcz_eEnergyCentroid")

tree1.Branch("nPMTs",nPMTs1, "nPMTs/I")
### PMT arrays
tree1.Branch("pmtPhi", pmtPhi1, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtCosTheta1, "pmtCosTheta[nPMTs]/F")
tree1.Branch("pmttime", pmttime1, "pmttime[nPMTs]/F")
tree1.Branch("timeRes", timeRes1, "timeRes[nPMTs]/F")
tree1.Branch("charge", charge1, "charge[nPMTs]/F")
tree1.Branch("chargePrompt", chargePrompt1, "chargePrompt[nPMTs]/F")

## just check the 6 PMTs for each event
#nHiQpmt1 = array('i',[0])
#tree1.Branch("nHiQpmt", nHiQpmt1, "nHiQpmt/I")
tree1.Branch("pmtPhi6", pmtPhi6, "pmtPhi6[6]/F")
tree1.Branch("pmtCosTheta6", pmtCosTheta6, "pmtCosTheta6[6]/F")
tree1.Branch("timeRes6", timeRes6, "timeRes6[6]/F")
tree1.Branch("charge6", charge6, "charge6[6]/F")
tree1.Branch("chargePrompt6", chargePrompt6, "chargePrompt6[6]/F")
tree1.Branch("pmtPortion6", pmtPortion6, "pmtPortion6[5]/F")
tree1.Branch("pmtPortion2", pmtPortion2, "pmtPortion2[1]/F")
tree1.Branch("dist_cher", dist_cherVal, "dist_cher[5]/F")

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
tree1.Branch("beta14", beta14Val, "beta14/F")
tree1.Branch("thetaij", thetaijVal, "thetaij/F")
tree1.Branch("beta14_6", beta14_6Val, "beta14_6/F")
tree1.Branch("thetaij_6", thetaij_6Val, "thetaij_6/F")
tree1.Branch("beta14_truePos", beta14_truePosVal, "beta14_truePos/F")
tree1.Branch("thetaij_truePos", thetaij_truePosVal, "thetaij_truePos/F")
tree1.Branch("beta14_6_truePos", beta14_6_truePosVal,"beta14_6_truePos/F")
tree1.Branch("thetaij_6_truePos", thetaij_6_truePosVal, "thetaij_6_truePos/F") 

#tree1.Branch("dist_cher", dist_cherVal, "dist_cher[nHiQpmt]/F")
tree1.Branch("dist_cher2", dist_cher2Val, "dist_cher2/F")

tree1.Branch("deltaT_cher", deltaT_cherVal, "deltaT_cher/F")
tree1.Branch("velocity_cher", velocity_cherVal, "velocity_cher/F")

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
            chargeTemp = {}
            chargePromptTemp = {}
            nPMT = CAL.GetPMTCount()
            nPMTs1[0] = nPMT
            #nHiQpmt1[0] = 10 ## dummy
            pmtTimeTemp = {} #np.zeros(255);
            timeResTemp = {}
            pmtDir = []
            pmtDir_truePos = []
            ### Event Reconstruction MB likelihood valid
            try:
                posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
            except:
                posx = -9999; posy = -9999; posz = -9999;

            evtx[0] = posx; evty[0] = posy; evtz[0] = posz;

            mc_posEvent = MC.GetMCParticle(0).GetPosition()
            if mc_posEvent.Mag() !=0: ### position at center simulations
                mc_posEvent.SetMag(851)

            mcx[0] = mc_posEvent.X(); mcy[0] = mc_posEvent.Y(); mcz[0] = mc_posEvent.Z();

            electronEnergyRMS = MC.GetMCSummary().GetElectronEnergyRMS();
            electronEnergyCentroid = MC.GetMCSummary().GetElectronEnergyCentroid();

            mcx_eEnergyRMSval[0] = electronEnergyRMS.X()
            mcy_eEnergyRMSval[0] = electronEnergyRMS.Y()
            mcz_eEnergyRMSval[0] = electronEnergyRMS.Z()

            mcx_eEnergyCentroidVal[0] = electronEnergyCentroid.X()
            mcy_eEnergyCentroidVal[0] = electronEnergyCentroid.Y()
            mcz_eEnergyCentroidVal[0] = electronEnergyCentroid.Z()

            posEvent =  TVector3(posx, posy, posz)
            for iPMT in range(nPMT):
            	pmt = CAL.GetPMT(iPMT)
            	pmt_id = pmt.GetID()
            	pmt_qPE = pmt.qPE
                # or qPrompt
                pmt_qprompt = pmt.qPrompt
                pmtpos = pmtPos[pmt_id]

                ###NOTE: remove 0 charge!
                ##if pmt_qPE == 0:
                ##   continue
                ##   #print "what's up?", pmt_id 

                pmt_charge[pmt_id] = pmt_qPE

                pmtPhi1[iPMT] = Phi[pmt_id]
                pmtCosTheta1[iPMT] = CosTheta[pmt_id]
                #pmttime1[pmt_id] = pmtTime[pmt_id]

                charge1[iPMT] = pmt_qPE ## for the charges of 6 pmt

                chargePrompt1[iPMT] = pmt_qprompt

                chargeTemp[pmt_id] = pmt_qPE
                chargePromptTemp[pmt_id] = pmt_qprompt
                saveSubpeakTime = []
                #timeRes1[iPMT] = -9999

                pmtTimeTemp[pmt_id] = 100000
                for j in range (pmt.GetPulseCount()): ## loop PMTs
                    PULSE = pmt.GetPulse(j)
                    for ipeak in range(PULSE.GetSubpeakCount()): ## loop subpeaks for a PMT
                        time = PULSE.GetSubpeakTime(ipeak) ##
                        if time < pmtTimeTemp[pmt_id]:
                            saveSubpeakTime.append(time)
                if len(saveSubpeakTime)>0:
                    pmtTimeTemp[pmt_id] = min(saveSubpeakTime)
                    pmttime1[iPMT] = pmtTimeTemp[pmt_id]
                    timeRes1[iPMT] = pmtTimeTemp[pmt_id] - CalTrigTime
                else: 
                    pmtTimeTemp[pmt_id] = 0
                    pmttime1[iPMT] = 0 
                    timeRes1[iPMT] = 0

                if pmtpos.X() != 0 and pmtpos.Y() != 0 and pmtpos.Z() != 0 and posEvent.Mag() !=0:
                    ## NOTE: prompt time cut here!
                    if timeRes1[iPMT]>-28 and timeRes1[iPMT]<60:
                        pmtDir.append( (pmtpos - posEvent).Unit() )
                        pmtDir_truePos.append( (pmtpos - mc_posEvent).Unit() )


            ### use MB position
            results = calcuIsotropy(pmtDir)
            beta14Val[0] = results[0]
            thetaijVal[0] = results[1]

            results_truePos = calcuIsotropy(pmtDir_truePos)
            beta14_truePosVal[0] = results_truePos[0]
            thetaij_truePosVal[0] = results_truePos[1]
            ################################################################################################################################################
            # Step 2: Sort the PMT ids based on their charge in an decreasing order; i.e first item is the brightest PMT
            pmt_charge_sorted = sorted(pmt_charge.items(), key=lambda x:x[1], reverse=True)
            pmtidSorted, chargeSorted = [], []
            for pmtdata in pmt_charge_sorted:
                ## print "!!! pmt data", pmtdata[0], pmtdata[1]
                pmtidSorted.append(pmtdata[0])
                chargeSorted.append(pmtdata[1])
            ## print "!!!!! sorted PMT", pmt_charge_sorted

            if len(pmtidSorted)<2: ## less than 2 PMTs, abort this event-->under-estimate total events!!
               print "too few PMTs, skip this event!"
               continue

            ################################################################################################################################################
            # Step3 : take the charge ratio and the distance of the 5 next brightest PMTs, to THE brightest PMT in the event
            
            # NOTE: **pmtPos** is a list of offline PMT vector positions where the zeroth element corresponds to PMTID=0
            pmtid_maxQ = pmtidSorted[0]
            maxQ = chargeSorted[0]
            pmtPos_maxq = TVector3(pmtPos[pmtid_maxQ][0], pmtPos[pmtid_maxQ][1], pmtPos[pmtid_maxQ][2])
            pmtPos_maxq.SetMag(851) ## Project PMT position to AV sphere!!!
            
            distances = []
            rangePMT = nHiQpmt 
            if len(pmtidSorted)<nHiQpmt:
                rangePMT = len(pmtidSorted)
            ## pmtidSorted saves the pmt ids after sorting
           
            portionsTemp = []
            distanceTemp = []
            pmtDir6= []
            pmtDir6_truePos = []

            ## for the 6 pmts after charge sorting
            
            for i in range(rangePMT):#i==0 is the brightest PMT itself.
                pmtid = pmtidSorted[i]
                temp_vector = TVector3(pmtPos[pmtid][0], pmtPos[pmtid][1], pmtPos[pmtid][2])
                temp_vector.SetMag(851)

                # print( pmt_id[i], int(chargeSorted[i]*100./chargeSorted[0]), int(chargeSorted[i]*100./qPE), round((pmtPos_maxq-temp_vector).Mag(),2))

                pmtPhi6[i] = temp_vector.Phi()
                pmtCosTheta6[i] = temp_vector.CosTheta()
                timeRes6[i] = pmtTimeTemp[pmtid] - CalTrigTime
                if pmtTimeTemp[pmtid] == 0:
                    timeRes6[i] = 0
                charge6[i] = chargeTemp[pmtid]
                chargePrompt6[i] = chargePromptTemp[pmtid] 
                pmtpos = temp_vector
                posEvent =  TVector3(posx, posy, posz)
                pmtDir6.append( (pmtpos - posEvent ).Unit() )
                pmtDir6_truePos.append( (pmtpos - mc_posEvent ).Unit() )

                # print( pmt_id[i], int(charge[i]*100./charge[0]), int(charge[i]*100./qPE), round(Distance(pmtPos_maxq, temp_vector),2))
                #distances.append(Distance(pmtPos_maxq, temp_vector))
                portion = chargeSorted[i]*100./maxQ
                portionsTemp.append(portion) ### or sort(chargeTemp)
                dist0 = (pmtPos_maxq - temp_vector).Mag()
                distanceTemp.append(dist0)
                if i>1:
                   pmtPortion6[i-1] = portion
                   dist_cherVal[i-1] = dist0
                   distances.append(dist0)#Distance(pmtPos_maxq, temp_vector))
                if i==1:
                   pmtPortion2[i-1] = portion
                   dist_cherVal[i-1] = dist0 
                   dist_cher2Val[0] = dist0## only save the 2nd Q PMT
                   distances.append(dist0)#Distance(pmtPos_maxq, temp_vector))
                #pmtPhi1[i] = Phi[pmtid]
                #pmtCosTheta1[i] = CosTheta[pmtid]
                #pmttime1[i] = pmtTime[pmtid]
                #charge1[i] = charge[i] ## for the charges of 6 pmt

                # tagCluster1[kk] = 1
                # print "pmtid", pmtid, kk, pmtPhi1[kk], pmtCosTheta1[kk], pmttime1[kk], charge1[kk], pmttime_unicor1[kk]
                ## print "!!!!! tree1.Fill count", countTrig_afterCuts       
         
            results6 = calcuIsotropy(pmtDir6)
            beta14_6Val[0] = results6[0]
            thetaij_6Val[0] = results6[1]

            #print "portions", portionsTemp
            results6_truePos = calcuIsotropy(pmtDir6_truePos)
            beta14_6_truePosVal[0] = results6_truePos[0]
            thetaij_6_truePosVal[0] = results6_truePos[1]

            ################################################################################################################################################
            # Step 4: apply the cut
            is_dbCherenkov_like = False

            record1stPMTid = -9999
            record2ndPMTid = -9999 # for 2nd highest charge PMT
            dist_temp = -9999
            dbCherenkov_pmtMaxpe = round(maxQ*100./qpe)
            #print "distance", len(distances), "pmtid", len(pmt_id)
            for jj in range(1, rangePMT):
                # if (portions[jj]>30 and distances[jj]>700) or (portions[jj]>60 and distances[jj]>900):
                # if (portions[jj]>35 and distances[jj]>700):# or (portions[jj]>25 and distances[jj]>900):
                # **qPE** is the qPE of the event
                ## if (portions[jj]>30 and distances[jj]>600) and (round(charge[0]*100./qPE)>7):
                if dbCherenkov_pmtMaxpe<7:
                    break
                if portionsTemp[jj]<fportionCut:
                    continue 
                else:
                  if distanceTemp[jj]>fdistCut: ## for all PMT distances
                      is_dbCherenkov_like = True
                      record1stPMTid = pmtidSorted[0]
                      record2ndPMTid = pmtidSorted[jj]
                      dist_temp = distances[jj-1]
                      break
                  else:
                      continue

            if is_dbCherenkov_like == True:
                #dist_cherVal[0] = dist_temp
                tres0 = pmtTimeTemp[record1stPMTid] - CalTrigTime 
                tres1 = pmtTimeTemp[record2ndPMTid] - CalTrigTime
                ## search for Cherenkov in [0, 10] ns
                t_lo = -28; t_hi = 60 
                if tres0>t_lo and tres0<t_hi and tres1>t_lo and tres1<t_hi: 
                    deltaT_cherVal[0] = tres1 - tres0
                    velocity_cherVal[0] =  dist_cher2Val[0]/deltaT_cherVal[0]  #dist_cherVal[0]/deltaT_cherVal[0]
            #else:
            #    #dist_cherVal[0] = -9999
            #    deltaT_cherVal[0] = 0
            #    velocity_cherVal[0] = 0

            ################################################################################################################################################
            check = 0
            if is_dbCherenkov_like == True:
                eventIDToRemove.append(eventID)
                check = 1

            ### Fill the tree
            runIDval[0] = DS.GetRunID(); 
            subrunIDval[0] = DS.GetSubrunID();
            evtID[0] = eventID

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
            pmtFmaxpeVal[0] = chargeSorted[0]*100./qpe

            #print eventID, checkROI, nSCBayesVal[0], rprompt60BayesVal[0], qpeVal[0], fpromptVal[0]
            tree1.Fill()
            #if eventID == int(check_eventID):
            #    break

fout.cd()
tree1.Write()
#fout.Close()
#print eventIDToRemove
