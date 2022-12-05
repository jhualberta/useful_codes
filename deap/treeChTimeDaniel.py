## for analyzing vacuum data to check Th-232 source, using Daniel's cuts
## Recalculate the fprompt (also qpe) with different prompt time windows
## NOTE: using the timefit2, a default fmaxpe<0.75 was applied automatically!!!
import os
import couchdb
import sys
import glob
from rat import *
#from rat import pmtInfo 
from ROOT import *
import numpy as np
from numpy import *
import math
from datetime import datetime, time as datetime_time, timedelta
from array import array

global qc

PI = np.pi
radius = 850

##### set new time-windows here
startIntegral = -28.0 ## time to start total integral (ns) relative to calibrated trigger time
endIntegral = 10000.0 #16000.0 ## time to end integral (ns) relative to calibrated trigger time

startPromptIntegral = -28.0 ## time to start prompt integral (ns) relative to calibrated trigger time
endPromptIntegral = 150.0 ## time to end prompt integral (ns) relative to calibrated trigger time
endPrompt60nsIntegral = 60.0 ## time to end optimized 60 ns prompt integral (ns) relative to calibrated trigger time
endPrompt75nsIntegral = 75.0


alphaStartIntegral = -28.0 ## time to start f_alpha integral relative to calibrated trigger time
alphaEndPromptIntegral = 900.0 ## time to end f_alpha prompt integral relative to calibrated trigger time
alphaEndIntegral = 13500.0 ## time to end f_alpha integral relative to calibrated trigger time

start_early_pulses = 0.0 ## time to start counting early pulses (ns) relative to start of waveform
end_early_pulses = 1600.0 ## time to end counting early pulses (ns) relative to start of waveform
end_early_pulses_extension = 2350.0 ## time to end counting early pulses (ns) relative to start of waveform. start of start_early_pulses_extension automatically defined through (end_early_pulses_extension-end_early_pulses)
## 
frontback_boundary = 8000.0 ## time to split charge into front/back for frontHalfFraction (ns) relative to start of waveform
## 
start400nsIntegral = -28.0 ## time to start 400ns integral (ns)
end400nsIntegral = 400 ## time to end 400ns integral (ns)
## 
start5000nsIntegral = -28.0 ## time to start 5000ns integral (ns)
end5000nsIntegral = 5000 ## time to end 5000ns integral (ns)

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
## pmtX = np.zeros(255)
## pmtY = np.zeros(255)
## pmtZ = np.zeros(255)

pmtInfo = RAT.PMTInfoUtil.GetPMTInfoUtil()

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
                 ## pmtX[ind] = radius*cos(phi)*sin(arccos(cost))
                 ## pmtY[ind] = radius*sin(phi)*sin(arccos(cost))
                 ## pmtZ[ind] = radius*cost

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
fout = TFile("DanielChTest_L7Cuts_"+fileName, "RECREATE")
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

evtxTF2bfVal = array('f',[0])
evtyTF2bfVal = array('f',[0])
evtzTF2bfVal = array('f',[0])

evtxTF2chVal = array('f',[0])
evtyTF2chVal = array('f',[0])
evtzTF2chVal = array('f',[0])
tf2_ch_tVal = array('f',[0])

fmaxpeVal = array('f',[0])

qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])
### fpromptSimple = array('f',[0])
nSCBayesVal = array('f',[0])
rprompt60BayesVal = array('f',[0])
### recalc with different window
qpe_vac = array('f',[0]) #unsigned double
fprompt_vac= array('f',[0])
fprompt60_vac = array('f',[0])
fprompt75_vac = array('f',[0])
nSCBayesRecalc = array('f',[0])
rprompt60BayesRecalc = array('f',[0])

dtmTrigSrcVal = array('f',[0])
subeventNVal = array('f',[0])
eventTimeVal = array('f',[0])
calCutsVal = array('f',[0])
deltatVal  = array('f',[0])
numEarlyPulsesVal = array('f',[0])

#### these values are reserved for master cuts, ChSource STR cuts
pulseindexfirstgarVal = array('f',[0]) ## EV.GetPulseIndexFirstGAr()
cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qPE<0.1

Nmax = 255 ## save how many PMTs for an event, for vacuum, should <200?
NpeakMax = Nmax*50 ## for Nsubpeaks; assume 50 pulses
NsubpeaksMax = Nmax*50*10 ## for pmt time, pulse time; assume 50 pulses and 10 subpeaks for each PMT

#### pmt by pmt information
nPMTs = array('i',[0])
pmtphi = array('f',Nmax*[0])
pmtcosTheta = array('f',Nmax*[0])

pmtposx = array('f',Nmax*[0])
pmtposy = array('f',Nmax*[0])
pmtposz = array('f',Nmax*[0])
nPulse = array('i',Nmax*[0])

#nSubpeaks = array('i',NpeakMax*[0]) ## counts of subpeaks in each pmt
pmttime = array('f',Nmax*[0]) # SubpeakTime
## pmttime1 = array('f',NsubpeaksMax*[0]) # SubpeakTime
## pmttime2 = array('f',NsubpeaksMax*[0]) # SubpeakTime
pulsetime = array('f',Nmax*[0])

##pmtq = array('f',Nmax*[0])

tree1.Branch("runID",runID,"runID/L")
tree1.Branch("subrunID",subrunID,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("fmaxpe", fmaxpeVal, "fmaxpe/F")
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
### needs to test, not correct
### tree1.Branch("fpromptSimple", fpromptSimple, "fpromptSimple/F") ## use pmt.qpe/pmt.q
tree1.Branch("nSCBayes", nSCBayesVal, "nSCBayes/F") # GetNSCBayes()
tree1.Branch("rprompt60Bayes", rprompt60BayesVal, "rprompt60Bayes/F")
######
tree1.Branch("qpe_vac", qpe_vac, "qpe_vac/F")
tree1.Branch("fprompt_vac", fprompt_vac, "fprompt_vac/F")
tree1.Branch("fprompt60_vac", fprompt60_vac, "fprompt60_vac/F")
tree1.Branch("fprompt75_vac", fprompt60_vac, "fprompt75_vac/F")
tree1.Branch("nSCBayes_vac", nSCBayesRecalc, "nSCBayes_vac/F") # GetNSCBayes()
tree1.Branch("rprompt60Bayes_vac", rprompt60BayesRecalc, "rprompt60Bayes_vac/F")

#####
tree1.Branch("subeventN",subeventNVal, "subeventN/F")
tree1.Branch("eventTime",eventTimeVal, "eventTime/F")
tree1.Branch("dtmTrigSrc", dtmTrigSrcVal, "dtmTrigSrc/F")
tree1.Branch("calCuts",calCutsVal, "calCuts/F")
tree1.Branch("deltat",deltatVal, "deltat/F")
tree1.Branch("numEarlyPulses",numEarlyPulsesVal, "numEarlyPulses/F")

tree1.Branch("pulseindexfirstgar", pulseindexfirstgarVal, "pulseindexfirstgar/F")
tree1.Branch("cft2r", cft2r,"cft2r/F")
tree1.Branch("cfb3r", cfb3r,"cfb3r/F")

tree1.Branch("evtx", evtx, "evtx/F") ## MBLikelihood positions
tree1.Branch("evty", evty, "evty/F") ## MBLikelihood positions
tree1.Branch("evtz", evtz, "evtz/F") ## MBLikelihood positions

tree1.Branch("evtxTF2bf", evtxTF2bfVal, "evtxTF2bf/F") ## timefit2 bf positions
tree1.Branch("evtyTF2bf", evtyTF2bfVal, "evtyTF2bf/F") ## timefit2 bf positions
tree1.Branch("evtzTF2bf", evtzTF2bfVal, "evtzTF2bf/F") ## timefit2 bf positions

tree1.Branch("evtxTF2ch", evtxTF2chVal, "evtxTF2ch/F") ## timefit2 ch positions
tree1.Branch("evtyTF2ch", evtyTF2chVal, "evtyTF2ch/F") ## timefit2 ch positions
tree1.Branch("evtzTF2ch", evtzTF2chVal, "evtzTF2ch/F") ## timefit2 ch positions
tree1.Branch("tf2_ch_time", tf2_ch_tVal, "tf2_ch_time/F") ## timefit2 ch time

tree1.Branch("nPMTs",nPMTs, "nPMTs/I")
#tree1.Branch("pmtq", pmtq, "pmtq[nPMTs]/F")

tree1.Branch("pmtPhi", pmtphi, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtcosTheta, "pmtCosTheta[nPMTs]/F")
#tree1.Branch("pmtPosX", pmtposx, "pmtPosX[nPMTs]/F")
#tree1.Branch("pmtPosY", pmtposy, "pmtPosY[nPMTs]/F")
#tree1.Branch("pmtPosZ", pmtposz, "pmtPosZ[nPMTs]/F")
#tree1.Branch("nPulse",nPulse, "nPulse[nPMT]/I") ## each pmt has nPulse

### nSubpeaks = nPMTs*nPulse
#tree1.Branch("nSubpeaks",nSubpeaks, "nSubpeaks[nPulse]/I") ## each pmt has nSubpeaks
### n subpeak time  = nPMTs*nPulse*nSubpeaks
#tree1.Branch("pmttime", pmttime, "pmttime[nPMTs]/F")
#tree1.Branch("pmttime1", pmttime1, "pmttime1[nSubpeaks]/F")
#tree1.Branch("pmttime2", pmttime2, "pmttime2[nSubpeaks]/F")
#tree1.Branch("pulsetime", pulsetime, "pulsetime[nPMTs]/F")

#tree1.Branch("pmttime_uni", pmttime_uni1, "pmttime_uni[nPMTs]/F")
#tree1.Branch("pmttime_unicor", pmttime_unicor1, "pmttime_unicor[nPMTs]/F")

HsubpeakTime = TH1F("HsubpeakTime","subpeak time, cut level 7", 20000, 0, 20000)
HpulseTime = TH1F("HpulseTime","pulse time, cut level 7", 20000, -4000, 16000)
HpulseTimeTF2 = TH1F("HpulseTimeTF2","pulse time, TF2, cut level 7", 20000, -4000, 16000)

HsubpeakTime_Qweight = TH1F("HsubpeakTime_Qweight","subpeak time, cut level 7", 20000, 0, 20000)
HpulseTime_Qweight = TH1F("HpulseTime_Qweight","pulse time, cut level 7", 20000, -4000, 16000)
HpulseTimeTF2_Qweight = TH1F("HpulseTimeTF2_Qweight","pulse time, TF2, q-weighted, cut level 7", 20000, -4000, 16000)

for event in range (nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."

    data.GetEntry(event); DS = data.ds;#(required)
    if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
        TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);

        ### test fitter time valid, cut on eventTime
        try:CalTrigTime = CAL.GetEventTime() #CalibratedTriggerTime()
        except:continue

        ### have to put this to reduce data !! Got to level 7 cuts from Danial and neckVeto
        if ( TS.dtmTrigSrc&0x82 ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
        if ( CAL.GetCuts().GetCutWord()&0x31f8 ): continue;
        if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
        if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
        if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
        if ( CalTrigTime <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
        if ( CalTrigTime >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
        if ( CAL.GetNeckVetoPMTCount() > 0 ): continue; ## neckVetoN == 1, STR cuts from ThSource
        if ( CAL.GetQPE() <= 60               ):  continue; # Skip events that have no charge information (primarily pre-scaled events)

        eventID = TS.GetEventID();
        ### ---- check for 1 event ----
        ## if eventID != 4970: continue
        ## print eventID

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
        eventTimeVal[0] = CalTrigTime # daniel 6th cut: eventTime > 2250 and eventTime < 2700
        qpeVal[0] = CAL.GetQPE() # daniel 7th cut: qPE>60
        fmaxpeVal[0] = EV.GetFmaxpe()

        ### timefit2 valid (fmaxpe<0.75 automatically applied)
        try: 
            time_tf2_ch = EV.timefit2.ch_t
            ### tf2 bf_pos are not valid in vacuum runs
            ##evtxTF2bf = EV.timefit2.bf_pos.fX; evtyTF2bf = EV.timefit2.bf_pos.fY; evtzTF2bf = EV.timefit2.bf_pos.fZ;
            ##evtxTF2ch = EV.timefit2.ch_pos.fX; evtyTF2ch = EV.timefit2.ch_pos.fY; evtzTF2ch = EV.timefit2.ch_pos.fZ;
        except:
            ## print "timefit2 not valid"
            continue

        ### MB likelihood valid
        try:
            posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
        except:
            continue

        ## print "yes !!! everything is ok"

        ## Fred Schuckman cuts
        qTopRing = EV.GetChargeTopRing(); qSecondRing = EV.GetChargeSecondRing();
        qBotRing = EV.GetChargeBottomRing(); qSecondBotRing = EV.GetChargeSecondBottomRing(); qThirdBotRing = EV.GetChargeThirdBottomRing(); 
        ### for liquid level, not used !!!
        pulseindexfirstgarVal[0] =EV.pulseindexfirstgar
        cft2r[0] = (qTopRing + qSecondRing)/qpeVal[0] 
        cfb3r[0] = (qBotRing + qSecondBotRing + qThirdBotRing)/qpeVal[0]

        evtx[0] = posx; evty[0] = posy; evtz[0] = posz;

        tf2_ch_tVal[0] = time_tf2_ch
        ## evtxTF2bfVal[0] = evtxTF2bf; evtyTF2bfVal[0] = evtyTF2bf; evtzTF2bfVal[0] = evtzTF2bf;
        ## evtxTF2chVal[0] = evtxTF2ch; evtyTF2chVal[0] = evtyTF2ch; evtzTF2chVal[0] = evtzTF2ch;

        charges = []; pmtid = []; 
        pmtTime = []; ### time of each PMT
        pmtTime1 = []; pmtTime2 = []; 
        pmtPhi = []; pmtCosTheta = []; pmtNpulse = []; pmtNsubpeaks = []
        ### pmtPosX = []; pmtPosY = []; pmtPosZ = [];

        ## print "!!!! pmtCount", CAL.GetPMTCount()
        countPMTzeroQ = 0
        pmtInfo.SetRunAndSubrun(runID[0], subrunID[0])

        totalChargeOverAllPMTsFullWaveform = 0.
        totalChargeOverAllPMTs = 0
        totalPEOverAllPMTs = 0
        totalPECombinedFullOverAllPMTs = 0

        partCharge_1_PE_AllPMTs = 0
        totalCharge_2_PE_AllPMTs = 0
        promptCharge_1_PE_AllPMTs = 0

        nSubpeaksOverAllPMTsFullWaveform = 0
        nSubpeaksOverAllPMTsFullWaveform_overthresh = 0

        qAPPEAllPMT = 0.;

        totalPromptPEOverAllPMTs = 0;
        totalAlphaPEOverAllPMTs = 0;
        totalAlphaPromptPEOverAllPMTs = 0;
        totalPyrenePEOverAllPMTs = 0;
        totalPyrenePromptPEOverAllPMTs = 0;

        total60nsPEOverAllPMTs = 0;
        total75nsPEOverAllPMTs = 0;
        total400nsPEOverAllPMTs = 0;
        total5000nsPEOverAllPMTs = 0;

        frontHalfChargeOverAllPMTs = 0;
        backHalfChargeOverAllPMTs = 0;

        earlyPulses = 0;
        earlyPulsesExtension = 0;

        ## Figure out if fittime is ok ... only make TOF correction
        ## if Timefit totally succeeded and
        ## multievent disabled (SubeventCount==0)
        ## OR
        ## multievent enabled and SubeventCount==1
        fittimeok = False
        forceCTT = False 
        if CAL.GetSubeventCount()<=1 and CAL.ExistTimeFit() and CAL.ExistCuts():
            if CAL.GetCuts().GetTimeFitStatus()==0 and not(forceCTT):
                 fittimeok = True
            else:
                 CAL.SetForceCTT(forceCTT)
          
        else: 
            CAL.SetForceCTT(forceCTT)

        nhit = 0;
        nhit_5000 = 0;
        promptNhit = 0;
        lateNhit = 0;

        firstLight = 1e6;
        lastLight = 0;
        pmtTime = np.zeros(255)
        ## pulseTime = np.zeros(255) ## PMT * subPeak

        fpromptRatio = 0 ## calculate by using [Sum_pmt (pmt.qpe)/(pmt.qprompt)]/nPMT
        countGoodPMT = 0
        sumQPEpmt = 0
        sumQpmt = 0
        # print "eventID=", eventID, "nPMTs=", CAL.GetPMTCount()
        for ipmt in range (CAL.GetPMTCount()): ### loop PMTs

             PMT = CAL.GetPMT(ipmt)
             PMTid = PMT.id
             nSubpeaksInWindow = 0

             pmtq = PMT.qPrompt
             pmtq1 = PMT.qTotal

             if pmtInfo.IsPMTGoodForAnalysis(PMT) == False: 
                 print "this pmt is bad", ipmt, "id=",PMTid
                 continue
             # if pmtInfo.IsPMTGood(PMTid) == False: continue

             qpePMT = PMT.qPE
             q = PMT.qPrompt #haojie

             ### Get the spe constant for this PMT only once!
             spe = pmtInfo.GetSPECharge(PMT)
             if spe <= 0: continue
             ### just keep PMTs with charge>0 ##???
             #if q == 0 or q<0:
             #    countPMTzeroQ += 1
             #    continue

             countGoodPMT += 1
             ## fpromptRatio += q/qpePMT
             sumQpmt += q
             sumQPEpmt += qpePMT

             charges.append(q)
             pmtPhi.append(Phi[PMTid])
             pmtCosTheta.append(CosTheta[PMTid])

             totalCharge = 0
             totalChargefull = 0
             total60nsCharge = 0
             total400nsCharge = 0
             total5000nsCharge = 0
             promptCharge = 0
             prompt60nsCharge = 0
             prompt75nsCharge = 0
             alphaCharge = 0
             alphaPromptCharge = 0
             ## pyreneCharge = 0
             ## pyrenePromptCharge = 0
             ## qAP = 0.

             partCharge_1 = 0
             totalCharge_2 = 0
             promptCharge_1 = 0

             pmt_has_prompt_pulse = False
             pmt_has_late_pulse = False
             pmt_has_pulse = False
             pmt_has_5000ns_pulse = False
             pmt_clipped = False

             #pmtid.append(PMTid)
             ### print "event", event, q, Phi[PMTid], PMTid
             timeThresh = 100000 # 100 us
             ### print "!!!! pulseCount", PMT.GetPulseCount()
             pmtTime[PMTid] = 100000
             pmt_clipped = False;

             ### LOOP subPulses

             ## pmt_info = rat.utility().GetPMTInfo()
             npulse = PMT.GetPulseCount()
             pmtNpulse.append(npulse)
             ## print eventID, "is there anything?? ", PMTid, "npulse1", npulse, 
             for ipulse in range (npulse):
                  PULSE = PMT.GetPulse(ipulse)
                  nsubpulse = PULSE.GetSubpeakCount()
                  pmtNsubpeaks.append(npulse*nsubpulse)
                  ## print eventID, "is there anything?? npulse", npulse
                  for ipeak in range(nsubpulse): ### loop subpeak
                       ### Check if PMT is saturated.
                       if PULSE.GetSubpeakMinimum(ipeak) == 0: pmt_clipped = True
                       time_subpeak = PULSE.GetSubpeakTime(ipeak) ##
                       subQ = PULSE.GetSubpeakCharge(ipeak)
                       #time1 = PULSE.GetSubpeakTimeUniform(ipeak)
                       #time2 = PULSE.GetSubpeakTimeUniformWithTOFCorrection(ipeak)

                       ### used by Jie: timeThresh put on subpeak time
                       if not(0<time_subpeak and time_subpeak<timeThresh):
                           ## print "subpeak time is outside threshold"
                           continue

                       pulseTime = time_subpeak - CalTrigTime ###
                       pulseTimeTF2 = time_subpeak - time_tf2_ch
                       if fittimeok and not(PMT.IsClipped()):
                           pulseTime = PULSE.GetSubpeakTimeTOFCorrected(ipeak) - CalTrigTime
                           pulseTimeTF2 = PULSE.GetSubpeakTimeTOFCorrected(ipeak) - time_tf2_ch

                           # pmtTime[PMTid] = time
                           # pulse_t = time - CalTrigTime ### 
                           # pulseTime[PMTid] = pulse_t
                           # if pulse_t == 0:
                           #     pulseTime[PMTid] = -9999
                           ## print "event", eventID, "pulse ", ipulse, "subpulse", ipeak, "pmtTime", time, "calTrigTime", CalTrigTime, "pulseTime", pulseTime 
                           ## pulsetime[0] = pulseTime
                       pmtNsubpeaks.append(nsubpulse)
                       
                       HsubpeakTime.Fill(time_subpeak)
                       HpulseTime.Fill(pulseTime)
                       HpulseTimeTF2.Fill(pulseTimeTF2)

                       HsubpeakTime_Qweight.Fill(time_subpeak, pmtq)
                       HpulseTime_Qweight.Fill(pulseTime, pmtq)
                       HpulseTimeTF2_Qweight.Fill(pulseTimeTF2, pmtq)

                       totalChargeOverAllPMTsFullWaveform += subQ/spe
                       nSubpeaksOverAllPMTsFullWaveform += 1
                       if subQ/spe > 0.8:
                          nSubpeaksOverAllPMTsFullWaveform_overthresh += 1

                       totalChargefull += subQ 
                       ### Early light is with respect to waveform not event time.
                       ### nominal large window

                       ## pulse_t1 = time_subpeak - CalTrigTime; ### PULSE.GetTime vs. PULSE.GetSubpeakTime, using time_subpeak gets less 0 peaks, but mostly are similar
                       # print "!!!! checking prompt charge", "pmt=", ipmt, "pulse",ipulse, "subpeak", ipeak, " pulseTime, ", pulseTime, "[",startPromptIntegral,endPromptIntegral,"]"
                       if pulseTime > startIntegral and pulseTime < endIntegral:
                           # print startIntegral, "<", pulse_t, "<", endIntegral
                           totalCharge += subQ  ### total charge in -28, 10000 ns window
                           pmt_has_pulse = True
                           nSubpeaksInWindow += 1
                           ## if pulseTime < (firstLight - CalTrigTime):
                           ##      firstLight = pulseTime + CalTrigTime
                           ## if pulseTime > (lastLight - CalTrigTime):
                           ##      lastLight = pulseTime + CalTrigTime;

                       ### calculate promptCharge
                       # print "!!!! checking prompt charge pulseTime", pulseTime, "[",startPromptIntegral,endPromptIntegral,"]"
                       if pulseTimeTF2 > startPromptIntegral and pulseTimeTF2 < endPromptIntegral:
                           ## print startPromptIntegral, "<", pulse_t, "<", endPromptIntegral, " ", endPrompt60nsIntegral
                           promptCharge += subQ
                           if pulseTime < endPrompt75nsIntegral: prompt75nsCharge += subQ
                           if pulseTime < endPrompt60nsIntegral: prompt60nsCharge += subQ

                           pmt_has_prompt_pulse = True
                           ## print "yes in window, promptCharge= ", promptCharge

                       if pulseTimeTF2 > start400nsIntegral and pulseTimeTF2 < end400nsIntegral:
                           total400nsCharge += subQ 

                       if pulseTimeTF2 > start5000nsIntegral and pulseTimeTF2 < end5000nsIntegral: 
                           total5000nsCharge += subQ
                           pmt_has_5000ns_pulse = True

                       if pulseTimeTF2 > endPromptIntegral and pulseTimeTF2 < endIntegral:
                           pmt_has_late_pulse = True

                       if pulseTimeTF2 > alphaStartIntegral and pulseTimeTF2 < alphaEndIntegral:
                           alphaCharge += subQ 
                           if pulseTimeTF2 < alphaEndPromptIntegral:
                                alphaPromptCharge += subQ

             ### end loop of PULSE                   
             ### Compute PE using charge division. Note this allows us
             ### to have a fractional number of PE!

             ### Int_{-28}^{60} Sum_{pmt = 0 to 254} {Q_i(t)/SPE_i} / Int_{-28}^{10000} Sum_{pmt = 0 to 254} {Q_i(t)/SPE_i}
             totalPE = totalCharge / spe;
             promptPE = promptCharge / spe;
             alphaPE = alphaCharge / spe;
             alphaPromptPE = alphaPromptCharge / spe;

             prompt60nsPE = prompt60nsCharge / spe;
             prompt75nsPE = prompt75nsCharge / spe;

             partCharge_1_PE = partCharge_1 / spe;
             totalCharge_2_PE = totalCharge_2 / spe;
             promptCharge_1_PE = promptCharge_1 / spe;

             totalChargeOverAllPMTs += totalCharge; ## total Q, adding all PMTs in this event
             totalPEOverAllPMTs += totalPE; ## total Q/spe, adding all PMTs in this event
             ## qAPPEAllPMT += qAP/spe;
          
             totalPromptPEOverAllPMTs += promptPE;
             totalAlphaPEOverAllPMTs += alphaPE;
             totalAlphaPromptPEOverAllPMTs += alphaPromptPE;
    
             total60nsPEOverAllPMTs += prompt60nsPE; ## for QPE_60ns
             total75nsPEOverAllPMTs += prompt75nsPE;

             total400nsPEOverAllPMTs += total400nsCharge / spe;
             total5000nsPEOverAllPMTs += total5000nsCharge / spe;
    
             partCharge_1_PE_AllPMTs += partCharge_1_PE;
             totalCharge_2_PE_AllPMTs += totalCharge_2_PE;
             promptCharge_1_PE_AllPMTs += promptCharge_1_PE;


        ## end of PMT loop
        if (countGoodPMT>0 and sumQpmt>0):
            fpromptRatio = sumQPEpmt/sumQpmt

        ## fpromptSimple[0] = fpromptRatio 

        ## Compute Fprompt using division of SPE-calibrated charges.
        ## fprompt_vac[0] = 0.
        if totalPEOverAllPMTs > 0.:
            fprompt_vac[0] = totalPromptPEOverAllPMTs / totalPEOverAllPMTs
            fprompt60_vac[0] = total60nsPEOverAllPMTs / totalPEOverAllPMTs 
            fprompt75_vac[0] = total75nsPEOverAllPMTs / totalPEOverAllPMTs

        #if fprompt_vac[0] == 1:
        #   print "look look", fprompt_vac[0], "normal", fpromptVal[0], "60ns", fprompt60_vac[0], " totalPromptPE ", totalPromptPEOverAllPMTs, " totalPEOverAllPMTs ", totalPEOverAllPMTs
 
        ###  SPE = pmtInfo.GetSPECharge(pmt);
        ###  pulse.GetSubpeakCharge(ipeak)/SPE
             
        qpe_vac[0] = totalPEOverAllPMTs
        #### Set event-level charge information in both CAL and EV branches
        ## cal->SetTotalQ(totalChargeOverAllPMTs);
        ## cal->SetQAP(qAPPEAllPMT);
        ## cal->SetQPEfullW(totalChargeOverAllPMTsFullWaveform);
        ## cal->SetQPEfullWCombined(totalPECombinedFullOverAllPMTs);
        ## cal->SetSubpeakFullWCount(nSubpeaksOverAllPMTsFullWaveform);
        ## cal->SetSubpeakFullWOverThresholdCount(nSubpeaksOverAllPMTsFullWaveform_overthresh);
        ## cal->SetQPE_60(total60nsPEOverAllPMTs);
        ## cal->SetQPE_400(total400nsPEOverAllPMTs);
        ## cal->SetQPE_5000(total5000nsPEOverAllPMTs);
        ## cal->SetNumEarlyPulses(earlyPulses);
        ## cal->SetNumEarlyPulsesExtension(earlyPulsesExtension);

        # pmtTime = list(pmtTime)
        # pulseTime = list(pulseTime)

        nPMTs[0] = len(charges)

        ### print "!!! nPMTs", nPMTs[0], "zero Q pmts:", countPMTzeroQ, "nSubpeaks", nSubpeaks[0]
        ## must swtich list to array to save tree

        if len(charges) != len(pmtPhi):
            print "pmt data dimension error, quit the event loop ..."
            continue

        pmtTime_filter = filter(lambda x: x != 0, pmtTime)
        # pulseTime_filter = filter(lambda x: x != 0, pulseTime)

        ##print "event", eventID, "nPMT", nPMTs[0], " pmtTime ", len(pmtTime_filter), " ", len(pulseTime_filter) 

        for i in range(nPMTs[0]):
              #indx = id_filter[i]
              ## pmtq[i] = charges[i]
              pmtphi[i] = pmtPhi[i]
              pmtcosTheta[i] = pmtCosTheta[i]
              pmttime[i] = pmtTime_filter[i]
              #if pulseTime_filter[i] == -9999:
              #    pulsetime[i] = 0
              #else:    
              #    pulsetime[i] = pulseTime_filter[i]
              #pmttime1[i] = pmtTime1[i]
              #pmttime2[i] = pmtTime2[i] 

        ### filter zeros from initialization the listers
        ## put conditions here; dimensions should be the same!!

        ## pmtq = filter(lambda x: x != 0, pmtq)
        ## pmttime = filter(lambda x: x!= 0, pmttime)
        ## pmtpos = zip(pmtphi,pmtcosTheta)
        ## pmtpos = filter(lambda x: x != 0, pmtpos)
        ## pmtphi = [x[0] for x in pmtpos]
        ## pmtcosTheta = [x[1] for x in pmtpos]

        ## if len(pmtq) != len(pmtphi) or len(pmtq) != len(pmtcosTheta):
        ##     print "After filter 0, pmt data dimension error, quit the event loop ..."
        ##     continue

        tree1.Fill()
#################################################################################################################################
#dstreeclone.Write()
#fout.Close()
#File.Close()

fout.cd()
tree1.Write()
HsubpeakTime.Write()
HpulseTime.Write()
HpulseTimeTF2.Write()

HsubpeakTime_Qweight.Write()
HpulseTime_Qweight.Write()
HpulseTimeTF2_Qweight.Write()

fout.Close()
