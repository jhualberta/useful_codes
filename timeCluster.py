#!/usr/bin/python2
# Print the location of each PMT
# @author Tina Pollmann, 2015
# updated 2017

''' This comes with RAT, so source rat.env first! '''
import couchdb
import sys
import glob
from rat import *
#from rat import PMTInfoUtil
from ROOT import *
import numpy as np
import math
from datetime import datetime, time as datetime_time, timedelta
from array import array

QclusterCut = 100 ### sum charge of cluster, was 114
QpmtCut = 2

global qc

PI = np.pi

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
fout = TFile("TimeCluster_"+file0,"RECREATE");
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])
server = couchdb.Server("https://ddaq3.snolab.ca:6984")
db = server["deapdb"]
Dir=sys.argv[1]
AllFiles=glob.glob(Dir)
AllFiles = sorted(AllFiles)

### save as different cluster events
tree1 = TTree("T1","single_cluster")
tree2 = TTree("T2","double_cluster")
tree3 = TTree("T3","triple_cluster")

### save cluster by cluster
tree4 = TTree("T4","cluster1")
tree5 = TTree("T5","cluster2")
tree6 = TTree("T6","cluster3")

### these are for each event, common for all cluster
runID = array('l',[0])
subrunID = array('l',[0])
evtID = array('l',[0])
evtx = array('f',[0])
evty = array('f',[0])
evtz = array('f',[0])
ncluster = array('i',[0])
qpeVal = array('f',[0]) #unsigned double 
fpromptVal = array('f',[0])

Nmax = 200

#### pmt by pmt, cluster 1
nPMTs1 = array('i',[0])
pmtPhi1 = array('f',Nmax*[0])
pmtCosTheta1 = array('f',Nmax*[0])
pmttime1 = array('f',Nmax*[0])
charge1 = array('f',Nmax*[0]) # SubpeakTime
pmttime_uni1 = array('f',Nmax*[0]) #uniform 
pmttime_unicor1 = array('f',Nmax*[0]) # corrected

#### pmt by pmt, cluster 2
nPMTs2 = array('i',[0])
pmtPhi2 = array('f',Nmax*[0])
pmtCosTheta2 = array('f',Nmax*[0])
pmttime2 = array('f',Nmax*[0])
charge2 = array('f',Nmax*[0])
pmttime_uni2 = array('f',Nmax*[0]) #uniform 
pmttime_unicor2 = array('f',Nmax*[0]) # corrected

#### pmt by pmt, cluster 3
nPMTs3 = array('i',[0])
pmtPhi3 = array('f',Nmax*[0])
pmtCosTheta3 = array('f',Nmax*[0])
pmttime3 = array('f',Nmax*[0])
charge3 = array('f',Nmax*[0])
pmttime_uni3 = array('f',Nmax*[0]) #uniform
pmttime_unicor3 = array('f',Nmax*[0]) # corrected

tree1.Branch("runID",runID,"runID/L")
tree1.Branch("subrunID",subrunID,"subrunID/L")
tree1.Branch("eventID",evtID,"eventID/L") #ULong64_
tree1.Branch("qpe", qpeVal, "qpe/F")
tree1.Branch("fprompt", fpromptVal, "fprompt/F")
tree1.Branch("evtx", evtx, "evtx/F")
tree1.Branch("evty", evty, "evty/F")
tree1.Branch("evtz", evtz, "evtz/F")
tree1.Branch("ncluster",ncluster,"ncluster/I")
tree1.Branch("nPMTs",nPMTs1, "nPMTs/I")
tree1.Branch("pmtPhi", pmtPhi1, "pmtPhi[nPMTs]/F")
tree1.Branch("pmtCosTheta", pmtCosTheta1, "pmtCosTheta[nPMTs]/F")
tree1.Branch("pmttime", pmttime1, "pmttime[nPMTs]/F")
tree1.Branch("pmttime_uni", pmttime_uni1, "pmttime_uni[nPMTs]/F") 
tree1.Branch("pmttime_unicor", pmttime_unicor1, "pmttime_unicor[nPMTs]/F")
tree1.Branch("charge", charge1, "charge[nPMTs]/F")

tree2.Branch("runID",runID,"runID/L")
tree2.Branch("subrunID",subrunID,"subrunID/L")
tree2.Branch("eventID",evtID,"eventID/L") #ULong64_
tree2.Branch("qpe", qpeVal, "qpe/F")
tree2.Branch("fprompt", fpromptVal, "fprompt/F")
tree2.Branch("evtx", evtx, "evtx/F")
tree2.Branch("evty", evty, "evty/F")
tree2.Branch("evtz", evtz, "evtz/F")
tree2.Branch("ncluster",ncluster,"ncluster/I")
tree2.Branch("nPMTs",nPMTs2, "nPMTs/I")
tree2.Branch("pmtPhi", pmtPhi2, "pmtPhi[nPMTs]/F")
tree2.Branch("pmtCosTheta", pmtCosTheta2, "pmtCosTheta[nPMTs]/F")
tree2.Branch("pmttime", pmttime2, "pmttime[nPMTs]/F")
tree2.Branch("pmttime_uni", pmttime_uni2, "pmttime_uni[nPMTs]/F")
tree2.Branch("pmttime_unicor", pmttime_unicor2, "pmttime_unicor[nPMTs]/F")
tree2.Branch("charge", charge2, "charge[nPMTs]/F")

tree3.Branch("runID",runID,"runID/L")
tree3.Branch("subrunID",subrunID,"subrunID/L")
tree3.Branch("eventID",evtID,"eventID/L") #ULong64_
tree3.Branch("qpe", qpeVal, "qpe/F")
tree3.Branch("fprompt", fpromptVal, "fprompt/F")
tree3.Branch("evtx", evtx, "evtx/F")
tree3.Branch("evty", evty, "evty/F")
tree3.Branch("evtz", evtz, "evtz/F")
tree3.Branch("ncluster",ncluster,"ncluster/I")
tree3.Branch("nPMTs",nPMTs3, "nPMTs/I")
tree3.Branch("pmtPhi", pmtPhi3, "pmtPhi[nPMTs]/F")
tree3.Branch("pmtCosTheta", pmtCosTheta3, "pmtCosTheta[nPMTs]/F")
tree3.Branch("pmttime", pmttime3, "pmttime[nPMTs]/F")
tree3.Branch("pmttime_uni", pmttime_uni3, "pmttime_uni[nPMTs]/F")
tree3.Branch("pmttime_unicor", pmttime_unicor3, "pmttime_unicor[nPMTs]/F")
tree3.Branch("charge", charge3, "charge[nPMTs]/F")

tree4.Branch("runID",runID,"runID/L")
tree4.Branch("subrunID",subrunID,"subrunID/L")
tree4.Branch("eventID",evtID,"eventID/L") #ULong64_
tree4.Branch("qpe", qpeVal, "qpe/F")
tree4.Branch("fprompt", fpromptVal, "fprompt/F")
tree4.Branch("evtx", evtx, "evtx/F")
tree4.Branch("evty", evty, "evty/F")
tree4.Branch("evtz", evtz, "evtz/F")
tree4.Branch("ncluster",ncluster,"ncluster/I")
tree4.Branch("nPMTs",nPMTs1, "nPMTs/I")
tree4.Branch("pmtPhi", pmtPhi1, "pmtPhi[nPMTs]/F")
tree4.Branch("pmtCosTheta", pmtCosTheta1, "pmtCosTheta[nPMTs]/F")
tree4.Branch("pmttime", pmttime1, "pmttime[nPMTs]/F")
tree4.Branch("pmttime_uni", pmttime_uni1, "pmttime_uni[nPMTs]/F")
tree4.Branch("pmttime_unicor", pmttime_unicor1, "pmttime_unicor[nPMTs]/F")
tree4.Branch("charge", charge1, "charge[nPMTs]/F")

tree5.Branch("runID",runID,"runID/L")
tree5.Branch("subrunID",subrunID,"subrunID/L")
tree5.Branch("eventID",evtID,"eventID/L") #ULong64_
tree5.Branch("qpe", qpeVal, "qpe/F")
tree5.Branch("fprompt", fpromptVal, "fprompt/F")
tree5.Branch("evtx", evtx, "evtx/F")
tree5.Branch("evty", evty, "evty/F")
tree5.Branch("evtz", evtz, "evtz/F")
tree5.Branch("ncluster",ncluster,"ncluster/I")
tree5.Branch("nPMTs",nPMTs2, "nPMTs/I")
tree5.Branch("pmtPhi", pmtPhi2, "pmtPhi[nPMTs]/F")
tree5.Branch("pmtCosTheta", pmtCosTheta2, "pmtCosTheta[nPMTs]/F")
tree5.Branch("pmttime", pmttime2, "pmttime[nPMTs]/F")
tree5.Branch("pmttime_uni", pmttime_uni2, "pmttime_uni[nPMTs]/F")
tree5.Branch("pmttime_unicor", pmttime_unicor2, "pmttime_unicor[nPMTs]/F")
tree5.Branch("charge", charge2, "charge[nPMTs]/F")

tree6.Branch("runID",runID,"runID/L")
tree6.Branch("subrunID",subrunID,"subrunID/L")
tree6.Branch("eventID",evtID,"eventID/L") #ULong64_
tree6.Branch("qpe", qpeVal, "qpe/F")
tree6.Branch("fprompt", fpromptVal, "fprompt/F")
tree6.Branch("evtx", evtx, "evtx/F")
tree6.Branch("evty", evty, "evty/F")
tree6.Branch("evtz", evtz, "evtz/F")
tree6.Branch("ncluster",ncluster,"ncluster/I")
tree6.Branch("nPMTs",nPMTs3, "nPMTs/I")
tree6.Branch("pmtPhi", pmtPhi3, "pmtPhi[nPMTs]/F")
tree6.Branch("pmtCosTheta", pmtCosTheta3, "pmtCosTheta[nPMTs]/F")
tree6.Branch("pmttime", pmttime3, "pmttime[nPMTs]/F")
tree6.Branch("pmttime_uni", pmttime_uni3, "pmttime_uni[nPMTs]/F")
tree6.Branch("pmttime_unicor", pmttime_unicor3, "pmttime_unicor[nPMTs]/F")
tree6.Branch("charge", charge3, "charge[nPMTs]/F")

H2_CosTheta_Phi = TH2F("H2_CosTheta_Phi", "Pmt Charge; #phi(rad) ; Cos(#theta)", 100  , 0 , 6, 100 , -1, 1);
T2_CosTheta_Phi = TH2D("T2_CosTheta_Phi", "Pmt Time;   #phi(rad) ; Cos(#theta)", 100  , 0 , 6, 100 , -1, 1);
### PMT positions of the cluster
Cluster_CosTheta_Phi1 = TH2D("Charge_CosTheta_PhiCluster1",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
Cluster_CosTheta_Phi2 = TH2D("Charge_CosTheta_PhiCluster2",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
Cluster_CosTheta_Phi3 = TH2D("Charge_CosTheta_PhiCluster3",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
Time_CosTheta_Phi1 = TH2D("Time_CosTheta_PhiCluster1",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
Time_CosTheta_Phi2 = TH2D("Time_CosTheta_PhiCluster2",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
Time_CosTheta_Phi3 = TH2D("Time_CosTheta_PhiCluster3",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);

### cluster positions/no weighting
#ClusterPos_CosTheta_Phi1 = TH2D("ClusterPos_CosTheta_Phi1",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
#ClusterPos_CosTheta_Phi2 = TH2D("ClusterPos_CosTheta_Phi2",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);
#ClusterPos_CosTheta_Phi3 = TH2D("ClusterPos_CosTheta_Phi3",";#phi(rad) ; Cos(#theta)" ,100, 0, 6, 100, -1, 1);

### Jie added
### save positions of valid cluster, fill all events
EventPos_CosTheta_Phi1 = TH2D("EventPos_CosTheta_Phi1",";#phi(rad);Cos(#theta)" ,100, 0, 6, 100, -1, 1);
EventPos_RhoZ1 = TH2D("EventPos_RhoZ1",";#rho [mm]; z [mm]" ,900, 0, 900, 1800, -900, 900);
EventFprompt_QPE = TH2D("EventFprompt_QPE",";frompt; QPE",100,0,1,1000,0,1000);
EventTimeResidual = TH1D("EventTimeResidual","tRes",1000,0,1000);
### save event PMT time
PmtTime_cluster1 = TH1D("PmtTime_cluster1","cluster1; PMT time [ns];",10000,0,10000)
PmtTime_cluster2 = TH1D("PmtTime_cluster2","cluster2; PMT time [ns]",10000,0,10000)
PmtTime_cluster3 = TH1D("PmtTime_cluster3","cluster3; PMT time [ns]",10000,0,10000)
TimeRes = TH1D("TimeRes","timeRes; ns;",1000,0,1000)

Multiplicity = TH1D("nclusters",";nclusters" ,5, 0, 5);

#### PMT results ####

## data/DEAP-3600/pmt_positions.py

results = db.view('WebView/PMTPos', include_docs=True)
Phi = np.zeros(255)
CosTheta = np.zeros(255)
## PMTx = np.zeros(255)
## PMTy = np.zeros(255)
## PMTz = np.zeros(255)

index = np.zeros(255,dtype = int)

#  x = math.sin(math.radians(row.doc["locationTheta"]))*math.cos(math.radians(row.doc["locationPhi"]))
#  y = math.sin(math.radians(row.doc["locationTheta"]))*math.sin(math.radians(row.doc["locationPhi"]))
#  z = math.cos(math.radians(row.doc["locationTheta"]))

for row in results:
	if row.doc["run_range"][0] == 0:
		 #print(int(row.key), row.doc["locationPhi"])
		 Phi[int(row.key)] = row.doc["locationPhi"]*PI/180
		 CosTheta[int(row.key)] = np.cos(row.doc["locationTheta"]*PI/180)
		 index[int(row.key)] = int(row.key)

PMTlocation = []  ## Cartesian coordination
for x in zip(CosTheta, Phi):
	PMTlocation.append(TVector3(np.sqrt(1.0-x[0]**2)*np.cos(x[1]),np.sqrt(1.0-x[0]**2)*np.sin(x[1]),x[0]))

zz = [] ## save PMT id if cos()<16.3 degree
for x in PMTlocation:
	zz.append([y[1] for y in zip(PMTlocation,index) if x.Dot(y[0])>np.cos(16.3*PI/180)])
#print(zz)
n1 = []
n2 = []
n3 = []

def addPMTs(id,cluster,q_histo,t_histo):
	global qc
	for k in zz[id]:
           if charges[k] > 0 and (k not in cluster):
              cluster.append(k)
              qc = qc + charges[k]
              t_histo.Fill(Phi[k],CosTheta[k],pmtTime[k])
              q_histo.Fill(Phi[k],CosTheta[k],charges[k])
              if charges[k] > QpmtCut: # find PMT q>2 pe, it was > 10 pe
                 addPMTs(k,cluster,q_histo,t_histo)

for Files in AllFiles:
	print "\n********************************** File " , Files , " loaded... **********************************\n"
	File = TFile(Files);                    data = File.Get("T_satCorr");#TTree
	nentries = data.GetEntries();           update = int(nentries/10);
        for event in range(nentries):#range (nentries)
                #if (event+1)%update == 0:
                #       print(event+1, "analyzed...")
		data.GetEntry(event);    DS = data.ds;#(required)
		if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
		        TS = DS.GetTS(0); EV = DS.GetEV(0); CAL = DS.GetCAL(0);
		        try:CalTrigTime = CAL.GetEventTime()
		        except:continue
                        eventID = TS.GetEventID();
                        evtID[0] = eventID
                        runID[0] = DS.GetRunID();
                        subrunID[0] = DS.GetSubrunID();

	                # Low Level cuts
			if ( CAL.GetQPE() <= 200               ):  continue; # Skip events that have no charge information (primarily pre-scaled events)
			if ( CAL.GetCuts().GetCutWord()&0x31f8 ):  continue; # Low level cut, NOTE: use 0x199 if looking at data between May-Aug 2016 processed in 2016
			if ( TS.dtmTrigSrc&0x82                ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
			if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
			if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
			if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
			if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
			if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
			if ( EV.GetFmaxpe() < 0.15             ):  continue; #
			if ( EV.GetFmaxpe() > 0.4              ):  continue; #
			if ( CAL.GetFprompt() > 0.58           ):  continue;
			H2_CosTheta_Phi.Reset()
                        T2_CosTheta_Phi.Reset();
                        ### position
                        posx = EV.mblikelihood.GetPosition().X(); posy = EV.mblikelihood.GetPosition().Y(); posz = EV.mblikelihood.GetPosition().Z();
                        posPhi = EV.mblikelihood.GetPosition().Phi(); posCosTheta = EV.mblikelihood.GetPosition().CosTheta();
                        evtx[0] = posx; evty[0] = posy; evtz[0] = posz;

                        charges = np.zeros(255); id = []; pmtTime = np.zeros(255); ### time of each PMT
                        pmtTime1 = np.zeros(255); pmtTime2 = np.zeros(255);
			for i in range (CAL.GetPMTCount()): ### loop PMTs
                                PMT = CAL.GetPMT(i)
                                PMTid = PMT.id
                                q = PMT.qPrompt #haojie
                                # specharge = PMTInfoUtil().GetSPECharge(PMTid) # James
                                charges[PMTid] = q
				id.append(PMTid)
                                pmtTime[PMTid] = 100000 #pulse subpeak time
                                pmtTime1[PMTid] = 100000 #pulse subpeak time uniform
                                pmtTime2[PMTid] = 100000 #pulse subpeak time TOF

                                for j in range (PMT.GetPulseCount()): ### loop pulses in a single PMT
                                        PULSE = PMT.GetPulse(j)
                                        for ipeak in range(PULSE.GetSubpeakCount()):
                                                time = PULSE.GetSubpeakTime(ipeak) ##
                                                time1 = PULSE.GetSubpeakTimeUniform(ipeak)
                                                time2 = PULSE.GetSubpeakTimeUniformWithTOFCorrection(ipeak)
						if time < pmtTime[PMTid]:
							pmtTime[PMTid] = time
                                                if time1 < pmtTime1[PMTid]:
                                                        pmtTime1[PMTid] = time1
                                                if time2 < pmtTime2[PMTid]:        
                                                        pmtTime2[PMTid] = time2
				H2_CosTheta_Phi.Fill(Phi[PMTid], CosTheta[PMTid], q)
				#print(PMTid,q,Phi[PMTid],Theta[PMTid])
                        ### loop PMTs and save into tuple         
			charges = list(charges)
			pmtTime = list(pmtTime)
                        ### fill this event
			for i in range(255):
                           if charges[i] > 0: ### Fill all the PMTs in this event
                                T2_CosTheta_Phi.Fill(Phi[i], CosTheta[i], pmtTime[i]); 

                        ### fprompt vs QPE
                        EventFprompt_QPE.Reset()
                        EventPos_CosTheta_Phi1.Reset()
                        EventPos_RhoZ1.Reset()
                        Cluster_CosTheta_Phi1.Reset()
                        Cluster_CosTheta_Phi2.Reset()
                        Cluster_CosTheta_Phi3.Reset()

                        Cluster_CosTheta_Phi1.Reset()
			Cluster_CosTheta_Phi2.Reset()
			Cluster_CosTheta_Phi3.Reset()
			Time_CosTheta_Phi1.Reset()
			Time_CosTheta_Phi2.Reset()
			Time_CosTheta_Phi3.Reset()
                        
                        ### look for cluster 1
                        cluster1 = []
                        qcluster1 = 0
                        array_pmtTime = np.array(pmtTime)
                        check_pmtTime = array_pmtTime[array_pmtTime != 0] ### remove zeros!!!
                        earliestTime = min(check_pmtTime)
                        id1 = pmtTime.index(earliestTime) ### !!!! Jie: find the first PMT with earliest hittime (min(time)) 
                        print "!!! cluster 1 min ", id1, earliestTime 
                        cluster1.append(id1)
                        qc = 0
			addPMTs(id1,cluster1,Cluster_CosTheta_Phi1,Time_CosTheta_Phi1)
                        print("cluster1 PMTid",cluster1)
                        qcluster1 = qc+charges[id1]
                        print("total Q cluster1",qcluster1)
                        count_ncluster = 0
                        if qcluster1 > QclusterCut: # was 100 for James
				count_ncluster = 1 ## count for single cluster
                                ncluster[0] = count_ncluster
                                nPMTs1[0] = len(cluster1)
                                kk = 0
                                #for pmtid in cluster1:
                                for kk in xrange(nPMTs1[0]):
                                    pmtid = cluster1[kk]
                                    #print pmtid, nPMTs1, Phi[pmtid], kk
                                    ## values below are filling trees
                                    pmtPhi1[kk] = Phi[pmtid]
                                    pmtCosTheta1[kk] = CosTheta[pmtid]
                                    pmttime1[kk] = pmtTime[pmtid]
                                    pmttime_uni1[kk] = pmtTime1[pmtid] 
                                    pmttime_unicor1[kk] = pmtTime2[pmtid]
                                    charge1[kk] = charges[pmtid]
                                tree4.Fill()
                                # continue
                        else:
				count_ncluster = 0
                                # continue

                        ### get event position and qpe
                        fprompt = CAL.GetFprompt() ; fpromptVal[0] = fprompt
                        qpe = CAL.GetQPE()         ; qpeVal[0] = qpe

                        EventFprompt_QPE.Fill(fprompt,qpe)
                        EventPos_CosTheta_Phi1.Fill(posPhi,posCosTheta)
                        EventPos_RhoZ1.Fill(sqrt(posx*posx+posy*posy),posz)

                        Cluster_CosTheta_Phi1.Fill(Phi[id1],CosTheta[id1],charges[id1])
                        Time_CosTheta_Phi1.Fill(Phi[id1],CosTheta[id1],pmtTime[id1])
                        # ClusterPos_CosTheta_Phi1.Fill(Phi[id1],CosTheta[id1])
                        print("cluster1 pmtTime",pmtTime[id1])
                        #PmtTime_cluster1.Fill(pmtTime[id1])

                        ### after finding, clear these PMTs
                        for pmtid in cluster1:
                                charges[pmtid] = 0 ## find the cluster1 and setting these PMTs' charge to 0, clear them
                                pmtTime[pmtid] = 0 ## find the cluster1 and setting these PMTs' time to 0, clear them

                        ### look for cluster 2 
                        cluster2 = []
			qcluster2 = 0
                        array_pmtTime = np.array(pmtTime)
                        check_pmtTime = array_pmtTime[array_pmtTime != 0] ### remove zeros!!!
                        earliestTime = min(check_pmtTime)
                        id2 = pmtTime.index(earliestTime) ### !!!! Jie: find the first PMT with earliest hittime (min(time))
                        print "!!! cluster 2 min ", id2, earliestTime
                        cluster2.append(id2)
			qc = 0
			addPMTs(id2,cluster2,Cluster_CosTheta_Phi2,Time_CosTheta_Phi2);
			Cluster_CosTheta_Phi2.Fill(Phi[id2],CosTheta[id2],charges[id2]);
			Time_CosTheta_Phi2.Fill(Phi[id2],CosTheta[id2],pmtTime[id2]);
			print("cluster2 PMTid",cluster2)
			qcluster2 = qc + charges[id2]
			print("total Q cluster2",qcluster2)
			if qcluster2 > QclusterCut: # was 100
				count_ncluster = count_ncluster + 1 ## count for double cluster
                                ncluster[0] = count_ncluster
                                #ClusterPos_CosTheta_Phi2.Fill(Phi[id2],CosTheta[id2]) ## PMT pos, no weighting
                                nPMTs2[0] = len(cluster2)
                                kk = 0
                                #for pmtid in cluster2:
                                for kk in xrange(nPMTs2[0]):
                                    pmtid = cluster2[kk]
                                    pmtPhi2[kk] = Phi[pmtid]
                                    pmtCosTheta2[kk] = CosTheta[pmtid]
                                    pmttime2[kk] = pmtTime[pmtid]
                                    pmttime_uni2[kk] = pmtTime1[pmtid]
                                    pmttime_unicor2[kk] = pmtTime2[pmtid]
                                    charge2[kk] = charges[pmtid]
                                    #print pmtid, nPMTs2[0], pmtCosTheta2[kk], kk
                                tree5.Fill() ## fill the trees in order (2nd), not the double-cluster events!

                        for pmtid in cluster2:
				charges[pmtid] = 0 ## find the cluster2 and setting these PMTs' charge to 0, clear them!
                                pmtTime[pmtid] = 0

                        ### look for cluster 3 
                        cluster3 = []
			qcluster3 = 0
                        array_pmtTime = np.array(pmtTime)
                        check_pmtTime = array_pmtTime[array_pmtTime != 0] ### remove zeros!!!
                        earliestTime = min(check_pmtTime)
                        id3 = pmtTime.index(earliestTime) ### !!!! Jie: find the first PMT with earliest hittime (min(time))
                        print "!!! cluster 3 min ", id3, earliestTime
			cluster3.append(id3)
			qc = 0
			addPMTs(id3,cluster3,Cluster_CosTheta_Phi3,Time_CosTheta_Phi3)
			print("cluster3 PMTid",cluster3)
			Cluster_CosTheta_Phi3.Fill(Phi[id3],CosTheta[id3],charges[id3]);
			Time_CosTheta_Phi3.Fill(Phi[id3],CosTheta[id3],pmtTime[id3]);
			qcluster3 = qc + charges[id3]
			print("total Q cluster3",qcluster3)
			if qcluster3 > QclusterCut: # was 100
				count_ncluster = count_ncluster + 1 ## count for tripple cluster
                                ncluster[0] = count_ncluster
                                # ClusterPos_CosTheta_Phi3.Fill(Phi[id3],CosTheta[id3])
                                nPMTs3[0] = len(cluster3)
                                kk = 0
                                #for pmtid in cluster3:
                                for kk in xrange(nPMTs3[0]):
                                    pmtid = cluster3[kk]
                                    pmtPhi3[kk] = Phi[pmtid]
                                    pmtCosTheta3[kk] = CosTheta[pmtid]
                                    pmttime3[kk] = pmtTime[pmtid]
                                    pmttime_uni3[kk] = pmtTime1[pmtid]
                                    pmttime_unicor3[kk] = pmtTime2[pmtid]
                                    charge3[kk] = charges[pmtid]
                                tree6.Fill() ## fill the trees in order (3rd), not the triple-cluster events!

                        for pmtid in cluster3:
				charges[pmtid] = 0 ## find the cluster3 and setting these PMTs' charge to 0, clear them!
                                pmtTime[pmtid] = 0

			Multiplicity.Fill(ncluster[0])
			if ncluster[0] == 1:
                                print "event",eventID, " cluster=1"
                                tree1.Fill() ## fill the tree for single-cluster event
				n1.append(ncluster[0])
			elif ncluster[0] == 2:
                                print "event",eventID, " cluster=2"
                                tree2.Fill() ## fill the tree for double-cluster event
				n2.append(ncluster[0])
			elif ncluster[0] == 3:
                                print "event",eventID, " cluster=3"
                                tree3.Fill() ## fill the tree for triple-cluster event
				n3.append(ncluster[0])
  
        ### end of looping events
                fout.cd()
                H2_CosTheta_Phi.Write("H2_CosTheta_Phi_%d"%eventID)
                Cluster_CosTheta_Phi1.Write("Charge_CosTheta_PhiCluster1_%d"%eventID)
                Cluster_CosTheta_Phi2.Write("Charge_CosTheta_PhiCluster2_%d"%eventID)
                Cluster_CosTheta_Phi3.Write("Charge_CosTheta_PhiCluster3_%d"%eventID)

                #ClusterPos_CosTheta_Phi1.Write("Clusterpos_CosTheta_Phi1_%d"%eventID)
                #ClusterPos_CosTheta_Phi2.Write("Clusterpos_CosTheta_Phi2_%d"%eventID)
                #ClusterPos_CosTheta_Phi3.Write("Clusterpos_CosTheta_Phi3_%d"%eventID)

                T2_CosTheta_Phi.Write("T2_CosTheta_Phi_%d"%eventID)
                Time_CosTheta_Phi1.Write("Time_CosTheta_Phi1_%d"%eventID)
                Time_CosTheta_Phi2.Write("Time_CosTheta_Phi2_%d"%eventID)
                Time_CosTheta_Phi3.Write("Time_CosTheta_Phi3_%d"%eventID)

                EventPos_CosTheta_Phi1.Write("EventPos_CosTheta_Phi1_%d"%eventID)
                EventPos_RhoZ1.Write("EventPos_RhoZ1_%d"%eventID)
                EventFprompt_QPE.Write("EventFprompt_QPE_%d"%eventID)
                #TimeRes.Write("TimeRes_%d"%eventID)

        fout.cd()
        tree1.Write()
        tree2.Write()
        tree3.Write()
        tree4.Write()
        tree5.Write()
        tree6.Write()
        Multiplicity.Write("nclusters")

	#_%d"%event

	##########################################################################################################################################
	##########################################################################################################################################
	File.Close()
fout.Close()