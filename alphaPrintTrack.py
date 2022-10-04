from ROOT import *
from rat import *
import couchdb

import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

#			EV s	   km/s	to conver to nm
hc = 4.1357*10e-15 * 300000 * 1e5
PI = np.pi
Malpha = 3727 #MeV
m_alpha_inKg = 6.64*1e-27 # kg
Colors = ["red", "blue", "green", "purple"]
### use charge-cluster/james
evtID_single = [6,8,14,16,17,24,27,29,32,51,54,56,70,73,75,77,79,86,87,88,92] #[7]  0,4,5,11,14,22,41]
evtID_double = [4,12,13,26,30,31,39,47,82,84]#7,8,10,51]  #54, 79, 80, 99,115]
evtID_triple = []

list_center = [(0,0,0)] # [(1,2,3),(-4,-5,6), (5,5,6)]
list_radius = [850] # [1,2,1]

filename = "trackMCusher_run30695_Po210_vacuum_avBulk_thin0p05_1e4evts_rmStepLimiter.root"
File = TFile(filename)
filenewname = "SaveTrack_"+filename

fsave = TFile(filenewname,"recreate")

### loading PMT info
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
# Dir=sys.argv[1]
#AllFiles=glob.glob(Dir)
#AllFiles = sorted(AllFiles)
results = db.view('WebView/PMTPos', include_docs=True)
Phi = np.zeros(255)
CosTheta = np.zeros(255)
#### PMTx = np.zeros(255); PMTy = np.zeros(255); PMTz = np.zeros(255)
index = np.zeros(255,dtype = int)
##  x = math.sin(math.radians(row.doc["locationTheta"]))*math.cos(math.radians(row.doc["locationPhi"]))
##  y = math.sin(math.radians(row.doc["locationTheta"]))*math.sin(math.radians(row.doc["locationPhi"]))
##  z = math.cos(math.radians(row.doc["locationTheta"]))
for row in results:
    if row.doc["run_range"][0] == 0:
        #print(int(row.key), row.doc["locationPhi"])
        Phi[int(row.key)] = row.doc["locationPhi"]*PI/180
        CosTheta[int(row.key)] = np.cos(row.doc["locationTheta"]*PI/180)
        index[int(row.key)] = int(row.key)
#############################################

### draw AV sphere
def plt_sphere(list_center, list_radius):
    for c, r in zip(list_center, list_radius):
        ax = fig.gca(projection='3d')
        # draw sphere
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
        x = r*np.cos(u)*np.sin(v)
        y = r*np.sin(u)*np.sin(v)
        z = r*np.cos(v)
        ax.plot_surface(x-c[0], y-c[1], z-c[2], color='blue', alpha=0.1)

### use half-charge
#evtID_double = [8,10,14,34,44,51]
#evtID_triple = [0,4,5,7,11,22,35,37,40,41,47]

### for TPB surface
#evtID_single = [31,40,41,42,54,77,81,85]
#evtID_double = [45,69]

## for alpha AV 500 simulations
# evtID_double = [21, 27, 44, 61, 74, 83,144,154,171,181,231,264,317,326]

##trackMCusher_Po210_vacuum_tpbSurface_pressure0p005_100evts_rmStepLimiter.root");

HdistVsTime = TH2F("HdistVsTime","", 100, 0, 100, 100, 0, 1700)
Hspeed = TH1F("Hspeed","",100,0,100)
Hspeed_vac = TH1F("Hspeed_vac","",100,0,100)

# hit point 1 (going out tpb)
H2D_Hit1_CosTheta_Phi = TH2D("H2D_Hit1_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# hit point 2 (going in tpb)
H2D_Hit2_CosTheta_Phi = TH2D("H2D_Hit2_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# mc pmts surrounding hit point 1
H2D_Hit1mcPMT_CosTheta_Phi = TH2D("H2D_Hit1mcPMT_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# mc pmts surrounding hit point 2
H2D_Hit2mcPMT_CosTheta_Phi = TH2D("H2D_Hit2mcPMT_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# all MC PMTs
H2D_mcPMT_CosTheta_Phi = TH2D("H2D_mcPMT_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# all trigger PMTs
H2D_trigPMT_CosTheta_Phi = TH2D("H2D_trigPMT_CosTheta_Phi", ";#phi(rad) ; Cos(#theta)", 50,0, 2*PI, 50 , -1, 1);
# all trigger PMTs vs charge
H2D_trigPMT_q_CosTheta_Phi = TH2D("H2D_trigPMT_q_CosTheta_Phi", "charge; #phi(rad) ; Cos(#theta);", 50,0, 2*PI, 50 , -1, 1);
# all trigger PMTs vs hitTime
H2D_trigPMT_t_CosTheta_Phi = TH2D("H2D_trigPMT_t_CosTheta_Phi", "time; #phi(rad) ; Cos(#theta);", 50,0, 2*PI, 50 , -1, 1);

## kinetic energy = 1/2mv^2
H_Ek_tot = TH1D("H_Ek_tot","kinetic energy (E_{k}=1/2mv^{2}) for double-cluster, total", 100,0,10)
H_Ek_vac = TH1D("H_Ek_vac","kinetic energy (E_{k}=1/2mv^{2}) in vacuum for double-cluster", 100,0,10)
H_Ei = TH1D("H_Ei","E at track start", 100,0,20)
H_Ef = TH1D("H_Ef","E at track end -1", 100,0,20)
H_Ehit = TH1D("H_Ehit","E at hitpoint in TPB", 100,0,20)
H_deltaEk = TH1D("H_deltaEk","E_{vac} - E_{total} for double-cluster", 100,0,20)

H_deltaEk_vac = TH1D("H_deltaEk_vac","E_{1,vac} - E_{2,vac} (1/2mv^{2}) for double-cluster", 100,0,10)
H_deltaEp_vac = TH1D("H_deltaEp_vac","E_{1,vac} - E_{2,vac} (#frac{p^{2}}{2m}) for double-cluster", 100,0,10)

H_energy_all = TH1D("H_energy_all", "continuous energy spectrum from all steps (for double-cluster)", 100,0,20)

## T = File.Get("T")
T = File.Get("T_satCorr")
nentries = T.GetEntries();
###############################################################################################################################################################################################
# DC 	   = RAT.DetectorConfig.GetDetectorConfig(0);#Connects to https://deapdb.physics.carleton.ca/deapdb
# nPMTs  = DC.GetNumPMTs();#Number of PMTs
# pmtPos = DC.GetPMTPositions();# A list of 255 TVector3's containing the position of the 255 PMTs
# print "\n\n\n\n\n\n\n\n************************************************* Detector Configuration loaded *************************************************\n\n\n\n\n\n\n\n"
###############################################################################################################################################################################################

countTrig = 0
countTrig_afterCuts = 0

print "Started reading events..."
for event in range (nentries):
	if (event)%10 == 0:
		print (event), "Analyzed..."
	T.GetEntry(event); DS = T.ds;
	if DS.GetTSCount() > 0 and DS.GetEVCount() > 0:
		TS = DS.GetTS(0);EV = DS.GetEV(0);MC = DS.GetMC();
                CAL = DS.GetCAL(0);### for calibrated, triggered
                eventID = TS.GetEventID()

                countTrig += 1
                ## if eventID not in evtID_single:
                #if eventID not in evtID_double:
                ### if eventID not in evtID_triple:
                #    continue
                #else:
                #    print "checking double-cluster events, eventID =", eventID, "event entry=",event

                # Low Level cuts
                #print "CAL.GetQPE()", CAL.GetQPE()
                if ( CAL.GetQPE() <= 200               ):  continue; # Skip events that have no charge information (primarily pre-scaled events)
                #### if ( CAL.GetCuts().GetCutWord()&0x31f8 ):  continue; # Low level cut, NOTE: use 0x199 if looking at data between May-Aug 2016 processed in 2016
                #### if ( TS.dtmTrigSrc&0x82                ):  continue; # Low level cut, remove internal/external/ periodic and muon veto triggers
                if ( CAL.GetSubeventCount() > 1        ):  continue; # pile-up cut, removes coincidence events
                if ( CAL.GetEventTime() <= 2250        ):  continue; # pre-trigger pile up cut, not typically identified by subeventN
                if ( CAL.GetEventTime() >= 2700        ):  continue; # post-trigger pile up cut, strongly correlated with subeventN
                if ( CAL.numEarlyPulses > 3            ):  continue; # Cut on a significant amount of early light pulses in the first 1600 ns of the waveform
                ### if ( EV.deltat <= 20000                ):  continue; # deltat cut on 20 us
                if ( EV.GetFmaxpe() < 0.15             ):  continue; #
                if ( EV.GetFmaxpe() > 0.4              ):  continue; #
                if ( CAL.GetFprompt() > 0.58           ):  continue;

                countTrig_afterCuts += 1

                H2D_Hit1_CosTheta_Phi.Reset()
                H2D_Hit2_CosTheta_Phi.Reset()
                H2D_Hit1mcPMT_CosTheta_Phi.Reset()
                H2D_Hit2mcPMT_CosTheta_Phi.Reset()
                H2D_trigPMT_CosTheta_Phi.Reset()
                H2D_mcPMT_CosTheta_Phi.Reset()
                H2D_trigPMT_q_CosTheta_Phi.Reset()
                H2D_trigPMT_t_CosTheta_Phi.Reset()

		# Save all the tracks in python dictionary, key is track ID and value is the track itself
		Track_Dictionary = {};
		for i in range(MC.GetMCTrackCount()):
			Track = MC.GetMCTrack(i); Track_Dictionary[ Track.trackID ] = Track;
                ## print Track_Dictionary
################################################################################################################################################################################################
                charges = np.zeros(255); id = []; pmtTime = np.zeros(255);
                for i in range (CAL.GetPMTCount()): ### loop PMTs
                        PMT = CAL.GetPMT(i)
                        pmtid = PMT.id
                        q = PMT.qPrompt
                        charges[pmtid] = q
                        H2D_trigPMT_CosTheta_Phi.Fill(Phi[pmtid], CosTheta[pmtid])
                        H2D_trigPMT_q_CosTheta_Phi.Fill(Phi[pmtid], CosTheta[pmtid], q)
                        pmtTime[pmtid] = 100000
                        ## print "!!! PMT get pulse count", PMT.GetPulseCount()
                        for j in range (PMT.GetPulseCount()):
                                PULSE = PMT.GetPulse(j)
                                for ipeak in range(PULSE.GetSubpeakCount()):
                                          time = PULSE.GetSubpeakTime(ipeak) ##
                                          if time < pmtTime[pmtid]:
                                                pmtTime[pmtid] = time
                charges = list(charges)
                pmtTime = list(pmtTime)

                aveTime = sum(pmtTime)/len(pmtTime)
                earliestTime = min(filter(lambda x: x>0, pmtTime))
                # print "earliest", earliestTime

                ### fill this event
                # print "pmt hit-time"
                for i in range(255):
                       if charges[i] > 0: ### Fill all the PMTs in this event
                            #print pmtTime[i],
                            H2D_trigPMT_t_CosTheta_Phi.Fill(Phi[i], CosTheta[i], pmtTime[i]-earliestTime) #2400)
                #print "\n pmt charges"
                #for i in range(255):
                #       if charges[i] > 0: ### Fill all the PMTs in this event
                            #print charges[i],
                #print "\n pmt position"
                #for i in range(255):
                #       if charges[i] > 0:
                            #print Phi[i], CosTheta[i], " " 

                ### Loop MC PMTs
                # print "!!!MC.GetPMTCount()", MC.GetPMTCount()
                for j in range(MC.GetPMTCount()):
                       mcpmt = MC.GetPMT(j)
                       H2D_mcPMT_CosTheta_Phi.Fill(Phi[mcpmt.id], CosTheta[mcpmt.id])

		for j in range(MC.GetPMTCount()):
                       # print "MC PMTcount", MC.GetPMTCount()                    
                       mcpmt = MC.GetPMT(j)
                       Cathodtimes = []; Photons = []
                       for l in range(mcpmt.GetMCPhotonCount()):
                       	        Cathodtimes.append(mcpmt.GetMCPhoton(l).GetCathodeTime());  Photons.append(mcpmt.GetMCPhoton(l));
                       
                       PMTid1 = []
                       PMTid2 = []
                       # print mcpmt.id
                       # Earliest Photons in each PMT!
                       # This is what I wanted to do! You choose your main track you want to trace back!
                       Index = Cathodtimes.index(min(Cathodtimes));  Earliest_photon = Photons[Index];
                       PhotonTrackID = Earliest_photon.trackID;  ParentTrackID = PhotonTrackID;
                       countPMT = 0
                       # print "earliest trackID", PhotonTrackID
                       try:
		        	# List of all the parent tracks starting with the chosen track!
		        	ParentTracks = [Track_Dictionary[PhotonTrackID]];
		        	# Get the parent ID of each parent until you reach the source
		        	while ParentTrackID > 1:# Not the event itself
		        		ParentTrackID = Track_Dictionary[ParentTrackID].parentID; ParentTracks.append(Track_Dictionary[ParentTrackID]);

		        	fig = plt.figure(figsize=(20,20));  ax = plt.axes(projection='3d'); Title = "";
		        	# ax.set_title( "T_res = " + str (T_RES) + "; CellID = " + str(CellID) )
		        	print "##############################################################################################################################################"
		        	Vectors = []; Time = [];
		        	for processes in range(len(ParentTracks)):
                                       # We need to loop backward since we traced the track backward; i.e the first track
                                       # in the list is the final track!
                                       Track = ParentTracks[-processes-1]
                                       ### Note: ensure it is looking for alpha!
                                       if Track.particleName != "alpha":
                                            continue
                                       print "\n New Track: ",Track.particleName,"\n"
                                       # Loop over steps and print some stuff nicely!
                                       checkVolume = []
                                       totalSteps = Track.GetMCTrackStepCount()
                                       Ei = 0 ## initial E at the start of steps
                                       Ef = 0 ## final E at the end of steps
                                       Ehit = 0 ## when the alpha hits the TPB
                                       Ehit1 = 0; # tpb->vacuum
                                       Ehit2 = 0; # vacuum->tpb
                                       Ek_hit1 = 0;# use velocity; Ek=1/2 m v^2, v = (x-x_prev)/(t-t_prev);
                                       Ek_hit2 = 0; 
                                       xhit1 = 0; yhit1 = 0; zhit1 = 0; thit1 = 0; ## tpb->vacuum hit point1
                                       xhit2 = 0; yhit2 = 0; zhit2 = 0; thit2 = 0; ## vacuum->tpb hit point2

                                       for steps in range(totalSteps): ## loop track steps
                                             momentum = Track.GetMCTrackStep(steps).GetMomentum().Mag()
                                             print "p=", momentum, "p^2=",momentum*momentum, "p^2/(2m)=",momentum*momentum/2/Malpha
                                             energyKinetic = (momentum*momentum)/2/Malpha
                                             step_pos = Track.GetMCTrackStep(steps).GetEndpoint()
                                       	     x = step_pos.x(); y = step_pos.y(); z = step_pos.z();
                                             pos_phi = step_pos.Phi(); pos_cost = step_pos.CosTheta();
                                             r = np.sqrt(x*x+y*y+z*z);
                                             Temp = TVector3(); 
                                             Temp.SetXYZ(x,y,z);
                                       	     Vectors.append(Temp);
                                             ## global time at this track step
                                             time_track = Track.GetMCTrackStep(steps).GetGlobalTime()
                                             Time.append(time_track);

                                             H_energy_all.Fill(energyKinetic)

                                             volume = Track.GetMCTrackStep(steps).GetVolume()
                                             if steps == 0: ## initial energy
                                                 Ei = energyKinetic
                                             if steps == totalSteps-1: ## final energy
                                                 Ef = energyKinetic
                                             # print "????? Energy",steps, Ei, Ef

                                             if volume == "tpb_bulk":
                                                 ## print "!!!!! step volume =", volume
                                                 if steps+1<totalSteps:
                                                     next_volume = Track.GetMCTrackStep(steps+1).GetVolume()
                                                     if next_volume == "cryoliquid": ### !!! for tpb->vacuum
                                                          Ehit = energyKinetic
                                                          Ehit1 = energyKinetic
                                                          print "????? energy hit point 1", Ehit
                                                          phi1 = pos_phi; cost1 = pos_cost
                                                          xhit1 = x; yhit1 = y; zhit1 = z; thit1 = time_track;
                                                          #### print "??????????? hit1=(",xhit1,yhit1,zhit1,"),",thit1
                                                          if phi1<0:
                                                              phi1 = phi1+2*PI
                                                          H2D_Hit1_CosTheta_Phi.Fill( phi1, cost1 )
                                                          #print "!!!! from tpb -> vacuum, going out"
                                                          print "!!!! tpb->vacuum: eventID=", eventID, "(phi,costheta)=(",phi1, ",", step_pos.CosTheta(),")", "Ek=", Ehit
                                                          ## calculate Ek = 1/2mv^2
                                                          step_posPrev = Track.GetMCTrackStep(steps-1).GetEndpoint()
                                                          deltaX = (step_pos-step_posPrev).Mag()
                                                          time_trackPrev = Track.GetMCTrackStep(steps-1).GetGlobalTime()
                                                          deltaT = time_track - time_trackPrev
                                                          if deltaT != 0:
                                                               vv = deltaX/deltaT
                                                               Ek_hit1 = 0.5*m_alpha_inKg*(vv*1e6)*(vv*1e6)/(1.6e-13)
                                                           
                                                 ## print "!!!!! step volume =", volume
                                                 if steps-1>0:
                                                     prev_volume = Track.GetMCTrackStep(steps-1).GetVolume()
                                                     if prev_volume == "cryoliquid": ### !!! for vacuum->tpb
                                                          Ehit = energyKinetic
                                                          Ehit2 = energyKinetic
                                                          phi2 = pos_phi; cost2 = pos_cost;
                                                          xhit2 = x; yhit2 = y; zhit2 = z; thit2 = time_track;
                                                          ### print "??????????? hit2=(",xhit2,yhit2,zhit2,"),",thit2
                                                          if phi2<0:
                                                              phi2 = phi2+2*PI
                                                          H2D_Hit2_CosTheta_Phi.Fill( phi2, cost2 )
                                                          #print "!!!! from vacuum -> tpb, going in"
                                                          print "!!!! vacuum->tpb: (phi,costheta)=", "(",phi2, ",", step_pos.CosTheta(),")", "Ek=", Ehit
                                                          ## calculate Ek = 1/2mv^2
                                                          step_posPrev = Track.GetMCTrackStep(steps-1).GetEndpoint()
                                                          deltaX = (step_pos-step_posPrev).Mag()
                                                          time_trackPrev = Track.GetMCTrackStep(steps-1).GetGlobalTime()
                                                          deltaT = time_track - time_trackPrev
                                                          if deltaT != 0:
                                                               vv = deltaX/deltaT
                                                               Ek_hit2 = 0.5*m_alpha_inKg*(vv*1e6)*(vv*1e6)/(1.6e-13)

                                             #if Track.particleName == "alpha" and len(Vectors) > 1 and (Time[-1]-Time[-2]) != 0 :#	process Volume time lambda  x   y   z
                                             #   print "process, volume, globalTime, x, y, z, r, diffPos, DeltaX/DeltaT" 
                                       	     #   print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
                                       	     #   Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), "    ------",
                                       	     #   round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), round((Vectors[-1]-Vectors[-2]).Mag()/(Time[-1]-Time[-2]),7)
                                       	     #   )
                                       
                                       	     #elif Track.particleName == "alpha" and len(Vectors) > 1 and (Time[-1]-Time[-2]) != 0 :#	process Volume time lambda  x   y   z
                                       	     #	print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
                                       	     #	Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), round(hc/Track.GetMCTrackStep(steps).GetKE(),1),
                                       	     #	round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), round((Vectors[-1]-Vectors[-2]).Mag()/(Time[-1]-Time[-2]),7)
                                       	     #	)
                                       	     #elif Track.particleName == "opticalphoton" and len(Vectors) > 1 and (Time[-1]-Time[-2]) == 0 :#	process Volume time lambda  x   y   z
                                       	     #   print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
                                       	     #   Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), round(hc/Track.GetMCTrackStep(steps).GetKE(),1),
                                       	     #   round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), "    ------"
                                       	     #   )
                                       ### Plot the points on the track                                       
                                       xdata = []; ydata = []; zdata = []; trackStepTime = []; trackLength = 0;
                                       Track = ParentTracks[-processes-1]
                                       Title += "eventID: "+str(eventID)+" "+Track.particleName + "(steps = " + str(Track.GetMCTrackStepCount()) + ") -->";
                                       ## go through steps
                                       for steps in range(Track.GetMCTrackStepCount()):
                                            xdata.append( Track.GetMCTrackStep(steps).GetEndpoint().x() );
                                            ydata.append( Track.GetMCTrackStep(steps).GetEndpoint().y() );
                                            zdata.append( Track.GetMCTrackStep(steps).GetEndpoint().z() );
                                            trackStepTime.append(Track.GetMCTrackStep(steps).GetGlobalTime());

                                       for i in range(len(xdata)):
                                           if (len(xdata) != len(ydata)) or (len(xdata) != len(ydata)) or (len(xdata) != len(ydata)):
                                               print "data lengths are wrong"
                                               break
                                           if i>0:
                                               #print "data points = ", len(xdata)
                                               deltaX = xdata[i] - xdata[i-1]
                                               deltaY = ydata[i] - ydata[i-1]
                                               deltaZ = zdata[i] - zdata[i-1]
                                               deltaR = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)
                                               #print deltaR, trackLength
                                               trackLength += deltaR
                                       posDiff = sqrt( (xdata[-1]-xdata[0])*(xdata[-1]-xdata[0]) + (ydata[-1]-ydata[0])*(ydata[-1]-ydata[0]) + (zdata[-1]-zdata[0])*(zdata[-1]-zdata[0]) )
                                       ax.plot3D(xdata, ydata, zdata, Colors[processes]);
                                       Text = "eventID="+str(eventID)+"; process="+str(processes) + ", " + Track.particleName;
                                       #print Text, "trackLength =", round(trackLength,2), "mm, t=", round(trackStepTime[-1],2), "ns, speed =", round(trackLength/trackStepTime[-1],2), "mm/ns", " |Xf-Xi|=", round(posDiff,2), "mm, speed2 =", round(posDiff/trackStepTime[-1],2), "mm/ns" 
                                       
                                       ## save info into histograms
                                       velocity = trackLength/trackStepTime[-1];
                                       Hspeed.Fill(velocity)
                                       Ek = 0.5*m_alpha_inKg*(velocity*1e6)*(velocity*1e6)/(1.6e-13) # mm/ns to m/s, J to MeV
                                       
                                       velocity_vacuum = 0
                                      
                                       if (thit2-thit1) != 0:
                                          Poshit1 = TVector3(xhit1,yhit1,zhit1)
                                          Poshit2 = TVector3(xhit2,yhit2,zhit2)
                                          #print "match??", xhit1, yhit1, zhit1," ", xhit2, yhit2, zhit2, " ", (Poshit1-Poshit2).Mag(), "mm", (Poshit1-Poshit2).Mag()/(thit2 - thit1)
                                          trackLengthVacuum = (Poshit1-Poshit2).Mag()
                                          velocity_vacuum = trackLengthVacuum/(thit2 - thit1)
                                        
                                       Ek_vac = 0.5*m_alpha_inKg*(velocity_vacuum*1e6)*(velocity_vacuum*1e6) # mm/ns to m/s
                                       Ek_vac = Ek_vac/(1.6*1e-13)

                                       ## H_Ek_tot.Fill(Ek)  ## here for E of total events, including single+double-cluster
                                       H_Ei.Fill(Ei)
                                       H_Ef.Fill(Ef)
                                       H_Ehit.Fill(Ehit) ## all the energies deposited at two hit points
                                       if velocity_vacuum != 0: 
                                            Hspeed_vac.Fill(velocity_vacuum)
                                            H_Ek_vac.Fill(Ek_vac)
                                            H_Ek_tot.Fill(Ek)
                                            H_deltaEk.Fill(Ek_vac - Ek)

                                       if Ek_hit1!=0 and Ek_hit2!=0:
                                            H_deltaEk_vac.Fill(Ek_hit1-Ek_hit2)
                                            H_deltaEp_vac.Fill(Ehit1 - Ehit2)

                                       HdistVsTime.Fill( trackStepTime[-1], trackLength )

                                       for steps in range(Track.GetMCTrackStepCount()):
                                            #print "drawing point for step", steps
                                            if processes > 0:
                                       		  ax.scatter( Track.GetMCTrackStep(steps).GetEndpoint().x(), Track.GetMCTrackStep(steps).GetEndpoint().y(), Track.GetMCTrackStep(steps).GetEndpoint().z() );
                                #print "*****************************************************************************************************************************************************************************************************"
                                #print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format( "Process", "Volume", "     Time", "   Lambda", "\tX", "  Y", "  Z", "  R", "  Distance", "Speed" )
                                #print "*****************************************************************************************************************************************************************************************************"
                                countPMT = 1
		                ### Wait for the user to read the track and then plot it!
		                # raw_input("Press Enter to continue...")
		                Text = str(processes) + "," + Track.particleName;
                                ax.text2D(0.01, 0.99, Title, transform=ax.transAxes);
                                Text2 = "track length = "+str(round(trackLength,2))+"mm; speed = "+str(round(trackLength/trackStepTime[-1],2))+"mm/ns"
                                Text3 = "position difference = "+str(round(posDiff,2))+"mm; speed2 = "+str(round(posDiff/trackStepTime[-1],2))+"mm/ns"
                                ax.text2D(0.01, 0.95, Text2, transform=ax.transAxes)
                                ax.text2D(0.01, 0.91, Text3, transform=ax.transAxes)
                                #ax.text2D(0.01, 0.95, round(Cathodtimes[Index],2), transform=ax.transAxes)
                                ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.set_zlabel('Z [mm]')

                                u = np.linspace(0, 2 * np.pi, 100)
                                v = np.linspace(0, np.pi, 100)
                                radius = 850.
                                x = radius * np.outer(np.cos(u), np.sin(v))
                                y = radius * np.outer(np.sin(u), np.sin(v))
                                z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
                                #for i in range(2):
                                # ax.plot_surface(x+random.randint(-5,5), y+random.randint(-5,5), z+random.randint(-5,5),  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.5)
                                
                                elev = 10.0
                                rot = 80.0 / 180 * np.pi
                                # ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.02) ## plot the sphere
                                #calculate vectors for "vertical" circle
                                a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
                                b = np.array([0, 1, 0])
                                b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
                                ax.plot(np.sin(u),np.cos(u),0,color='k', linestyle = 'dashed')
                                horiz_front = np.linspace(0, np.pi, 100)
                                ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color='k')
                                vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
                                ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed')
                                ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k')

                                #fig = plt.figure()
                                #### plot sphere and track together ######
                                plt_sphere(list_center, list_radius)
                                #plt.show()

		                print "done with reading track"
		                # exit()
		       except Exception as e:
		               print(e)

                       if countPMT == 1:
                            break

# Hspeed.Draw()
# HdistVsTime.Draw()

                for i in PMTid1:
                    H2D_Hit1mcPMT_CosTheta_Phi.Fill(Phi[i], CosTheta[i])

                for i in PMTid2:
                    H2D_Hit2mcPMT_CosTheta_Phi.Fill(Phi[i], CosTheta[i])


                #fsave.cd()
                #H2D_Hit1_CosTheta_Phi.Write("H2D_Hit1_CosTheta_Phi_evt%d"%eventID)
                #H2D_Hit2_CosTheta_Phi.Write("H2D_Hit2_CosTheta_Phi_evt%d"%eventID)
                #H2D_Hit1mcPMT_CosTheta_Phi.Write("H2D_Hit1mcPMT_CosTheta_Phi_evt%d"%eventID)
                #H2D_Hit2mcPMT_CosTheta_Phi.Write("H2D_Hit2mcPMT_CosTheta_Phi_evt%d"%eventID)
                #H2D_trigPMT_CosTheta_Phi.Write("H2D_trigPMT_CosTheta_Phi_evt%d"%eventID)
                #H2D_mcPMT_CosTheta_Phi.Write("H2D_mcPMT_CosTheta_Phi_evt%d"%eventID)
                #H2D_trigPMT_q_CosTheta_Phi.Write("H2D_trigPMT_q_CosTheta_Phi_evt%d"%eventID)
                #H2D_trigPMT_t_CosTheta_Phi.Write("H2D_trigPMT_t_CosTheta_Phi_evt%d"%eventID)


print "***************************"
print "total trigger count", countTrig
print "trigger after cuts", countTrig_afterCuts

fsave.cd()
Hspeed.Write()
Hspeed_vac.Write()
H_Ek_vac.Write()
H_deltaEk.Write()
H_Ek_tot.Write()
H_Ei.Write()
H_Ef.Write()
H_Ehit.Write()
H_energy_all.Write()
H_deltaEp_vac.Write()
H_deltaEk_vac.Write()


#c1 = TCanvas("c1","hit1",800,600)
#c1.cd()
#H2D_Hit1_CosTheta_Phi.Draw("colz")
#
#c2 = TCanvas("c2","hit2",800,600)
#c2.cd()
#H2D_Hit2_CosTheta_Phi.Draw("colz")
#
#c3 = TCanvas("c3","MC pmts",800,600)
#c3.cd()
#H2D_mcPMT_CosTheta_Phi.Draw("colz")
##H2D_Hit1mcPMT_CosTheta_Phi.Draw("colz")
#
#c4 = TCanvas("c4","trig pmts",800,600)
#c4.cd()
#H2D_trigPMT_CosTheta_Phi.Draw("colz")
#
#raw_input("Press Enter to continue...")
