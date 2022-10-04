from ROOT import *
from rat import *
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

#			EV s	   km/s	to conver to nm
hc = 4.1357*10e-15 * 300000 * 1e5
Colors = ["red", "blue", "green", "purple"]

File = TFile("trackMC_Run30681_Alpha_vacuum_avBulk_-xCenter.root");

#T = File.Get("T");#TTree
T = File.Get("T_satCorr")
nentries = T.GetEntries();
############################################################################################################################################################################################################################
# DC 	   = RAT.DetectorConfig.GetDetectorConfig(0);#Connects to https://deapdb.physics.carleton.ca/deapdb
# nPMTs  = DC.GetNumPMTs();#Number of PMTs
# pmtPos = DC.GetPMTPositions();# A list of 255 TVector3's containing the position of the 255 PMTs
# print "\n\n\n\n\n\n\n\n************************************************* Detector Configuration loaded *************************************************\n\n\n\n\n\n\n\n"
############################################################################################################################################################################################################################

print "Started reading events..."
for event in range (nentries):
	if (event+1)%10 == 0:
		print (event+1), "Analyzed..."
	T.GetEntry(event);			  DS = T.ds;
	if DS.GetTSCount() > 0 and DS.GetEVCount() > 0:
		TS = DS.GetTS(0);EV = DS.GetEV(0);MC = DS.GetMC();

		# Save all the tracks in python dictionary, key is track ID and value is the track itself
		Track_Dictionary = {};
		for i in range(MC.GetMCTrackCount()):
			Track = MC.GetMCTrack(i);	Track_Dictionary[ Track.trackID ] = Track;
###########################################################################################################################################################################################################

		for j in range(MC.GetPMTCount()):
			pmt = MC.GetPMT(j)
			Cathodtimes = []; Photons = []
			for l in range(pmt.GetMCPhotonCount()):
				Cathodtimes.append(pmt.GetMCPhoton(l).GetCathodeTime());	Photons.append(pmt.GetMCPhoton(l));

			# Earliest Photons in each PMT!
			# This is what I wanted to do! You choos your main track you want to trace back!
			Index = Cathodtimes.index(min(Cathodtimes));		Earliest_photon = Photons[Index];
			PhotonTrackID = Earliest_photon.trackID;			ParentTrackid = PhotonTrackID;

			try:
				# List of all the parent tracks starting with the chosen track!
				ParentTracks = [Track_Dictionary[PhotonTrackID]];
				# Get the parent ID of each parent until you reach the source
				while ParentTrackid > 1:# Not the event itself
					ParentTrackid = Track_Dictionary[ParentTrackid].parentID;		ParentTracks.append(Track_Dictionary[ParentTrackid]);
				fig = plt.figure(figsize=(20,20));		ax = plt.axes(projection='3d');				Title = "";
				# ax.set_title( "T_res = " + str (T_RES) + "; CellID = " + str(CellID) )
				print "##############################################################################################################################################"
				Vectors = [];					Time = [];
				for processes in range(len(ParentTracks)):
					# We need to loop backward since we traced the track backward; i.e the first track
					# in the list is the final track!
					Track = ParentTracks[-processes-1]
					print "\n New Track: ",Track.particleName,":\n"
					# Loop over steps and print some stuff nicely!
					for steps in range(Track.GetMCTrackStepCount()):
						x = Track.GetMCTrackStep(steps).GetEndpoint().x();	y = Track.GetMCTrackStep(steps).GetEndpoint().y();	z = Track.GetMCTrackStep(steps).GetEndpoint().z();	r = np.sqrt(x**2+y**2+z**2);
						Temp = TVector3();			Temp.SetXYZ(x,y,z);
						Vectors.append(Temp);		Time.append(Track.GetMCTrackStep(steps).GetGlobalTime());

						if Track.particleName != "opticalphoton" and len(Vectors) > 1 and (Time[-1]-Time[-2]) != 0 :#	process Volume time lambda  x   y    z
							print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
							Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), "    ------",
							round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), round((Vectors[-1]-Vectors[-2]).Mag()/(Time[-1]-Time[-2]),7)
							)

						elif Track.particleName == "opticalphoton" and len(Vectors) > 1 and (Time[-1]-Time[-2]) != 0 :#	process Volume time lambda  x   y    z
							print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
							Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), round(hc/Track.GetMCTrackStep(steps).GetKE(),1),
							round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), round((Vectors[-1]-Vectors[-2]).Mag()/(Time[-1]-Time[-2]),7)
							)
						elif Track.particleName == "opticalphoton" and len(Vectors) > 1 and (Time[-1]-Time[-2]) == 0 :#	process Volume time lambda  x   y    z
							print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format(
							Track.GetMCTrackStep(steps).GetProcess(),Track.GetMCTrackStep(steps).GetVolume(), round(Track.GetMCTrackStep(steps).GetGlobalTime(),8), round(hc/Track.GetMCTrackStep(steps).GetKE(),1),
							round(x,2),round(y,2),round(z,2),round(r,2), round((Vectors[-1]-Vectors[-2]).Mag(),7), "    ------"
							)

					xdata = [];     ydata = [];     zdata = [];
					Track = ParentTracks[-processes-1]
					Title += Track.particleName + "(" + str(Track.GetMCTrackStepCount()) + ")" + "-->";
					for steps in range(Track.GetMCTrackStepCount()):
						xdata.append( Track.GetMCTrackStep(steps).GetEndpoint().x() );
						ydata.append( Track.GetMCTrackStep(steps).GetEndpoint().y() );
						zdata.append( Track.GetMCTrackStep(steps).GetEndpoint().z() );
					ax.plot3D(xdata, ydata, zdata, Colors[processes]);		Text = str(processes) + "," + Track.particleName;

					for steps in range(Track.GetMCTrackStepCount()):
						if processes > 0:
							ax.scatter( Track.GetMCTrackStep(steps).GetEndpoint().x(), Track.GetMCTrackStep(steps).GetEndpoint().y(), Track.GetMCTrackStep(steps).GetEndpoint().z() );

				print "*****************************************************************************************************************************************************************************************************"
				print "{:40} {:15} {:15} {:10} {:10} {:10} {:10} {:10} {:15} {:15}".format( "Process", "Volume", "     Time", "   Lambda", "\tX", "  Y", "  Z", "  R", "  Distance", "Speed" )
				print "*****************************************************************************************************************************************************************************************************"

				# Wait for the user to read the track and then plot it!
				raw_input("Press Enter to continue...")
				Text = str(processes) + "," + Track.particleName;			ax.text2D(0.01, 0.99, Title, transform=ax.transAxes);	ax.text2D(0.01, 0.95, round(Cathodtimes[Index],2) , transform=ax.transAxes)
				ax.set_xlabel('X [mm]');									ax.set_ylabel('Y [mm]');
				ax.set_zlabel('Z [mm]');									plt.show();
				print ""
				# exit()

			except Exception as e:
				print(e)
