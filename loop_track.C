#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
//#include <RAT/TrackNode.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNav.hh>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "TH3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
void loop_track()
{
    using namespace std;
    RAT::DU::DSReader dsReader("FitMPW1p40_tightTres_modeCut_Water_6MeV_100evt_onePosition_trackOn.root");
    int countBetaOnly =0,countNone=0, count2gamma=0, count1gamma=0, countGammaInside = 0;
    TFile *ftrack = new TFile("ftrack.root","recreate");
    TH1F *histMomentum = new TH1F("histMomentum","track momentum",100,0,100);
    TH1F *histPos = new TH1F("histPos","track position",1000,0,9000);

    for (size_t iEntry = 0; iEntry<1500;iEntry++)//dsReader.GetEntryCount();iEntry++)//
    {
     const RAT::DS::Entry& ds = dsReader.GetEntry(iEntry);
     const RAT::DS::MC& rmc= ds.GetMC();
     int particleNum = rmc.GetMCParticleCount();
       //std::cout<<"event "<<iEntry<<"  "<<particleNum<<std::endl;
       //for the 1st initial gamma

     for(size_t iev=0;iev<ds.GetEVCount(); iev++)
     {
      // Get the Event information    
      const RAT::DS::EV& rev = ds.GetEV(iev);

      if(rev.GetGTID() ==  587)
      {
 	vector<unsigned int> trackInfo = rmc.GetMCTrackIDs();
	cout<< " Number of Tracks "<< rmc.GetMCTrackCount() <<endl;
	for(vector<unsigned int>::iterator it = trackInfo.begin(); it != trackInfo.end(); ++it)
        {
	 unsigned int iTrack = *it; 
	 //cout<<rmc.GetMCTrack(iTrack).GetTrackID()<<endl;
	 if(rmc.GetMCTrack(iTrack).GetTrackID()==1)
	 {
	  cout<<"track length "<<rmc.GetMCTrack(iTrack).GetLength()<<endl;
	  TVector3 xi = rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition();
	  TVector3 xf = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition();
	  //vector<unsigned int> stepInfo = rmc.GetMCTrackStep(iTrack);
  	  int trackSteps = rmc.GetMCTrack(iTrack).GetMCTrackStepCount() ;	  
	  for(int k = 0; k<trackSteps; k++)
	  {
	   TVector3 stepPos = rmc.GetMCTrack(iTrack).GetMCTrackStep(k).GetPosition();
	   TVector3 stepMomentum = rmc.GetMCTrack(iTrack).GetMCTrackStep(k).GetMomentum();
	   cout<<stepPos.X()<<", "<<stepPos.Y()<<", "<<stepPos.Z()<<endl;
	   cout<<stepMomentum.X()<<", "<<stepMomentum.Y()<<", "<<stepMomentum.Z()<<endl;
	  }
	 }
       }
      }//eventID     
     }//triggered event
    }//loop events

    ftrack->cd();
    histMomentum->Write(); 
    histPos->Write();
    ftrack->Write();
    ftrack->Close();
}
