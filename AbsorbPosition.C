#ifndef __CINT__
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <TCanvas.h>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#endif
using namespace std ;

void  AbsorbPosition()
{
   RAT::DU::DSReader dsReader("WLS_beta_5MeV.root");
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();

   TH1D *timeDiff = new TH1D("timeDiff","time difference between absorbed and reemitted photon", 1600, -100, 300);
   //TH1D *hxAbsPos = new TH1D("hxAbsPos","x position where cherenkov light is absorbed", 2000, -6000, 6000);  
   //TH1D *hyAbsPos = new TH1D("hyAbsPos","y position where cherenkov light is absorbed", 2000, -6000, 6000);  
   //TH1D *hzAbsPos = new TH1D("hzAbsPos","z position where cherenkov light is absorbed", 2000, -6000, 6000); 
     
   TH1D *AbsDistance = new TH1D("AbsDistance","distance between cherenkov light starts and absorbed", 1000, 0, 6000);
   TH1D *CosTheta = new TH1D("cosTheta"," ", 1000, 0, 1.0);
 
   TH2D *DistanceCosTheta = new TH2D("AbsDistance-cosTheta","CosTheta",1000,0,6000,1000,0,1.0);

   TH2D *hXYzAbsPosbin1 = new TH2D("sqrt(x^2+y^2) vs z position where cherenkov light is absorbed","z",1000, 0, 6000, 2000, -6000, 6000);
   TH2D *hYZxAbsPosbin1 = new TH2D("yzXpos_bi_Xrange","x",1000, 0, 6000, 2000, -6000, 6000);
   
   TH2D *hXYzAbsPosbin2 = new TH2D("sqrt(x^2+y^2) vs z position where cherenkov light is absorbed","z",1000, 0, 6000, 1000, 0, 6000);
   TH2D *hYZxAbsPosbin2 = new TH2D("yzXpos_single_Xrange","x",1000, 0, 6000, 1000, 0, 6000);
   
   TVector3 PMTVector, calpmtVector, MomentumVec, eventPosition, momentumCerenkov ;
   
   eventPosition.SetXYZ(0,0,0);

   double grVelocity = 2.17554021555098529e+02 ;
   double eventTime = 0.0 ;
   double KECerenkov = 0.0 ;
   double cerenkovangle = 0.0 ;
   double CerenkovWavelength = 0.0 ;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   string process4 = "Reemission_from_comp1";
   
   double energy, wavelength ;	

   UInt_t UTrackID = 999999 ;
   Int_t nPhotons_r =0;
   Int_t nPhotons_c =0;
   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )//dsReader.GetEntryCount()
   {
     const RAT::DS::Entry& rDS = dsReader.GetEntry(iEntry);
     const RAT::DS::MC& rmc= rDS.GetMC();
     nPhotons_r =nPhotons_r+ rmc.GetNReemittedPhotons() ;//the number of reemitted photons
     nPhotons_c = nPhotons_c + rmc.GetNCherPhotons() ;// the number of Cherenkov photons
     
     TVector3 startPos_track2, endPos_track1, startPos_track1 ;
     double startTime_track2 =0.0, endTime_track1=0.0 , startTime_track1= 0.0 ;
     //cout<<rmc.GetMCTrackCount()<<endl;
     for(size_t iTrack=0; iTrack< rmc.GetMCTrackCount(); iTrack++)
     {       
       if(rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetProcess()==process3 &&  rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()==process2) 
       {
       	    startPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().z()) ;
	        endPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().z());        
	        UTrackID = rmc.GetMCTrack(iTrack).GetTrackID();
	                //hxAbsPos->Fill(endPos_track1.X());
			//hyAbsPos->Fill(endPos_track1.Y());
			//hzAbsPos->Fill(endPos_track1.Z());
			
			double distance;
			
    			distance=(endPos_track1-startPos_track1).Mag();		

			AbsDistance->Fill(distance);
		        
 			CosTheta->Fill(endPos_track1.X()/endPos_track1.Mag());

                        DistanceCosTheta->Fill(distance,endPos_track1.X()/endPos_track1.Mag());				
			hXYzAbsPosbin1->Fill(TMath::Sqrt(pow(endPos_track1.X(),2)+pow(endPos_track1.Y(),2)),endPos_track1.Z()); 
			hYZxAbsPosbin1->Fill(TMath::Sqrt(pow(endPos_track1.Z(),2)+pow(endPos_track1.Y(),2)),endPos_track1.X()); 
			  
			hXYzAbsPosbin2->Fill(TMath::Sqrt(pow(endPos_track1.X(),2)+pow(endPos_track1.Y(),2)),endPos_track1.Z()); 
			hYZxAbsPosbin2->Fill(TMath::Sqrt(pow(endPos_track1.Z(),2)+pow(endPos_track1.Y(),2)),endPos_track1.X()); 
	    }
     } 
   }



  TFile *f1=new TFile("histo_absorb_position_11.root","RECREATE");
  f1->cd();  
  //hxAbsPos->Write();
  //hyAbsPos->Write();
  //hzAbsPos->Write();

  hXYzAbsPosbin1->Write();
  hYZxAbsPosbin1->Write();

  AbsDistance->Write();
  CosTheta->Write();  

  DistanceCosTheta->Write();

  hXYzAbsPosbin2->Write();
  hYZxAbsPosbin2->Write();
  f1->Close();

std::cout<< " Cerenkov photons "<< nPhotons_c << std::endl;
std::cout<< " Remitted photons " << nPhotons_r <<std::endl ;

}
