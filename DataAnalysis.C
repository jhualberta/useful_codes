#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"

using namespace std ;

void  DataAnalysis()
{
   RAT::DU::DSReader dsReader("WLS_beta_5MeV.root");
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();

   TH2D* PMTIDHitTime = new TH2D("PMTIDHitTime", "PMTID vs Hit Time", 10000,1,10000,600,0,600);
   //TH2D* AngleHitTime = new TH2D("AngleHitTime", "Cherenkov angle vs time residual distribution", 100,-1,1,600,0,600);
   //TH1D* ChrenkovDis = new TH1D("ChrenkovDis","PMTPosition*MomentumVector",100,-1,1);
   //TH1D* PMThitTime = new TH1D("PMThitTime","pmt hit time distribution",400,0,400);
   TH1D *timeDiff = new TH1D("timeDiff","time difference between absorbed and reemitted photon", 1600, -100, 300);
   TH2D *timewavelength = new TH2D("timewavelength"," time difference vs wavelength of Cerenkov photon", 1000,-50,50, 50000, 0, 1000) ;
   TH1D *CerenkovAngle = new TH1D("CerenkovAngle","Cerenkov Cone angle distribution", 200,-1, 1);
   TH1D *distanceHist = new TH1D("distanceHist","distance between cerenkovv and scintillated photons", 10000, 0, 100);

//   TH1D *hWavelength = new TH1D("hWavelength","wavelength of the photons", 10000,0,1000 );
//   TH1D *hWavelengthPMT = new TH1D("hWavelengthPMT","wavelength of the photons reaching PMT", 10000,0,1000 );

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
   double energy, wavelength ;	

   UInt_t UTrackID = 999999 ;

     Int_t nPhotons_r =0;
     Int_t nPhotons_c =0;

   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
  {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     const RAT::DS::MC& rmc= rDS.GetMC();
     nPhotons_r =nPhotons_r+ rmc.GetNReemittedPhotons() ;
     nPhotons_c = nPhotons_c + rmc.GetNCherPhotons() ;
    
     if(rmc.GetMCParticle(0).GetPDGCode()==11)  MomentumVec.SetXYZ(rmc.GetMCParticle(0).GetMomentum().x(),rmc.GetMCParticle(0).GetMomentum().y(),rmc.GetMCParticle(0).GetMomentum().z()) ;
     TVector3 startPos_track2, endPos_track1, startPos_track1 ;
     double startTime_track2 =0.0, endTime_track1=0.0 , startTime_track1= 0.0 ;
     for(size_t iTrack=0; iTrack< 15; iTrack++)  // rmc.GetMCTrackCount(); iTrack++)
       {
	 if(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()==process2)//Cerenkov
          {
             momentumCerenkov.SetXYZ( rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().x(), rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().y(), rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().z()) ;
             KECerenkov = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetKineticEnergy() ;
             CerenkovWavelength = ((6.625e-34)*(3e+8)/(KECerenkov*1e+6*1.6e-19))*1e+9 ;
             cerenkovangle = (momentumCerenkov.Unit()*MomentumVec.Unit()) ;
             CerenkovAngle->Fill(cerenkovangle) ;
          }
         if(rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetProcess()==process3) 
	    {
       	        startPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().z()) ;
	        endPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().z());        
                startTime_track1 = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetGlobalTime();
                endTime_track1 = rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetGlobalTime();
	        UTrackID = rmc.GetMCTrack(iTrack).GetTrackID() ;
		for(size_t iTr=0; iTr< rmc.GetMCTrackCount(); iTr++)
		{
		  if(rmc.GetMCTrack(iTr).GetParentID() == UTrackID && rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetEndVolume()=="inner_av") 
		     {
			startPos_track2.SetXYZ(rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().z());
			distanceHist->Fill((endPos_track1-startPos_track2).Mag());
			startTime_track2 = rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetGlobalTime() ;
                	for(int i=1; i< 9727; i++ ){
	                      PMTVector = pmtInfo.GetPosition(i) ;
        	              if(PMTVector.Mag()==0 || pmtInfo.GetType(i)!=1) {
                	           double timeFill = -(PMTVector-startPos_track2).Mag()/grVelocity + (PMTVector-startPos_track1).Mag()/grVelocity + (startTime_track2 - startTime_track1) ;
                        	   timewavelength->Fill(timeFill, CerenkovWavelength);
	                           if(startTime_track2 >= endTime_track1) timeDiff->Fill(timeFill) ;
        	                }
                	   }
	    	      }
		}
             } 
	 }
      
          

     /*
     for(size_t iTrack=0; iTrack< rmc.GetMCTrackCount(); iTrack++)
       {
          startPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().z()) ;
	  endPos_track1.SetXYZ(rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetPosition().z());	  
	  startTime_track1 = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetGlobalTime();
	  endTime_track1 = rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetGlobalTime();

          if(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()==process2) {
		momentumCerenkov.SetXYZ( rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().x(), rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().y(), rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetMomentum().z()) ;
		KECerenkov = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetKineticEnergy() ;
		CerenkovWavelength = ((6.625e-34)*(3e+8)/(KECerenkov*1e+6*1.6e-19))*1e+9 ;
		cerenkovangle = (momentumCerenkov.Unit()*MomentumVec.Unit()) ;
		CerenkovAngle->Fill(cerenkovangle) ;
	    }
          //std::cout<< " process name of first track "<< rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess() << std::endl ;

          for(size_t iTr=0; iTr< rmc.GetMCTrackCount(); iTr++)	  
	    {
		startPos_track2.SetXYZ(rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().x(),rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().y(),rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetPosition().z());
	        if(rmc.GetMCTrack(iTrack).GetPDGCode()==22 && (rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()==process2 && rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetProcess()==process1) && (endPos_track1-startPos_track2).Mag()< 1.0) {
			distanceHist->Fill((endPos_track1-startPos_track2).Mag());
			std::cout << "step length of first track "<< rmc.GetMCTrack(iTrack).GetLastMCTrackStep().GetLength() << " step length of second track " << rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetLength() << "end position of the first track "<< endPos_track1.x()<<","<< endPos_track1.y()<<","<< endPos_track1.z() << " end position of the second track :"<< startPos_track2.x()<<","<< startPos_track2.y()<<","<< startPos_track2.z()<<std::endl ;
			}
		//std::cout<< " process name of second track "<< rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetProcess() <<  " distance "<< (endPos_track1-startPos_track2).Mag()<< std::endl ;
		if((endPos_track1-startPos_track2).Mag()< 30.0 && (rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()==process2 && rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetProcess()==process1))                       {
			startTime_track2 = rmc.GetMCTrack(iTr).GetFirstMCTrackStep().GetGlobalTime() ;
			for(int i=1; i< 9727; i++ ){
    			    PMTVector = pmtInfo.GetPosition(i) ;
			    if(PMTVector.Mag()==0 || pmtInfo.GetType(i)!=1) {
					double timeFill = -(PMTVector-startPos_track2).Mag()/grVelocity + (PMTVector-startPos_track1).Mag()/grVelocity + (startTime_track2 - startTime_track1) ;
					timewavelength->Fill(timeFill, CerenkovWavelength);
					if(startTime_track2 >= endTime_track1) timeDiff->Fill(timeFill) ;
				}
		            }
	                 }
                  }
            }
          */
      /*
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
            {
              const RAT::DS::CalPMTs& calibratedPMTs = rDS.GetEV( iEV ).GetCalPMTs();
	      for( size_t iPMT = 0; iPMT < calibratedPMTs.GetAllCount(); iPMT++ )
                {
                      const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetAllPMT( iPMT );
                      int PMTID = pmtCal.GetID() ;
                      PMTVector = pmtInfo.GetPosition(PMTID) ;
		
                      double TOF = (PMTVector - eventPosition).Mag()/grVelocity ;     

		      double angle = MomentumVec*PMTVector.Unit();
                      double hitTime = pmtCal.GetTime();
        	      double timeRes = hitTime - TOF - eventTime ;

		      //PMThitTime->Fill(hitTime);
                      PMTIDHitTime->Fill(PMTID, hitTime);
		      //ChrenkovDis->Fill(angle);		
		      //AngleHitTime->Fill(angle,hitTime );

		}
           }
    */

   }


  TFile *f=new TFile("WLS_process_5MeV.root","RECREATE");
  f->cd();  
  timeDiff->Write();
  CerenkovAngle->Write();
  timewavelength->Write();
  distanceHist->Write();
//  hWavelength->Write();
//  hWavelengthPMT->Write();	
  //AngleHitTime->Write();
  //ChrenkovDis->Write();
  //PMThitTime->Write();
  PMTIDHitTime->Write();
  
//  f->Close();

std::cout<< " Cerenkov photons "<< nPhotons_c << " Remitted photons " << nPhotons_r <<std::endl ;

}
