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
#include "TH3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
//-------RAT 6.1.1-----------//
//New track extraction !
//UNIT: cm
//PDG code==11  //electron
//PDG code==22  //gamma

using namespace std ;

void  DataAnalysisN16()
{
   RAT::DU::DSReader dsReader("RatN16_center_20000evts-5.root");
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   
   TH1D *hE_totalE = new TH1D("hE_totalE","total electron energy", 100,0,10);

   TH1D *hE_Compton = new TH1D("hE_Compton","Compton scattering", 100,0,10);
   //WavelengthEmission->GetXaxis()->SetTitle("Wavelength [nm]");

   TH1D *hE_PairProduct = new TH1D("hE_PairProduct","pair production", 100,0,10);
   //WavelengthAbsorption->GetXaxis()->SetTitle("Wavelength [nm]");
   
   TH1D *hE_PEabsorb = new TH1D("hE_PEabsorb","Photoelectric absorption", 100,0,10);
   //WavelengthAbsorption->GetXaxis()->SetTitle("Wavelength [nm]");
   
   TH1D *hPos_Cerenkov= new TH1D("hPos_Cerenkov","Cerenkov point",1000,-200,200);
   TH1D *hEbeta = new TH1D("hEbeta","initial beta Energy",100,0,100);
   TH1D *hEgamma1 = new TH1D("hEgamma1","initial gamma1 Energy",100,0,100);
   TH1D *hEgamma2 = new TH1D("hEgamma2","initial gamma2 Energy",100,0,100);
   TH1D *hEgammaFirstCompton = new TH1D("hEgammaFirstCompton","gamma first compton energy",100,0,100);
   
   TH1D *hEventPos1 = new TH1D("hEventPos1","gamma1 event position",1000,-1000,1000);
   TH1D *hEventPos2 = new TH1D("hEventPos2","gamma2 event position",1000,-1000,1000);
     
   TH1D *hTrackLength1 = new TH1D("hTrackLength1","gamma1 track Length",1000,0,6000);  
   TH1D *hTrackLength2 = new TH1D("hTrackLength2","gamma2 track Length",1000,0,6000);  
   TH1D *hPos_firstPoint = new TH1D("hPos_first","first interaction point",100,-200,200);
   TH1D *hPos_firstPointCut = new TH1D("hPos_firstCut","first interaction point with 200cm cut",1000,-200,200);
   TH1D *hPos_firstPointCut2 = new TH1D("hPos_firstCut2","parent ID ==1",100,-200,200);
   TH1D *hPos_firstPointCut3= new TH1D("hPos_firstCu3t","remove Cerenkov",100,-200,200);
   
   TVector3 PMTVector, calpmtVector, MomentumVec, eventPosition, momentumCerenkov, posCherenkov ;
   
   eventPosition.SetXYZ(0,0,0);

   double grVelocity = 2.17554021555098529e+02;//184.513002875  ;// 2.17554021555098529e+02 ;
   double eventTime = 0.0 ;
   double KECerenkov = 0.0 ;
   double KEWLS = 0.0 ;
   
   double cerenkovangle = 0.0 ;
   double CerenkovWavelength = 0.0 ;
   double WLSShiftedWavelength =0.0 ;
  
   double energy, wavelength ;	

   UInt_t UTrackID = 999999 ;
   UInt_t RTrackID = 999999 ;
   UInt_t PTrackID = 999999 ;//Parent TrackID
   UInt_t CTrackID = 999999 ;//Child TrackID
   UInt_t TrackID = 999999 ;   
   Int_t nPhotons_r =0;
   Int_t nPhotons_c =0;

   TVector3 startPos ;
   double startTime_track2 =0.0, endTime_track1=0.0 , startTime_track1= 0.0 ;
  
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
  {
     std:cout << " event ID "<< iEntry <<std::endl ;
     eventPosition.SetXYZ(0,0,0);
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     const RAT::DS::MC& rmc= rDS.GetMC();
     nPhotons_r = nPhotons_r+ rmc.GetNReemittedPhotons() ;
     nPhotons_c = nPhotons_c + rmc.GetNCherPhotons() ;
     int particleNum = rmc.GetMCParticleCount() ;
     //cout<<particleNum<<endl;
     if(particleNum>1 && rmc.GetMCParticle(1).GetPDGCode()==22){
     eventPosition= rmc.GetMCParticle(1).GetPosition();
     hEventPos1->Fill(eventPosition.X());
	 } 
	 
     if(particleNum>2 && rmc.GetMCParticle(2).GetPDGCode()==22){
     eventPosition= rmc.GetMCParticle(2).GetPosition();
     hEventPos2->Fill(eventPosition.X());
	 }   
	 
	 hEbeta->Fill(rmc.GetMCParticle(0).GetKineticEnergy());
	 
     if(particleNum>1) {hEgamma1->Fill(rmc.GetMCParticle(1).GetKineticEnergy());}//std::cout<<rmc.GetMCParticle(1).GetPDGCode()<<endl;}
     if(particleNum>2) {hEgamma2->Fill(rmc.GetMCParticle(2).GetKineticEnergy());}//std::cout<<rmc.GetMCParticle(2).GetPDGCode()<<endl;}
     
     vector<unsigned int> trackInfo = rmc.GetMCTrackIDs();
     cout<< " Number of Tracks "<< rmc.GetMCTrackCount() <<endl;
  
     //----------Electron interaction----------------------------------------------------------------  	
     for(vector<unsigned int>::iterator it = trackInfo.begin(); it != trackInfo.end(); ++it)  
     {

       if(rmc.GetMCTrack(*it).GetParticleName() == "e-" && rmc.GetMCTrack(*it).GetFirstMCTrackStep().GetProcess() =="compt")
       { 
	     hE_Compton->Fill(rmc.GetMCTrack(*it).GetFirstMCTrackStep().GetDepositedEnergy());  
       }
      }
     
 
    //-----------Gamma interaction----------------------------------------------------------------  
    //cout<<"all tracks "<<rmc.GetMCTrackCount()<<endl;
    // cut off the events with beta- only   
      for(vector<unsigned int>::iterator it = trackInfo.begin(); it != trackInfo.end(); ++it)
      {
	   unsigned int iTrack = (*it);	    
       if(particleNum>1 && rmc.GetMCParticle(1).GetPDGCode()==22) 
		 hTrackLength1->Fill(rmc.GetMCTrack(iTrack).GetLength());
	   if(particleNum>2 && rmc.GetMCParticle(2).GetPDGCode()==22) 
		 hTrackLength2->Fill(rmc.GetMCTrack(iTrack).GetLength());

	   if(rmc.GetMCTrack(iTrack).GetParticleName() == "gamma" && rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetProcess()=="compt") 
	   {	 
	     int comptonTrackID =rmc.GetMCTrack(iTrack).GetTrackID();
	     cout<<comptonTrackID<<endl;    
	    	     
	     int parentID = rmc.GetMCTrack(comptonTrackID).GetParentID();
	     
	     double gammaEnergy = rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetKineticEnergy() ;
	     hEgammaFirstCompton->Fill(gammaEnergy);
	          
	     hPos_firstPoint->Fill(rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetPosition().X()/10);  
	     
	     cout<<rmc.GetMCTrack(parentID).GetFirstMCTrackStep().GetProcess()<<endl;
	     
	     double parentEgamma = rmc.GetMCTrack(parentID).GetFirstMCTrackStep().GetKineticEnergy() ;
	      
	     if(parentEgamma>6.0 && rmc.GetMCTrack(comptonTrackID).GetFirstMCTrackStep().GetPosition().Mag()<2000)
	     {hPos_firstPointCut->Fill(rmc.GetMCTrack(comptonTrackID).GetFirstMCTrackStep().GetPosition().X()/10);}   
	     
	     if(parentID == 1)
	     {hPos_firstPointCut2->Fill(rmc.GetMCTrack(comptonTrackID).GetFirstMCTrackStep().GetPosition().X()/10);}    
	     
	     if(rmc.GetMCTrack(parentID).GetFirstMCTrackStep().GetProcess()!="Cerenkov" && rmc.GetMCTrack(parentID).GetFirstMCTrackStep().GetProcess()!="compt")
	       {hPos_firstPointCut3->Fill(rmc.GetMCTrack(comptonTrackID).GetFirstMCTrackStep().GetPosition().X()/10);}//rmc.GetMCTrack(iTrack).GetFirstMCTrackStep().GetStartVolume()=="inner_av" &    
        }//end of iTrack loop
       }
   
   }
  

 
TFile *f=new TFile("AnalyTrack_N16-5.root","RECREATE");
f->cd();  
hE_totalE->Write();hEgamma1->Write();hEgamma2->Write();hEgammaFirstCompton->Write();
hEbeta->Write();
hTrackLength1->Write();hTrackLength2->Write();
hE_Compton->Write();hEventPos1->Write();hEventPos2->Write(); 
hE_PairProduct->Write();
hE_PEabsorb->Write(); 
hPos_firstPoint->Write(); hPos_firstPointCut->Write();hPos_firstPointCut2->Write();hPos_firstPointCut3->Write();
 
}// end of the DataAnalysis_new function

