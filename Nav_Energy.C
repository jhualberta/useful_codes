/*
Rat6.4.1, Updated 9 Feb 2018

Rat6.5.2, Updated 8 Mar 2018
*/
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/TrackNode.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNav.hh>
#include <iostream>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "TH3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
//NOTE: in Boulay thesis, he calculates total energy of e-, E = sqrt(0.511^2+Ek^2) MeV
void Nav_Energy()
{
    using namespace RAT;
    //using namespace std;
    RAT::DU::DSReader dsReader("trackOn_N16_r107055_10000evt_1.root");
//  RAT::DU::DSReader dsReader("/media/vphys/seagate/N16/RatSimuData/trackOn_D2O_N16Source_center_10000evts.root");
    TH1F *hSx = new TH1F("hSx","Sx",2000,-10000,10000);
    TH1F *hSy = new TH1F("hSy","Sy",2000,-10000,10000);
    TH1F *hSz = new TH1F("hSz","Sz",2000,-10000,10000);
    TH1F *hTimeElectron = new TH1F("hTimeElectron","child e- global time",100,0,50);
    TH1F *hE_initialgamma = new TH1F("hE_initialgamma","E gamma at track start",100,0,10);
    TH1F *hE_Compton = new TH1F("hE_Compton","compton energy",100,0,10); 
    TH1F *hE_PairProduct = new TH1F("hE_PairProduct","pair production",100,0,10);
    TH1F *hE_PEabsorb = new TH1F("hE_PEabsorb","PEabsorb",100,0,10);
    double Eelectron = 0.511; 
    TVector3 pos0; pos0.SetXYZ(0,0,0);
    
    int countBetaOnly =0,countNone=0, count2gamma=0, count1gamma=0, countGammaInside = 0;
    int countCompt = 0, countProd = 0, countPhot = 0;
    unsigned int evtNum = 4;//dsReader.GetEntryCount();

    for (size_t iEntry = 3; iEntry<evtNum; iEntry++)//dsReader.GetEntryCount();iEntry++)
    {
       const RAT::DS::Entry& ds = dsReader.GetEntry(iEntry);
       const RAT::DS::MC& rmc= ds.GetMC();
       int particleNum = rmc.GetMCParticleCount();
       std::cout<<"event "<<iEntry<<"  "<<particleNum<<std::endl;
       //check for e- caused gammas
       //RAT::TrackNav nav(&ds);
       //RAT::TrackCursor c0 = nav.Cursor(false);   
       //c0.GoChild(0);//e- track
       //RAT::TrackNode *n0 = c0.Here();
       //bool flag= 1;bool checkLeak = 0;
       //while(flag)
       //{
         //TString volume = n0->GetStartVolume();
         //TString process = n0->GetProcess();
         //if(process=="Cerenkov")
         //{flag = 0;cout<<"found!"<<endl;}
         //else {
			 //if(c0.IsTrackEnd()) {checkLeak = 1 ; break;} //if hit the bottom, break in case of overflow
			 //else {c0.GoNext(); *n0 = c0.Here();}
			 //}
		//}

       //for the 1st initial gamma
       TVector3 PosFirstGamma;
       if (particleNum>1 )//if particleNum == 1, only have electrons, pass
       {
        RAT::TrackNav nav1(&ds);
	RAT::TrackCursor c = nav1.Cursor(true);//false);   
	c.GoChild(1);//1st gamma
	RAT::TrackNode *n = c.Here();
        double Egamma1 = n->GetKineticEnergy();
        if(n->GetParticleName()=="gamma")
        {
         if((6<Egamma1 && Egamma1<6.2) || (7<Egamma1 && Egamma1<7.2) )
         {
          hE_initialgamma->Fill(n->GetKineticEnergy());
	  bool flag= 1;bool checkLeak = 0;
	  while(flag)
	  {
	     TString volume = n->GetStartVolume();
             TString process = n->GetProcess();
		     
             if(volume=="inner_av" && (process =="compt" || process =="conv" ||process =="phot" ))
             {flag = 0;}//break;
             else {
        		 if(c.IsTrackEnd()) {checkLeak = 1 ;  countGammaInside++; break;} //if hit the bottom, break in case of overflow
        		 else {c.GoNext();  RAT::TrackNode *n_next = c.Here(); n = n_next;}
                  }
             //std::cout<<n->GetProcess()<<" "<<n->GetStartVolume()<<std::endl;		      
          }	   
          //std::cout<<n->GetPosition().X()<<", "<<std::endl;
          if(!checkLeak)
          {
           int childTrack =  n->fChild.size();
           double Eprevious = n->GetKineticEnergy();
           for(int iChild=0;iChild<childTrack;iChild++)
           {
               if((n->fChild[iChild])->GetParticleName()=="e-")
               { 
                //posFirstGamma.SetXYZ(n->GetPosition().X(),n->GetPosition().Y(),n->GetPosition().Z());
                if((n->fChild[iChild])->GetProcess()=="compt")//Only record Compton
                {
                  hSx->Fill(n->GetPosition().X());hSy->Fill(n->GetPosition().Y());hSz->Fill(n->GetPosition().Z());
                }
                hTimeElectron->Fill(n->GetGlobalTime());
                if(n->GetProcess()=="compt") 
                {  double Ecompt = (n->fChild[iChild])->GetKineticEnergy();
                   if(Ecompt>7) cout<<"warning  "<<Ecompt<<endl;
                   if (sqrt(Ecompt*Ecompt + Eelectron*Eelectron)>7) std::cout<<" watch out !!"<<std::endl;
                   hE_Compton->Fill(sqrt(Ecompt*Ecompt + Eelectron*Eelectron));
                   countCompt++;}
                else if(n->GetProcess()=="conv") 
                {
                  double Epairprod = (n->fChild[iChild])->GetKineticEnergy();
                  hE_PairProduct->Fill(sqrt(Epairprod*Epairprod + Eelectron*Eelectron));
                  countProd++;}
                else if(n->GetProcess()=="phot") 
                {
                  double Eabsorb = (n->fChild[iChild])->GetKineticEnergy();
                  hE_PEabsorb->Fill(sqrt(Eabsorb*Eabsorb + Eelectron*Eelectron));
                  countPhot++;}       
              }//child is electron   
           }//loop child
          }//check leak
          count1gamma++;
         }//energy cut, ensure 6 MeV, 7 MeV
        }//ensure it is gamma11
         delete n;
       }
      ////for the 2nd initial gamma
      if(particleNum>2)
      {
         RAT::TrackNav nav2(&ds);
         RAT::TrackCursor c2 = nav2.Cursor(false);   
         c2.GoChild(2);//2nd gamma
         RAT::TrackNode *n2 = c2.Here();
         double Egamma2 = n2->GetKineticEnergy();
         if(n2->GetParticleName()=="gamma")
         {  
          if((6<Egamma2 && Egamma2<6.2) || (7<Egamma2 && Egamma2<7.2))
          {
           hE_initialgamma->Fill(n2->GetKineticEnergy());

           bool flag2= 1;bool checkLeak2 = 0;
           while(flag2)
           {
                TString volume2 = n2->GetStartVolume();
                TString process2 = n2->GetProcess();
                 
                 if(volume2 =="inner_av" && (process2 =="compt" || process2 =="conv" ||process2 =="phot" ))
                 {flag2 = 0;}
                 else {
                   	if(c2.IsTrackEnd()) {checkLeak2 = 1;  countGammaInside++; break;} //if hit the bottom, break in case of overflow
            	else {c2.GoNext(); RAT::TrackNode *n2_next = c2.Here(); n2 = n2_next;}
                 }
                  //std::cout<<n->GetProcess()<<" "<<n->GetStartVolume()<<std::endl;
                  
           }   
             //  std::cout<<n2->GetPosition().X()<<", "<<std::endl;
           if(!checkLeak2)
           {
            int childTrack2 =  n2->fChild.size();
            double Eprevious = n2->GetKineticEnergy();
            for(int iChild=0;iChild<childTrack2;iChild++)
            {
               if((n2->fChild[iChild])->GetParticleName()=="e-")
               { 
                 if((n2->fChild[iChild])->GetProcess()=="compt")//Only record Compton
                 {  
                   hSx->Fill(n2->GetPosition().X());hSy->Fill(n2->GetPosition().Y());hSz->Fill(n2->GetPosition().Z());
                 }
                 hTimeElectron->Fill(n2->GetGlobalTime());
                if(n2->GetProcess()=="compt")
                {  double Ecompt = (n2->fChild[iChild])->GetKineticEnergy();
                   hE_Compton->Fill(sqrt(Ecompt*Ecompt + Eelectron*Eelectron));
                   countCompt++;}
                else if(n2->GetProcess()=="conv")
                {
                  double Epairprod = (n2->fChild[iChild])->GetKineticEnergy();
                  hE_PairProduct->Fill(sqrt(Epairprod*Epairprod + Eelectron*Eelectron));
                  countProd++;}
                else if(n2->GetProcess()=="phot")
                {
                  double Eabsorb = (n2->fChild[iChild])->GetKineticEnergy();
                  hE_PEabsorb->Fill(sqrt(Eabsorb*Eabsorb + Eelectron*Eelectron));countPhot++;}
              }
           }   
          }
           count2gamma++;
         }//energy cut, ensure 6 MeV, 7 MeV
        }//ensure it is gamma
           delete n2;  
      }
      if(particleNum==1) countBetaOnly++;
      if(particleNum==0) countNone++;
    }//for entries     
    std::cout<<"BetaOnly  "<<countBetaOnly<<" 1-gamma "<<count1gamma<<" 2-gamma "<<count2gamma<<" gamma_inside-events "<<countGammaInside<<std::endl;
    std::cout<<"compt "<<countCompt<<" pair prod "<<countProd<<" pe absorb "<<countPhot<<std::endl;
    
    TFile *ff=new TFile("fSpatial_NewRat_test_6and7MeVgamma.root","RECREATE");
    ff->cd(); hE_initialgamma->Write();hSx->Write();  hSy->Write();hSz->Write();hE_Compton->Write();hE_PairProduct->Write();hE_PEabsorb->Write();
    hTimeElectron->Write();
    //hSx->Draw();
}
