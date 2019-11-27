///*Extract 1st gamma interaction point, using optizmized method*/
///13 Mar 2017
//#include <RAT/DU/DSReader.hh>
//#include <RAT/DU/Utility.hh>
//#include <RAT/DS/Entry.hh>
//#include <RAT/DS/MC.hh>
//#include <RAT/DS/EV.hh>
//#include <RAT/DS/PMT.hh>
//#include <RAT/DU/PMTInfo.hh>
//#include <RAT/TrackNode.hh>
//#include <RAT/TrackCursor.hh>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "TH3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include <iostream>
#include <fstream>
void Nav_New1()
{   
    double ECh_thresh = 0.26577;//[MeV] e- Cherenkov thresh. in heavy water 
    double grVelocity = 218.255; //Effective overall speed of light //n = 1.328, c/n=225.75
    //using namespace std;
    RAT::DU::DSReader dsReader("/media/vphys/seagate/N16/RatSimuData/trackOn_D2O_N16Source_center_10000evts.root");
    const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
    TH1F *hSx = new TH1F("hSx","Sx",1200,-6000,6000);
    TH1F *hSy = new TH1F("hSy","Sy",1200,-6000,6000);
    TH1F *hSz = new TH1F("hSz","Sz",1200,-6000,6000);
    
    //all Cherenkov steps
    TH1F *hSx_all = new TH1F("hSx_all","Sx all",1200,-6000,6000);
    
    TH1F *hTimeElectron = new TH1F("hTimeElectron","child e- global time",1000,0,100);
    TH1F *hE_Compton = new TH1F("hE_Compton","compton energy",100,0,10); 
    TH1F *hE_PairProduct = new TH1F("hE_PairProduct","pair production",100,0,10);
    TH1F *hE_PEabsorb = new TH1F("hE_PEabsorb","PEabsorb",100,0,10);
    //for each step gamma transfers energy to electrons. E> 0.26577 Cherenkov threshold
    TH1F *hE_scatter_e = new TH1F("hE_scatter_e","gamma scattered electron energy (>Ech threshold)",100,0,10); 
    
    TH1F *htime_differ = new TH1F("htime_resol","time resolution",2000,-100,100);
    TH1F *hToF_differ = new TH1F("hToF_differ","time of flight differ",2000,-100,100);
    
    TH1F *hGamma1startPos = new TH1F("hGamma1startPos","1st gamma pos",2000,-1000,1000);
    TH1F *hGamma2startPos = new TH1F("hGamma2startPos","2nd gamma pos",2000,-1000,1000);
    TH1F *hCerenkovTime = new TH1F("hCerenkovTime","cerenkov time",1000,0,100);
    TH1F *hCerenkovPos = new TH1F("hCerenkovPos","cerenkov pos",2000,-1000,1000);
    TH1F *hTimeRes_firstInter = new TH1F("hTimeRes_firstInter", "time Res of first Gamma interaction", 2000,-300,300);
    TH1F *hTimeRes_source = new TH1F("hTimeRes_source", "time Res of source(exactly center) event", 2000,-300,300);
    TH1F *hTacTime = new TH1F("hTacTime","Tac time distr. for N16 event",200,0,200);
    TH1F *hEventTime = new TH1F("hEventTime","Cherenkov happen time",1000,0,100);

    
    TVector3 PMTVector; TVector3 PosCenter; PosCenter.SetXYZ(0,0,0);
    TVector3 posGamma1, posGamma2;
    TVector3 startPos_gamma1, startPos_gamma2;
    
    TH1F *hSumOfTres = new TH1F("hSumOfTres","sum of tRes for every hit PMTs",1000,0,1000);
    TH2F *hPosVsTres = new TH2F("hPosVsTres","Tres vs PosX",1000,0,1000,2000,-2000,2000);
    TH1F *hSx_min = new TH1F("hSx_min","hSx by the min value of tEvent+ToF",1200,-6000,6000);
    TH1F *hSx_optimize = new TH1F("hSx_optimize","hSx by minimize tEvent+ToF",1200,-6000,6000);
    
    
  //  startPos_gamma1.SetXYZ(0,0,0);startPos_gamma2.SetXYZ(0,0,0);
    double startTime_gamma1, startTime_gamma2;//always = 0
 
    double delta_i = 0;
    vector<double> sumOfTres_store;
    
    int countBetaOnly =0,countNone=0, count2gamma=0, count1gamma=0, countGammaInside = 0, count6MeVGamma=0, count7MevGamma=0;
    int countCompt = 0, countProd = 0, countPhot = 0;
    
    vector<int> eventIDtemp;
    vector<double> hitPMTposTemp;
    vector<double> hitTimeTemp; 
    vector<int> pmtNumberRecord; 
    vector<double> I0temp;//store first photon step
    vector<double> eventTime0temp;

    for (int iEntry =0; iEntry<1; iEntry++)//dsReader.GetEntryCount();iEntry++)
    {
       int ITER = 0;
       TGraph *gMin = new TGraph();//To calculate the minimum tRes
       const RAT::DS::Entry& ds = dsReader.GetEntry(iEntry);
       const RAT::DS::MC& rmc= ds.GetMC();
       int particleNum = rmc.GetMCParticleCount();
       std::cout<<"event "<<iEntry<<"  "<<particleNum<<std::endl;
       PMTVector.SetXYZ(0,0,0);
       //for the 1st initial gamma
       TVector3 PosFirstGamma;
       if(ds.GetEVCount()!=0)///first check whether this simulation caused an event 
       { 
        //get hitted PMT informations
        const RAT::DS::CalPMTs& calpmts = ds.GetEV(0).GetCalPMTs();
        int eventID = iEntry;
        int pmtNumber = 0;     
        for( size_t ipmt = 0; ipmt < calpmts.GetNormalCount();ipmt++)
		{
          const RAT::DS::PMTCal& pmtCal = calpmts.GetNormalPMT(ipmt);
          int PMTID = pmtCal.GetID() ;
          PMTVector = pmtInfo.GetPosition(PMTID);
          hitPMTposTemp.push_back(PMTVector.X());hitPMTposTemp.push_back(PMTVector.Y());hitPMTposTemp.push_back(PMTVector.Z());
          double hitTime = (calpmts.GetPMT(ipmt)).GetTime();
		  hitTimeTemp.push_back(hitTime);
		  pmtNumber++;
		}

		pmtNumberRecord.push_back(pmtNumber);///for tagging the PMTvector and hitTime
        //hitPMTposTemp.clear(); vector<double>().swap(hitPMTposTemp); 
              
        if (particleNum>1 )//some events only have electron, remove these events
        {
           //TH2F *hTest = new TH2F("hTest","Tres vs PosX",1000,0,500,2000,-2000,2000);	
           count1gamma++;
           RAT::TrackNav nav1(&ds);
	       RAT::TrackCursor c = nav1.Cursor(false);   
	       c.GoChild(1);//1st gamma
	       RAT::TrackNode *n = c.Here();
	       startPos_gamma1.SetXYZ(n->GetPosition().X(),n->GetPosition().Y(),n->GetPosition().Z());
	       hGamma1startPos->Fill(startPos_gamma1.X());
	       startTime_gamma1 = n->GetGlobalTime();  
	       double initialE = n->GetKineticEnergy();
	       if(initialE>6 && initialE<6.3)
	       { count6MeVGamma++;}//count 6.1 MeV gamma
	       else if(initialE>7 && initialE<7.5){count7MevGamma++;}//count 7.1 MeV gamma
	       //cout<<initialE<<" MeV "<<endl;
	       bool flag= 1;bool checkLeak = 0;
	       ///mark the interaction point c, *n->c
	       while(flag)
	       {
		     TString volume = n->GetStartVolume();
		     TString process = n->GetProcess();	     
		     if(volume=="inner_av" && (process =="compt" || process =="conv" ||process =="phot" ))
		     {flag = 0;break;}
		     else {
				 if(c.IsTrackEnd()) {checkLeak = 1; countGammaInside++; break;} //if hit the bottom, break in case of overflow
				 else {c.GoNext(); *n = c.Here();}
				 }	      
		   }	   
		  ///extract node information from *n, check no overflow
		  if(!checkLeak)//Only compton process has multiple steps
		  {
		   //the beginning of the step
		    double position0 = n->GetPosition().X();
		    hSx->Fill(n->GetPosition().X());
		    hSx_all->Fill(n->GetPosition().X());///for hSx_all, first store the first step	
		    //********************************//

		    I0temp.push_back(n->GetPosition().X());I0temp.push_back(n->GetPosition().Y());I0temp.push_back(n->GetPosition().Z());
		    eventTime0temp.push_back(n->GetGlobalTime());     
		    eventIDtemp.push_back(eventID);      
		    ///go through only the Compton steps which create Cherenkov             
           } 
           
          ///if no multiple Compton-scattering Cherenkov electrons, save the initial point
          delete n;   
		     
         ///if 2nd gamma exists, for the 2nd initial gamma
	      if(particleNum>2)
	      {          
		       count2gamma++;
			   RAT::TrackNav nav2(&ds);
		       RAT::TrackCursor c2 = nav2.Cursor(false);   
		       c2.GoChild(2);//2nd gamma
		       RAT::TrackNode *n2 = c2.Here();
		       startPos_gamma2.SetXYZ(n2->GetPosition().X(),n2->GetPosition().Y(),n2->GetPosition().Z());
	           hGamma2startPos->Fill(startPos_gamma2.X());
	           startTime_gamma2 = n2->GetGlobalTime();  
	           double initialE = n2->GetKineticEnergy();
	           if(initialE>7 && initialE<7.5) count7MevGamma++;
	           //cout<<"2nd Gamma Energy "<<initialE<<" MeV"<<endl;
		       bool flag2= 1;bool checkLeak2 = 0;
		       while(flag2)
		       {
			     TString volume2 = n2->GetStartVolume();
			     TString process2 = n2->GetProcess();
			     
			     if(volume2 =="inner_av" && (process2 =="compt" || process2 =="conv" ||process2 =="phot" ))
			     {flag2 = 0;}
			     else {
					 if(c2.IsTrackEnd()) {checkLeak2 = 1;  countGammaInside++; break;} //if hit the bottom, break in case of overflow
					 else {c2.GoNext(); *n2 = c2.Here();}
					 }      
			   }
   
			 //  std::cout<<n2->GetPosition().X()<<", "<<std::endl;
	      if(!checkLeak2)//Only compton process has multiple steps
		  {
		   //the beginning of the step
		   //last step of gamma energy
		    double position0_gamma2 = n2->GetPosition().X();
		    hSx->Fill(n2->GetPosition().X());
		    hSx_all->Fill(n2->GetPosition().X());//for hSx_all, store the first step
		    I0temp.push_back(n2->GetPosition().X());I0temp.push_back(n2->GetPosition().Y());I0temp.push_back(n2->GetPosition().Z());
		    eventTime0temp.push_back(n2->GetGlobalTime());
		    eventIDtemp.push_back(eventID);            
		   }    
		   
         delete n2;   
		} 
	  }   
      if(particleNum==1) countBetaOnly++;
      if(particleNum==0) countNone++;     
      
     }///check event exists       
    }///loop for iEntry 
    
    ///save to TTree
   TFile *ftree = new TFile("treeNavTrack.root", "update");
   TTree *tree = new TTree("navTrack", "nav_track Info");
   TVector3 eventPos, hitPMTpos;
   int pmtNumber; int ev; double pmtHitTime,eventTime; 
   tree->Branch("eventPos",&eventPos);
   tree->Branch("eventTime",&eventTime);
   tree->Branch("hitPMTpos",&hitPMTpos);
   tree->Branch("pmtNumber",&pmtNumber);
   tree->Branch("pmtHitTime",&pmtHitTime);
   tree->Branch("eventID",&ev);
   for(std::vector<double>::iterator it = eventTime0temp.begin() ; it != eventTime0temp.end(); ++it)
   {
     eventTime=(*it);
     tree->Fill();  
   } 
   
   for(std::vector<int>::iterator itb = eventIDtemp.begin() ; itb != eventIDtemp.end(); ++itb)
   {
     ev=(*itb);
     tree->Fill();  
   } 
    
   
   for (std::vector<double>::iterator it = I0temp.begin() ; it != I0temp.end(); ++it)
   {  
       double x = *it; double y = *(it+1); double z = *(it+2);
       eventPos.SetXYZ(x,y,z);
       tree->Fill();  
    } 

   ////how many hitted PMTs for one event
   for (std::vector<int>::iterator iit = pmtNumberRecord.begin() ; iit != pmtNumberRecord.end(); ++iit)
   {
     pmtNumber=(*iit);
     tree->Fill(); 
   }  

   ////PMT hit times for one event
   int count = 0; 
   for (std::vector<double>::iterator it = hitTimeTemp.begin() ; it != hitTimeTemp.end(); ++it)
   {
		 pmtHitTime = (*it); 
		 tree->Fill(); 
   }
 
   int count = 0; 
   for (std::vector<double>::iterator it = hitPMTposTemp.begin() ; it != hitPMTposTemp.end(); ++it)
   {
       double x = *it; double y = *(it+1); double z = *(it+2);
       hitPMTpos.SetXYZ(x,y,z);
	   tree->Fill(); 
   }
    
  ftree->cd();
  tree->Write();
  ftree->Close();
  
   //ofstream fEventTime, fEventPos,fpmtNumber,fHitTime,fPMTpos;
   //fEventTime.open ("fEventTime.txt"); 
   //fEventPos.open ("fEventPos.txt"); 
   //fpmtNumber.open ("fpmtNumber.txt"); 
   //fHitTime.open ("fHitTime.txt"); 
   //fPMTpos.open ("fPMTpos.txt"); 
    
   ////event time of first interaction position
   //for (std::vector<double>::iterator it = eventTime0temp.begin() ; it != eventTime0temp.end(); ++it)
   //{
     //fEventTime<<(*it)<<' ';
   //}

   ////3-D vector of first interaction position
   //for (std::vector<double>::iterator it = I0temp.begin() ; it != I0temp.end(); ++it)
   //{
     //fEventPos<<(*it)<<' ';
   //} 

   ////how many hitted PMTs for one event
   //for (std::vector<int>::iterator iit = pmtNumberRecord.begin() ; iit != pmtNumberRecord.end(); ++iit)
   //{
     //fpmtNumber<<(*iit)<<' ';
   //}  

   ////PMT hit times for one event
   //int count = 0; 
   //for (std::vector<int>::iterator iit = pmtNumberRecord.begin() ; iit != pmtNumberRecord.end(); ++iit)
   //{
      //for(int iter=0;iter<(*iit);iter++)//PMT number
      //{
		 //fHitTime<<hitTimeTemp[iter+count]<<' ';
		 //count = (*iit);  
	  //}
   //}
 
   //count = 0; 
   //for (std::vector<int>::iterator iit = pmtNumberRecord.begin() ; iit != pmtNumberRecord.end(); ++iit)
   //{
      //for(int iter=0;iter<(*iit)*3;iter++)//PMT number
      //{
		 //cout<<(*iit)<<endl;
		 //fPMTpos<<hitPMTposTemp[iter+count*3]<<' ';//hit PMT positions
		 //count = (*iit);  
	  //}
   //}
    
    //fEventTime.close();fEventPos.close();fpmtNumber.close();fHitTime.close();fPMTpos.close();
    
    std::cout<<"BetaOnly  "<<countBetaOnly<<" 1-gamma "<<count1gamma<<" 2-gamma "<<count2gamma<<" gamma inside events "<<countGammaInside<<std::endl;
    std::cout<<"compt "<<countCompt<<" pair prod "<<countProd<<" pe absorb "<<countPhot<<std::endl;
    std::cout<<"6MeV gamma "<<count6MeVGamma<<" 7MeV gamma "<<count7MevGamma<<endl;
    TFile *ff=new TFile("fSpatial_timeDiff_Optimize_test.root","RECREATE");
    ff->cd(); hSx->Write();hSx_all->Write();  hSy->Write();hSz->Write();hE_Compton->Write();hE_PairProduct->Write();hE_PEabsorb->Write();
    hTimeElectron->Write();hGamma1startPos->Write();hGamma2startPos->Write(); hEventTime->Write();
    hCerenkovPos->Write();hCerenkovTime->Write();hTimeRes_source->Write(); hTimeRes_firstInter->Write(); hTacTime->Write();
    hSumOfTres->Write();hPosVsTres->Write();hSx_min->Write();
    
    hitPMTposTemp.clear(); vector<double>().swap(hitPMTposTemp); 
    hitTimeTemp.clear(); vector<double>().swap(hitTimeTemp); 
    pmtNumberRecord.clear(); vector<int>().swap(pmtNumberRecord); 
    I0temp.clear(); vector<double>().swap(I0temp); 
    eventTime0temp.clear(); vector<double>().swap(eventTime0temp); 

}
//2000events, BetaOnly  552 1 gamma 1311 2 gamma 13 gamma inside events 138
