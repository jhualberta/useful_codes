//2017-7-13
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <vector>
#include "TH2.h"
#include "TH1.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include <TRotation.h>
#include <numeric>
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = MultiWater
fit9 = FTP
*/
double PMTsep = 500;//50 cm
using namespace std ;
void analyCalTime()
{  
   double correctParam[]={1,1,1,1,0};//0.9838,-63.375};
   double driveP0 = correctParam[3];
   double driveP1 = correctParam[4];
   TString filename = "FitMulti_TimeCor_Analysis_r0000100934_s000_p000.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader("FitMulti_TimeCor_Analysis_r0000100934_s000_p000.root");
   TVector3 sourcePos; //!!Note: use the position get from reflection fit
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TH1F* hPMTpos = new TH1F("hPMTpos","PMT position",1000,0,9000);
   TH1F* hTheta = new TH1F("hcosThetaSNOP_MultiWater", "SNO+ MPW (fitPos-sourcePos)*fitDirec, with corrected Pos and tResCut", 1000,-1,1);
   TH1F* hThetaCut1 = new TH1F("hcosThetaCut1SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.5 m", 1000,-1,1);
   TH1F* hThetaCut2 = new TH1F("hcosThetaCut2SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.0 m", 1000,-1,1);
   TH1F* hThetaCut3 = new TH1F("hcosThetaCut3SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.75m", 1000,-1,1);
   TH1F* hThetaCut4 = new TH1F("hcosThetaCut4SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.5 m", 1000,-1,1);
   TH1F* hThetaCut5 = new TH1F("hcosThetaCut5SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.2 m", 1000,-1,1);
   TH2F* hPMTphiTheta = new TH2F("hPMTphiTheta","hitted PMT positions theta vs phi",628,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);
   TH2F* hPMTphiTheta_inCone = new TH2F("hPMTphiTheta_inCone","hitted PMT positions theta vs phi, in prompt light cone",628,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);
   TH2F* hPMTphiTheta_oppCone = new TH2F("hPMTphiTheta_oppCone","hitted PMT positions theta vs phi, opposite prompt light cone",628,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);
 
   TH1F* hXcor = new TH1F("hXcor","Xcor",2000,-9000,9000);
   TH1F* hYcor = new TH1F("hYcor","Ycor",2000,-9000,9000);
   TH1F* hZcor = new TH1F("hZcor","Zcor",2000,-9000,9000);
//drive cor with trig cut
   TH1F* hXcor_trig = new TH1F("hXcor_trig","Xcor trigger",2000,-9000,9000);
   TH1F* hYcor_trig = new TH1F("hYcor_trig","Ycor trigger",2000,-9000,9000);
   TH1F* hZcor_trig = new TH1F("hZcor_trig","Zcor trigger",2000,-9000,9000);
   TH1F* hPosMagCor_trig = new TH1F("hPosMagCor_trig","SNO+ H_{2}O MultiWater corrected fitter Mag(), trigger",1000,0,10000);
   TH1F* hDeltaXcor_trig = new TH1F("hDeltaXcor_trig","Xcor-Xtrue, trigger",2000,-9000,9000);
   TH1F* hDeltaYcor_trig = new TH1F("hDeltaYcor_trig","Ycor-Ytrue, trigger",2000,-9000,9000);
   TH1F* hDeltaZcor_trig = new TH1F("hDeltaZcor_trig","Zcor-Ztrue, trigger",2000,-9000,9000);
   TH1F* hDeltaRcor_trig = new TH1F("hDeltaRcor_trig","Rcor-Rtrue, trigger",2000,-9000,9000);
   TH2F* hXZcor_trig = new TH2F("hXZcor_trig","Xcor vs Zcor trigger",2000,-9000,9000,2000,-9000,9000);
   TH2F* hRZcor_trig = new TH2F("hRZcor_trig","Rcor vs Zcor trigger",1000,0,9000,2000,-9000,9000);

   TH2F* hPosXvsTime = new TH2F("hPosXVsTime","hit time vs posX",100,0,500,2000,-6000,6000);
   TH1F *hMfit = new TH1F("hMfit","Mfit: pos_fit to PSUP sphere",1000,1000,15000);  
   TH1F* htRes = new TH1F("htRes", "Time Residue", 400,-100,300);
   TH1F* htRes_trig = new TH1F("htRes_trig", "Time Residue, trig", 400,-100,300);
      
   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O MultiWater fitter X", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O MultiWater fitter Y", 2000,-9000,9000);   
   TH2F* hfitXY = new TH2F("hfitXY", "SNO+ H_{2}O MultiWater fitter X vs Y", 2000,-9000,9000,2000,-9000,9000); 
   TH2F* hfitXZ = new TH2F("hfitXZ", "SNO+ H_{2}O MultiWater fitter X vs Z", 2000,-9000,9000,2000,-9000,9000);

   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);
   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "fitted time,FECD cuts, triggered", 800,0,800); 
   TH2F* hfitRvsZ = new TH2F("hfitRvsZ","SNO+ H_{2}O MultiWater fitter R vs Z",1000,0,9000,2000,-9000,9000);   
   TH2F* hFlatVsZ = new TH2F("hFlatVsZ","fitted (x^2+y^2)/6000^2 vs Z",5000,0,5,2000,-9000,9000);
   TH1F* hPosMag = new TH1F("hPosMag","SNO+ H_{2}O MultiWater fitter Mag()",1000,0,10000);
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O MultiWater fitter Z", 2000,-9000,9000);

   TH2F *htResVsPMTposx = new TH2F("htResVsPMTposx","htRes vs PMTposX, trig",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTposy = new TH2F("htResVsPMTposy","htRes vs PMTposY",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTposz = new TH2F("htResVsPMTposz","htRes vs PMTposZ",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTphi = new TH2F("htResVsPMTphi","htRes vs PMTphi",628,-TMath::Pi(),TMath::Pi(),400,-100,300);
   TH2F *htResVsPMTtheta = new TH2F("htResVsPMTtheta","htRes vs PMTtheta",1000,-1,1,400,-100,300);
   
   TH1F *htRes_inPMTsep = new TH1F("htRes_inPMTsep","tRes within 50 cm",400,-100,300);
   TH1F *htRes_oppPMTsep = new TH1F("htRes_oppPMTsep","tRes out 50 cm",400,-100,300);
   TH1F *htPrompt_inPMTsep = new TH1F("htPrompt_inPMTsep","",400,-100,300);
   TH1F *htPrompt_oppPMTsep = new TH1F("htPrompt_oppPMTsep","",400,-100,300);
   
   TH2F *hSumTOFvsPMTposx = new TH2F("hSumTOFvsPMTposx","sum of TOF vs PMTposX",1800,-9000,9000,200,-100,100); 
   TH2F *hSumTOFvsPMTposy = new TH2F("hSumTOFvsPMTposy","sum of TOF vs PMTposY",1800,-9000,9000,200,-100,100); 
   TH2F *hSumTOFvsPMTposz = new TH2F("hSumTOFvsPMTposz","sum of TOF vs PMTposZ",1800,-9000,9000,200,-100,100);
     
   TH2F *hDiffTOFvsPMTposx = new TH2F("hDiffTOFvsPMTposx","difference of TOF vs PMTposX",1800,-9000,9000,500,0,500); 
   TH2F *hDiffTOFvsPMTposy = new TH2F("hDiffTOFvsPMTposy","difference of TOF vs PMTposY",1800,-9000,9000,500,-100,300); 
   TH2F *hDiffTOFvsPMTposz = new TH2F("hDiffTOFvsPMTposz","difference of TOF vs PMTposZ",1800,-9000,9000,500,-100,300);  
  
   TH2F *htResVsPMTposx_trig = new TH2F("htResVsPMTposx_trig","htRes vs PMTposX, trig",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTposy_trig = new TH2F("htResVsPMTposy_trig","htRes vs PMTposY, trig",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTposz_trig = new TH2F("htResVsPMTposz_trig","htRes vs PMTposZ, trig",1800,-9000,9000,400,-100,300);
   TH2F *htResVsPMTphi_trig = new TH2F("htResVsPMTphi_trig","htRes vs PMTphi, trig",628,-TMath::Pi(),TMath::Pi(),400,-100,300);
   TH2F *htResVsPMTtheta_trig = new TH2F("htResVsPMTtheta_trig","htRes vs PMTtheta, trig",1000,-1,1,400,-100,300);

   TVector3 u_fit,pos_fit;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02;//2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02 ; MPW in h2o: 3.0/1.383=216.77
   Double_t rPSUP = 8890;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   Double_t energy, wavelength ;
 
   Double_t countFitValid = 0;//total fitted events
   Double_t countFECDtotal =0;//total FECD==9188 events
   Double_t countSuccess = 0;
   Double_t countTimeWindow = 0;//fitValid && FECD && tFit cuts
   Double_t countClean = 0;
   Double_t countTrig = 0;
   Double_t countTrigFECD = 0;
   Double_t countNoCleanTotal=0; Double_t countNoCleanTrig=0;Double_t countNoCleanFECD=0;
   int trigWord = 2;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   int nhitCut = 0;// 15;
   unsigned int fecdID = 9188;
   TString specFit = "MultiPathProcessor";
   //TString specFit = "fit10";
   ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   
   for( size_t iEntry = 0; iEntry <dsReader.GetEntryCount();iEntry++ )// dsReader.GetEntryCount(); iEntry++ )//46814
   {
     //std:cout << " event ID "<< iEntry <<std::endl ;
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    Int_t nevC = rDS.GetEVCount();
    for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     bool specFitValid = false;//check specific fitter works
     double nhits = rev.GetNhits();
     int trig = rev.GetTrigType();//get triggerType
     for (it=fitname.begin(); it<fitname.end(); it++)
     {if(*it==specFit) specFitValid = true;}
 
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     countNoCleanTotal++;
     if(trig&(trigWord)){countNoCleanTrig++;}

     const RAT::DS::DataQCFlags& dataQCFlags = rev.GetDataCleaningFlags();
     Int_t passNumber = dataQCFlags.GetLatestPass();
     const RAT::DS::BitMask& bitMaskApplied = dataQCFlags.GetApplied(passNumber);
     const RAT::DS::BitMask& bitMaskFlags = dataQCFlags.GetFlags(passNumber);

     if( RAT::EventIsClean( rev, dcAnalysisWord ) )
     {
      countClean++;//count cleaned total events
      if(trig&(2)){countTrig++;}
      if(specFitValid){
       RAT::DS::FitVertex fitVertex =rev.GetFitResult("MultiPathProcessor").GetVertex(0);
       std::string FOMname ="FitterWater";// (rev.GetFitResult("MultiPathProcessor").GetFOMNames())[0];
       double fomValue =  rev.GetFitResult("MultiPathProcessor").GetFOM(FOMname);
       
       pos_fit=fitVertex.GetPosition(); 
       if(pos_fit.Mag()<9000){  //fit valid: inside cavity, specified by MPW fitter
	    countFitValid++;//count events of fitValid 
		u_fit=fitVertex.GetDirection().Unit();
		double tfit = fitVertex.GetTime();
		double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());//drive correction
		TVector3 pos_cor; 
		//!!!! Do drive correction here
        pos_cor = pos_fit;//driveP0*pos_fit+driveP1*u_fit; //cout<<pos_cor.Unit()*sourcePos.Unit()<<endl;
		double Xcor = pos_cor.X();double Ycor = pos_cor.Y();double Zcor = pos_cor.Z();
        
        for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
        {
	      if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID)
	      {
		   //std::cout<<"get fitted "<<calpmts.GetFECDCount()<<std::endl;
		   countSuccess++;//count events of fitValid && FECD 
		   hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
		   hFitTime->Fill(tfit);
		   hPosMag->Fill(pos_fit.Mag());
		   hfitXY->Fill(posX,posY);
		   hfitXZ->Fill(posX,posZ); 
		   hXcor->Fill(Xcor);hYcor->Fill(Ycor);hZcor->Fill(Zcor);
 		   //  std::cout<<"fitted time  "<<tfit<<std::endl;
		   Double_t radius = pos_fit.Mag();
		   //Double_t fEnergy = fitVertex.GetEnergy();
		   //hfEnergy->Fill(fEnergy);
		   double temp1fit = pow(u_fit*pos_fit,2)-pos_fit.Mag2()+rPSUP*rPSUP;
		   double dPSUPfit=-u_fit*pos_fit+sqrt(temp1fit);
		   hMfit->Fill(dPSUPfit); //distance

		   for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
		   {
			vector<double> tPrompt1, tPrompt2; 
			double t1prompt=0, t2prompt=0; 
			TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
			hPMTpos->Fill(pmtpos.Mag());
			double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
			double tRes = hitTime-fitVertex.GetTime()-(pmtpos-pos_cor).Mag()/grVelocity;
			htResVsPMTposx->Fill(pmtpos.X(),tRes);htResVsPMTposy->Fill(pmtpos.Y(),tRes);htResVsPMTposz->Fill(pmtpos.Z(),tRes);
			htResVsPMTphi->Fill(pmtpos.Phi(),tRes);htResVsPMTtheta->Fill(cos(pmtpos.Theta()),tRes);
			hPMTphiTheta->Fill(pmtpos.Phi(),cos(pmtpos.Theta()));
		   //pmt cuts
			for(unsigned int jpmt=0;jpmt<calpmts.GetCount();jpmt++)
			{
			  TVector3 pmtposj = pmtInfo.GetPosition(calpmts.GetPMT(jpmt).GetID());
			  double hitTimej =(calpmts.GetPMT(jpmt)).GetTime();
			  double tResj= hitTimej-fitVertex.GetTime()-(pmtposj-pos_cor).Mag()/grVelocity;
			  		 
			  if(pmtposj.Unit()*pmtpos.Unit()>1-PMTsep*PMTsep/(2*8400*8400) && pmtposj.Unit()*pmtpos.Unit()!=1)
			  {
				htRes_inPMTsep->Fill(tResj);
				hPMTphiTheta_inCone->Fill(pmtposj.Phi(),cos(pmtposj.Theta()));
				if(tResj<1 && tResj>0)
				{
				 htPrompt_inPMTsep->Fill(hitTimej-tResj-fitVertex.GetTime());
				 tPrompt1.push_back(hitTimej-tResj-fitVertex.GetTime());
				} 
			  }
			  if(-pmtposj.Unit()*pmtpos.Unit()>1-PMTsep*PMTsep/(2*8400*8400) && pmtposj.Unit()*pmtpos.Unit()!=1)
			  {
				htRes_oppPMTsep->Fill(hitTimej-tResj-fitVertex.GetTime());
				hPMTphiTheta_oppCone->Fill(pmtposj.Phi(),cos(pmtposj.Theta()));
				if(tResj<1 && tResj>0)
				{
				 htPrompt_oppPMTsep->Fill(hitTimej-tResj-fitVertex.GetTime());
				 tPrompt2.push_back(hitTimej-tResj-fitVertex.GetTime());
				 }
			  }
			}//loop PMT j

			 if(tPrompt1.size()!=0)
			 {
			  double sum = 0;
			  for(std::vector<double>::iterator it = tPrompt1.begin(); it != tPrompt1.end(); ++it)
			  {double ss = *it; 
 			   sum = sum+ss;
			  }
			  t1prompt = sum/tPrompt1.size();
		         } 	
			 if(tPrompt2.size()!=0)
			 {
			  double sum2 = 0;
			  for(std::vector<double>::iterator it1 = tPrompt2.begin(); it1 != tPrompt2.end(); ++it1)
			  {double ss2 = *it1; 
			   sum2 = sum2+ss2;
			  }
			  t2prompt = sum2/tPrompt2.size();
		     	 } 
			if(t1prompt!=0 && t2prompt!=0)
			{
				double tSum = t1prompt + t2prompt;
				hSumTOFvsPMTposx->Fill(pmtpos.X(),tSum);  
				hSumTOFvsPMTposy->Fill(pmtpos.Y(),tSum);    
				hSumTOFvsPMTposz->Fill(pmtpos.Z(),tSum);    
			}		
	
		    if(trig&(trigWord))
		    { 
		    htResVsPMTposx_trig->Fill(pmtpos.X(),tRes);htResVsPMTposy_trig->Fill(pmtpos.Y(),tRes);htResVsPMTposz_trig->Fill(pmtpos.Z(),tRes);
          	    htResVsPMTphi_trig->Fill(pmtpos.Phi(),tRes);htResVsPMTtheta_trig->Fill(cos(pmtpos.Theta()),tRes);
			}
		   }//loop calibrated PMTs
		   if(trig&(trigWord)) {hFitTime_trig->Fill(tfit);}
			   
		   if(trig&(trigWord))
		   { 
		    if(nhits>nhitCut)
		    {
		     countTimeWindow++;//count events of fitValid && FECD && tFit Cut
		     hXcor_trig->Fill(Xcor);hYcor_trig->Fill(Ycor);hZcor_trig->Fill(Zcor);
 		     hXZcor_trig->Fill(Xcor,Zcor);hRZcor_trig->Fill(sqrt(Xcor*Xcor+Ycor*Ycor),Zcor);
		     hDeltaXcor_trig->Fill(Xcor-sourcePos.X());hDeltaYcor_trig->Fill(Ycor-sourcePos.Y());hDeltaZcor_trig->Fill(Zcor-sourcePos.Z());
		     hDeltaRcor_trig->Fill(pos_cor.Mag()-sourcePos.Mag());
		     hTheta->Fill((pos_cor-sourcePos).Unit()*u_fit);
 		     double source2pos = (pos_cor-sourcePos).Mag();
		     if(source2pos>1500) hThetaCut1->Fill((pos_cor-sourcePos).Unit()*u_fit);
		     if(source2pos>1000) hThetaCut2->Fill((pos_cor-sourcePos).Unit()*u_fit);
		     if(source2pos>750) hThetaCut3->Fill((pos_cor-sourcePos).Unit()*u_fit);
		     if(source2pos>1500) hThetaCut4->Fill((pos_cor-sourcePos).Unit()*u_fit);
		     if(source2pos>1200) hThetaCut5->Fill((pos_cor-sourcePos).Unit()*u_fit);    
	    	  }//nhitCuts  
	        }//trigger cuts; trigger2 =  trigger2>>1;
     	countFECDtotal++;     
		   }//FECD==9188
         } //FECD loop
       }//if fitValid <9m
     }//if fitValid
   }//if passed dataCleaning 
  }//if triggered events
 }//event loop
 
//  TTree *get = (TTree*)f1->Get("T");
//  TH1F *hFECD = new TH1F("hFECD","temp",100,9100,9400);
//  get->Draw("ds.evs.calPMTs.fecd.GetID()>>hFECD");
 std::cout<<"before clean total "<<countNoCleanTotal<<" triggered "<<countNoCleanTrig<<" FECD "<<countNoCleanFECD<<std::endl;
 std::cout<<"cleaned data "<<countClean<<" clean+triggered "<<countTrig<<" clean+FECD "<<countFECDtotal<<" clean+triggered+FECD "<<countTrigFECD<<std::endl;
 std::cout<<"fitValid/total events: "<<countFitValid<<"/"<<countClean<<" fitRate: "<<countFitValid/countClean*100<<std::endl;
 //  std::cout<<"fitValid && FECD/FECD total(9188+9027): "<<countSuccess<<"/"<<hFECD->GetEntries()<<" ratio: "<<countSuccess/hFECD->GetEntries()*100<<std::endl; //fitValid && FECD/ FECD, FECD means FECD == 9188, FECD total means 9188+9207
 std::cout<<"fitValid && FECD/FECD: "<<countSuccess<<"/"<<countFECDtotal<<" ratio: "<<countSuccess/countFECDtotal*100<<std::endl;
 std::cout<<"fitValid && FECD && trigger/FECD: "<<countTimeWindow<<"/"<<countFECDtotal<<" ratio: "<<countTimeWindow/countFECDtotal*100<<std::endl;
 std::cout<<"fitValid && FECD/fitValid: "<<countSuccess<<"/"<<countFitValid<<" ratio: "<<countSuccess/countFitValid*100<<std::endl;
 std::cout<<"fitValid && FECD && trigger/fitValid: "<<countTimeWindow<<"/"<<countFitValid<<" ratio: "<<countTimeWindow/countFitValid*100<<std::endl;
 std::cout<<"fitValid && FECD && trigger/fitValid && FECD: "<<countTimeWindow<<"/"<<countSuccess<<" ratio: "<<countTimeWindow/countSuccess*100<<std::endl;
  
  TString newfilename = "ResolTest_"+filename; 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  
  hPMTpos->Write(); 
  hfitX->Write();hfitY->Write();hfitZ->Write();hfitXY->Write();hfitXZ->Write();hfitRvsZ->Write();hPosMag->Write(); 
  hPosXvsTime->Write();
  hFitTime->Write();hHitTime->Write();htRes->Write();htRes_trig->Write();
  hTheta->Write();hThetaCut1->Write();hThetaCut2->Write();hThetaCut3->Write();hThetaCut4->Write(); hThetaCut5->Write();
  hPMTphiTheta->Write();hPMTphiTheta_inCone->Write();hPMTphiTheta_oppCone->Write();

  hXcor_trig->Write();hYcor_trig->Write();hZcor_trig->Write();
  hXZcor_trig->Write();hRZcor_trig->Write();
  hDeltaXcor_trig->Write();hDeltaYcor_trig->Write();hDeltaZcor_trig->Write();hDeltaRcor_trig->Write();

  hFitTime_trig->Write();
  htRes_inPMTsep->Write();htRes_oppPMTsep->Write();htPrompt_inPMTsep->Write();htPrompt_oppPMTsep->Write();
  htResVsPMTposx->Write();htResVsPMTposy->Write();htResVsPMTposz->Write();htResVsPMTphi->Write();htResVsPMTtheta->Write();
  htResVsPMTposx_trig->Write();htResVsPMTposy_trig->Write();htResVsPMTposz_trig->Write();htResVsPMTphi_trig->Write();htResVsPMTtheta_trig->Write();
  hSumTOFvsPMTposx->Write();hSumTOFvsPMTposy->Write();hSumTOFvsPMTposz->Write();    
  fp->Close();
  
}
