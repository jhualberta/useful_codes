//For Rat-v6.2.8 SNOP data
//2017-6-25
//#include <RAT/DataCleaningUtility.hh>
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
#include "TF1.h" 
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = Rat
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;
double trigDelay = 250;
Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
//Double_t grVelocity = 183.503698029;//labppo 
Double_t rPSUP = 8900;

void analyRatWaterFOM()
{  
   const char* filename = "FitRat_MC_rat6168_realisticTrigger_run106904_Water_gamma_2p2MeV_isoFill.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hFOM = new TH1F("hFOM", "FOM MPW", 1000,0,1000);
   TH1F* hFOM_nhitScaled = new TH1F("hFOM_nhitScaled", "FOM scaled to nhit", 1000,0,1000);
   TH1F* hFOM_nhitSelected = new TH1F("hFOM_nhitSelected", "FOM scaled to nhit", 1000,0,1000);
   TH2F* hFOMvsNhit = new TH2F("hFOMvsNhit","FOM vs Nhit",150,0,150,1500,0,1500);

   TH1F* hTheta = new TH1F("hcosThetaSNO_Rat", "SNO+ MPW (fitPos-sourcePos)*fitDirec", 1000,-1,1);
   TH1F* hThetaCut1 = new TH1F("hcosThetaCut1_Rat", "SNO+ H_{2}O Rat fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.5 m", 1000,-1,1);
   TH1F* hThetaCut2 = new TH1F("hcosThetaCut2_Rat", "SNO+ H_{2}O Rat fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.0 m", 1000,-1,1);
   TH1F* hThetaCut3 = new TH1F("hcosThetaCut3_Rat", "SNO+ H_{2}O Rat fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.75m", 1000,-1,1);
   TH1F* hThetaCut4 = new TH1F("hcosThetaCut4_Rat", "SNO+ H_{2}O Rat fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.5 m", 1000,-1,1);
   TH1F* hThetaCut5 = new TH1F("hcosThetaCut5_Rat", "SNO+ H_{2}O Rat fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.2 m", 1000,-1,1);
   TH2F* hFOM_nhitSelectedVsFitT = new TH2F("hFOM_nhitSelectedVsFitT","logL_selectedNhit vs Fitted Time",300,0,300,150,0,150); 
   TH2F* hFOM_nhitSelectedVsFitX = new TH2F("hFOM_nhitSelectedVsFitX","logL_selectedNhit vs Fitted X",2000,-6000,6000,150,0,150);
   TH2F* hFOM_nhitSelectedVsFitY = new TH2F("hFOM_nhitSelectedVsFitY","logL_selectedNhit vs Fitted Y",2000,-6000,6000,150,0,150);
   TH2F* hFOM_nhitSelectedVsFitZ = new TH2F("hFOM_nhitSelectedVsFitZ","logL_selectedNhit vs Fitted Z",2000,-6000,6000,150,0,150);
   TH2F* hPhiTheta = new TH2F("hPhiTheta","theta vs phi",1000,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);

//check pos.Mag() vs fitTime

  //for N16 source 
   TH1F* htRes_oddT_trig = new TH1F("htRes_oddT_trig","time residue trigger, tFit<100ns,trig", 400,-100,300);
   //corrected values with non FECD cuts
   TH1F* hFitTime_noFECD = new TH1F("hFitTime_noFECD","fitted time, noFECDs",800,0,800);
   TH1F* htRes_oddT_noFECD = new TH1F("htRes_oddT_noFECD","time residue trigger, tFit<100ns,noFECD", 400,-100,300);

   TH1F* hFitTime_noFECD_trig = new TH1F("hFitTime_noFECD_trig","fitted time, noFECD but trig",800,0,800);
   TH1F* htRes_oddT_noFECD_trig = new TH1F("htRes_oddT_noFECD_trig","time residue trigger, tFit<100ns,noFECD but trig", 400,-100,300);

   TH2F* hPosXvsTime = new TH2F("hPosXVsTime","hit time vs posX",100,0,500,2000,-6000,6000);
   TH1F *hMfit = new TH1F("hMfit","Mfit: pos_fit to PSUP sphere",1000,1000,15000);  
   TH1F* htRes = new TH1F("htRes", "Time Residue", 400,-100,300);
   TH1F* htRes_trig = new TH1F("htRes_trig", "Time Residue, trig", 400,-100,300);
   TH1F* htResTrue_trig = new TH1F("htResTrue_trig", "MCtrue Time Residue, trig", 1600,-100,300);
      
   TH1F* hfitX_noFECD = new TH1F("hfitX_noFECD", "SNO+ H_{2}O Rat fitter X, no FECD", 2000,-9000,9000);
   TH1F* hfitY_noFECD = new TH1F("hfitY_noFECD", "SNO+ H_{2}O Rat fitter Y, no FECD", 2000,-9000,9000);
   TH1F* hfitZ_noFECD = new TH1F("hfitZ_noFECD", "SNO+ H_{2}O Rat fitter Z, no FECD", 2000,-9000,9000);
   TH2F* hfitXY_noFECD = new TH2F("hfitXY_noFECD", "SNO+ H_{2}O Rat fitter X vs Y, no FECD", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_noFECD = new TH2F("hfitXZ_noFECD", "SNO+ H_{2}O Rat fitter X vs Z, no FECD", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_noFECD = new TH2F("hfitRZ_noFECD", "SNO+ H_{2}O Rat fitter R vs Z, no FECD", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_noFECD = new TH1F("hPosMag_noFECD","SNO+ H_{2}O Rat fitter Mag(), no FECD",1000,0,10000);

   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O Rat fitter X, FECD cut", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O Rat fitter Y, FECD cut", 2000,-9000,9000);  
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O Rat fitter Z, FECD cut", 2000,-9000,9000);
   TH2F* hfitXY = new TH2F("hfitXY", "SNO+ H_{2}O Rat fitter X vs Y, FECD cut", 2000,-9000,9000,2000,-9000,9000); 
   TH2F* hfitXZ = new TH2F("hfitXZ", "SNO+ H_{2}O Rat fitter X vs Z, FECD cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ = new TH2F("hfitRZ", "SNO+ H_{2}O Rat fitter R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag = new TH1F("hPosMag","SNO+ H_{2}O Rat fitter Mag()",1000,0,10000);

   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+ H_{2}O Rat fitter X, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+ H_{2}O Rat fitter Y, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+ H_{2}O Rat fitter Z, FECD+trig cut", 2000,-9000,9000);
   TH2F* hfitXY_trig = new TH2F("hfitXY_trig", "SNO+ H_{2}O Rat fitter X vs Y, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_trig = new TH2F("hfitXZ_trig", "SNO+ H_{2}O Rat fitter X vs Z, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "SNO+ H_{2}O Rat fitter R vs Z, FECD+trig cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_trig = new TH1F("hPosMag_trig","SNO+ H_{2}O Rat fitter Mag(), FECD+trig cut",1000,0,10000);
   
   TH1F* hDeltaX = new TH1F("hDeltaX", "Rat fitted X - mc X, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaY = new TH1F("hDeltaY", "Rat fitted Y - mc Y, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaZ = new TH1F("hDeltaZ", "Rat fitted Z - mc Z, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaPosMag = new TH1F("hDeltaPosMag","MP fitted - MC Mag()",1000,0,10000);

   TH1F* hDeltaX_trig = new TH1F("hDeltaX_trig", "Rat fitted X - mc X, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaY_trig = new TH1F("hDeltaY_trig", "Rat fitted Y - mc Y, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaZ_trig = new TH1F("hDeltaZ_trig", "Rat fitted Z - mc Z, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaPosMag_trig = new TH1F("hDeltaPosMag_trig","Rat fitted - MC Mag(), trigCut",1000,0,10000);

   TH1F* hfitX_trigITR_beta14 = new TH1F("hfitX_trigITR_beta14", "SNO+ H_{2}O Rat fitter X, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH1F* hfitY_trigITR_beta14  = new TH1F("hfitY_trigITR_beta14", "SNO+ H_{2}O Rat fitter Y, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH1F* hfitZ_trigITR_beta14  = new TH1F("hfitZ_trigITR_beta14", "SNO+ H_{2}O Rat fitter Z, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH2F* hfitXY_trigITR_beta14  = new TH2F("hfitXY_trigITR_beta14", "SNO+ H_{2}O Rat fitter X vs Y, FECD,trig,ITR,beta14", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_trigITR_beta14  = new TH2F("hfitXZ_trigITR_beta14", "SNO+ H_{2}O Rat fitter X vs Z, FECD,trig,ITR,beta14", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_trigITR_beta14  = new TH2F("hfitRZ_trigITR_beta14", "SNO+ H_{2}O Rat fitter R vs Z, FECD,trig,ITR,beta14", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_trigITR_beta14 = new TH1F("hPosMag_trigITR_beta14","SNO+ H_{2}O Rat fitter Mag(), FECD,trig,ITR,beta14",1000,0,10000);

   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);
   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "fitted time,FECD cuts, triggered", 800,0,800); 
   TH2F* hFlatVsZ = new TH2F("hFlatVsZ","fitted (x^2+y^2)/6000^2 vs Z",5000,0,5,2000,-9000,9000);
   TH1F* hfEnergy = new TH1F("hfEnergy","fitted Energy",100,0,100); 
   TH1F* hRSPfEnergy = new TH1F("hRSPfEnergy","RSPs fEnergy, 0<E<100",1000,0,100);
   TH1F* hRSPfNhits = new TH1F("hRSPfNhits","RSPs fNhits",100,0,100); 
   TH1F* hNhits = new TH1F("hNhits","Nhits",100,0,100); 
   TH1F* hNhitsFECD = new TH1F("hNhitsFECD_SNO_Rat","Nhits FECD",100,0,100);
   TH2F* hNhitsVsMag = new TH2F("hNhitsVsMag","pos.Mag() vs Nhits",1000,0,10000,100,0,100);
   TH2F* hNhitVsFlat = new TH2F("hNhitsVsFlat","Nhits vs  (x^2+y^2)/6000^2",100,0,100,5000,0,5);
   TVector3 u_fit,pos_fit, pos_fit_ItrBeta14, u_fit_pos_ItrBeta14;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
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
   Double_t countITRbeta14 = 0;
   Double_t countNoCleanTotal=0; Double_t countNoCleanTrig=0;Double_t countNoCleanFECD=0;
   int trigWord = 0x11;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   int nhitCut = 0;//15;
   unsigned int fecdID = 9190;
   TString specFit = "waterFitter";
   //TString specFit = "fit10";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   bool trigCut = true; 
   if(trigCut) grVelocity = 183.503698029;
   TVector3 pos_mc;
   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
   {
     //std:cout << " event ID "<< iEntry <<std::endl ;
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    Int_t nevC = rDS.GetEVCount();
    const RAT::DS::MC& rmc= rDS.GetMC();  
    Int_t nMCEV = rDS.GetMCEVCount(); 
 
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
    pos_mc =rmcparticle.GetPosition();
    double tmc = rmcparticle.GetTime();
    for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     const RAT::DS::MCEV& mcEvent = rDS.GetMCEV(iev);  
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     double nhits = rev.GetNhits();
     hNhits->Fill(nhits);
     int trig = rev.GetTrigType();//get triggerType
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     
     countNoCleanTotal++;
     if(trig&(trigWord)){countNoCleanTrig++;}
     //if( RAT::EventIsClean( rev, dcAnalysisWord ) )
     {
      countClean++;//count cleaned total events
      if(trig&(trigWord)){countTrig++;}
  
      for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)
      {
	   if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID)
           {
		countFECDtotal++;
	     if(trig&(trigWord)){countTrigFECD++;}	
	   }
      }//count total FECD events

      if(rev.FitResultExists("waterFitter"))
      {//find waterFitter exists
       RAT::DS::FitVertex fitVertex = rev.GetFitResult("waterFitter").GetVertex(0);
       if(rev.GetFitResult("waterFitter").GetValid())//Global Validity of the fitter !!! NOTE: this is different
       {
	pos_fit=fitVertex.GetPosition(); 
	//u_fit=fitVertex.GetDirection().Unit();
 	double tfit = fitVertex.GetTime();
	if(fitVertex.ValidPosition())//NOTE: check fitPosition Valid !!
 	{ countFitValid++;//count events of fitValid 
	  double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
 	  hFitTime_noFECD->Fill(tfit);
          //if(trig&(trigWord))//trig cuts
	  {
	     hFitTime_noFECD_trig->Fill(tfit);
   	      for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
              {
                    TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                    double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                    double tRes = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/grVelocity;
		    htRes_oddT_noFECD_trig->Fill(tRes);
	      }	
	}//trigger = trigger>>1;
        
        //for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
        {
	      //if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID)
	      {
		   //std::cout<<"get fitted "<<calpmts.GetFECDCount()<<std::endl;
		   countSuccess++;//count events of fitValid && FECD 
		   //double knhit=rev.GetNhits(),kFOM = fomValue;
                  //vector<string> fomname = rev.GetFitResult("waterFitter").GetFOMNames();
                  //for(vector<string>::iterator i = fomname.begin(); i != fomname.end(); ++i) std::cout<<*i<<" ";
                 //cout<<endl;
                 /*DirectionLogL DirectionSelectedNHit EnergyGtest EnergyMedianDev EnergyMedianDevHit EnergyMedianProb EnergyMedianProbHit EnergyRSP EnergyUtest EnergynCerPhotons EnergypromptHits PositionLogL PositionSelectedNHit*/
                   double fomValue = rev.GetFitResult("waterFitter").GetFOM("PositionLogL");
                   double logLnhitSelected = rev.GetFitResult("waterFitter").GetFOM("PositionSelectedNHit");
                   hFOM->Fill(fomValue); hFOM_nhitScaled->Fill(fomValue/nhits);
                   hFOM_nhitSelected->Fill(logLnhitSelected);
                   hFOMvsNhit->Fill(fomValue,nhits);
                   hFOM_nhitSelectedVsFitT->Fill(tfit, logLnhitSelected);
                   hFOM_nhitSelectedVsFitX->Fill(posX-pos_mc.X(), logLnhitSelected);
                   hFOM_nhitSelectedVsFitY->Fill(posY-pos_mc.Y(), logLnhitSelected);
                   hFOM_nhitSelectedVsFitZ->Fill(posZ-pos_mc.Z(), logLnhitSelected);

		   hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
		   hFitTime->Fill(tfit);
		   hPosMag->Fill(pos_fit.Mag());
		   hfitXY->Fill(posX,posY);
		   hfitXZ->Fill(posX,posZ); 
		   hfitRZ->Fill(sqrt(posX*posX+posY*posY),posZ);
		   hFlatVsZ->Fill((posX*posX+posY*posY)/(6000*6000),posZ);
		   hPhiTheta->Fill(pos_fit.Phi(),pos_fit.CosTheta());
                   hDeltaX->Fill(posX-pos_mc.X()); hDeltaY->Fill(posY-pos_mc.Y()); hDeltaZ->Fill(posZ-pos_mc.Z());
                   hDeltaPosMag->Fill((pos_mc-pos_fit).Mag());

		   //hFitXYvsFOM->Fill(posX,posY,fomValue);
		   //hFOMvsFitT->Fill(tfit,fomValue);
		   //hFOMvsFitX->Fill(posX,fomValue);
		   //hFOMvsFitY->Fill(posY,fomValue);
		   //hFOMvsFitZ->Fill(posZ,fomValue);
		   ////
 		   //  std::cout<<"fitted time  "<<tfit<<std::endl;
		   Double_t radius = pos_fit.Mag();
		   //Double_t fEnergy = fitVertex.GetEnergy();
		   hNhitsFECD->Fill(rev.GetNhits());
		   hNhitsVsMag->Fill(pos_fit.Mag(),rev.GetNhits());
		   hNhitVsFlat->Fill(rev.GetNhits(),(posX*posX+posY*posY)/(6000*6000));
		   //hfEnergy->Fill(fEnergy);
		   //double temp1fit = pow(u_fit*pos_fit,2)-pos_fit.Mag2()+rPSUP*rPSUP;
		   //double dPSUPfit=-u_fit*pos_fit+sqrt(temp1fit);
		   //hMfit->Fill(dPSUPfit); //distance
		   for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++) 
		   {
		   	TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
			double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
			hHitTime->Fill(hitTime);
                        double mcHitTime = (mcEvent.GetMCHits().GetAllPMT( ipmt )).GetTime(); 
			double tRes = hitTime - fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/grVelocity;
                        double tResTrue = hitTime - (hitTime- mcHitTime)-(pmtpos-pos_mc).Mag()/grVelocity;
			hPosXvsTime->Fill(hitTime,posX);
			htRes->Fill(tRes);
			if(!trigCut || trig == 1)
			{
                          htRes_trig->Fill(tRes);
                          htResTrue_trig->Fill(tResTrue);
                        }
		   }
			   
		   if(!trigCut || trig == 1)
                   { 
		    if(nhits>nhitCut)
		    {
		      //if(rev.GetClassifierResult("ITR:waterFitter").GetValid() && rev.GetClassifierResult("isotropy:waterFitter").GetValid())
		      {
          		countITRbeta14++;
          		pos_fit_ItrBeta14 = fitVertex.GetPosition();
			double posXcuts = pos_fit_ItrBeta14.X();double posYcuts = pos_fit_ItrBeta14.Y();double posZcuts = pos_fit_ItrBeta14.Z();
			hfitX_trigITR_beta14->Fill(posXcuts);hfitY_trigITR_beta14->Fill(posYcuts); hfitZ_trigITR_beta14->Fill(posZcuts);  
			hfitXY_trigITR_beta14->Fill(posXcuts,posYcuts);hfitXZ_trigITR_beta14->Fill(posXcuts,posZcuts);
 			hfitRZ_trigITR_beta14->Fill(sqrt(posXcuts*posXcuts+posYcuts*posYcuts),posZcuts); 
 			hPosMag_trigITR_beta14->Fill(pos_fit_ItrBeta14.Mag());
    		      }

		     hTrig->Fill(0);
		     countTimeWindow++;//count events of fitValid && FECD && tFit Cut
		     hfitX_trig->Fill(posX);hfitY_trig->Fill(posY);hfitZ_trig->Fill(posZ);hfitXZ_trig->Fill(posX,posZ);
		     hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);hPosMag_trig->Fill(pos_fit.Mag());
                     hDeltaX_trig->Fill(posX-pos_mc.X()); hDeltaY_trig->Fill(posY-pos_mc.Y()); hDeltaZ_trig->Fill(posZ-pos_mc.Z());
                     hDeltaPosMag_trig->Fill((pos_mc-pos_fit).Mag());
		     //hTheta->Fill((pos_fit-sourcePos).Unit()*u_fit);
 		     double source2pos = (pos_fit-sourcePos).Mag();
		     //if(source2pos>1500) hThetaCut1->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1000) hThetaCut2->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>750) hThetaCut3->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1500) hThetaCut4->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1200) hThetaCut5->Fill((pos_fit-sourcePos).Unit()*u_fit);
 		     hFitTime_trig->Fill(tfit);
		    }//nhitCuts  
		   }//trigger cuts// trigger2 =  trigger2>>1;
	     }//FECD==9188
           } //FECD loop
	}//fit position valid //!!!NOTE: this is necessary 
       }//if water fitter global valid //!!!NOTE: this is necessary
     }//if water fitter exists
   }//if passed dataCleaning 
    }//if triggered events
  }
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
  std::cout<<"fitValid && FECD && trigger && ITR && beta14 "<<countITRbeta14<<std::endl;
  TString newfilename = "ResolRat_"+TString(filename); 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  
  TF1 *gx = new TF1("gx","gaus",hfitX_trig->GetMean()-2000,hfitX_trig->GetMean()+2000);
  TF1 *gy = new TF1("gy","gaus",hfitY_trig->GetMean()-2000,hfitY_trig->GetMean()+2000);
  TF1 *gz = new TF1("gz","gaus",hfitZ_trig->GetMean()-2000,hfitZ_trig->GetMean()+2000);
  gx->SetLineColor(kRed);gy->SetLineColor(kRed);gz->SetLineColor(kRed);
  if(trigCut) {
  hDeltaX_trig->Fit(gx);
  hDeltaY_trig->Fit(gy);
  hDeltaZ_trig->Fit(gz);
  }
  else {
  hDeltaX->Fit(gx,"R");
  hDeltaY->Fit(gy,"R");
  hDeltaZ->Fit(gz,"R");
  }

  hTrig->Write(); 
  hfitX->Write();hfitY->Write();hfitZ->Write();hfitXY->Write();hfitXZ->Write();hfitRZ->Write();hPosMag->Write();
  hfitX_noFECD->Write();hfitY_noFECD->Write();hfitZ_noFECD->Write();hfitXZ_noFECD->Write();hfitRZ_noFECD->Write();hPosMag_noFECD->Write();
  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hfitXZ_trig->Write();hfitRZ_trig->Write();hPosMag_trig->Write();

  hfitX_trigITR_beta14->Write();hfitY_trigITR_beta14->Write();hfitZ_trigITR_beta14->Write();hfitXY_trigITR_beta14->Write();
  hfitXZ_trigITR_beta14->Write();hfitRZ_trigITR_beta14->Write();hPosMag_trigITR_beta14->Write();

  hDeltaX->Write(); hDeltaY->Write(); hDeltaZ->Write();
  hDeltaPosMag->Write();
  hDeltaX_trig->Write();hDeltaY_trig->Write(); hDeltaZ_trig->Write();
  hDeltaPosMag_trig->Write();

  hNhitsFECD->Write();hPosXvsTime->Write();
  hFOM->Write();hFOMvsNhit->Write();hFOM_nhitScaled->Write();hFOM_nhitSelected->Write();
  hFitTime->Write();hHitTime->Write();htRes->Write();htRes_trig->Write();htResTrue_trig->Write();
  hFOM_nhitSelectedVsFitT->Write();hFOM_nhitSelectedVsFitX->Write();hFOM_nhitSelectedVsFitY->Write();hFOM_nhitSelectedVsFitZ->Write();

  hTheta->Write();hThetaCut1->Write();hThetaCut2->Write();hThetaCut3->Write();hThetaCut4->Write(); hThetaCut5->Write();
  hFitTime_noFECD_trig->Write();
  fp->Close();
  
}
