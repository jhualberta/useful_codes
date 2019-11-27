//For Rat-v6.2.8 SNOP data
//2017-6-25
//#include <RAT/DataCleaningUtility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DataCleaningUtility.hh>
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
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = 
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;
void analyPartialMPtRes()
{  
   const char* filename = "FitMP_waterLevelCorrection_PartialScintNeck_FillC14_Scint_r1_s0_p1.root";
   double driveP0 = 0.995765, driveP1 = -63.826;
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity

   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hTheta = new TH1F("hcosThetaSNO_", "SNO+ MPW (fitPos-sourcePos)*fitDirec", 1000,-1,1);
   TH1F* hThetaCut1 = new TH1F("hcosThetaCut1SNO_", "SNO+ H_{2}O fitted (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.5 m", 1000,-1,1);
   TH1F* hThetaCut2 = new TH1F("hcosThetaCut2SNO_", "SNO+ H_{2}O fitted (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.0 m", 1000,-1,1);
   TH1F* hThetaCut3 = new TH1F("hcosThetaCut3SNO_", "SNO+ H_{2}O fitted (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.75m", 1000,-1,1);
   TH1F* hThetaCut4 = new TH1F("hcosThetaCut4SNO_", "SNO+ H_{2}O fitted (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.5 m", 1000,-1,1);
   TH1F* hThetaCut5 = new TH1F("hcosThetaCut5SNO_", "SNO+ H_{2}O fitted (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.2 m", 1000,-1,1);
   TH1F* hFOM = new TH1F("hFOM", "FOM MPW", 1000,0,1000); 
   TH2F* hFOMvsNhit = new TH2F("hFOMvsNhit","FOM vs Nhit",150,0,150,1500,0,1500);
   TH2F* hFOMvsFitT = new TH2F("hFOMvsFitT","FOM vs Fitted Time",300,0,300,1500,0,1500); 
   TH2F* hFOMvsFitX = new TH2F("hFOMvsFitX","FOM vs Fitted X",2000,-6000,6000,1500,0,1500);
   TH2F* hFOMvsFitY = new TH2F("hFOMvsFitY","FOM vs Fitted Y",2000,-6000,6000,1500,0,1500);
   TH2F* hFOMvsFitZ = new TH2F("hFOMvsFitZ","FOM vs Fitted Z",2000,-6000,6000,1500,0,1500);
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
   // use lightPathCalculator
   TH1F* htRes_lightPath = new TH1F("htRes_lightPath", "Time Residue, trig", 400,-100,300);
   // use lightPathCalculator by default, check boundary
   TH1F* htRes_interface = new TH1F("htRes_interface", "Time Residue near boundary", 400,-100,300);
   TH1F* htRes_av = new TH1F("htRes_av", "Time Residue near AV", 400,-100,300);

   TH1F* htResMC_lightPath = new TH1F("htResMC_lightPath", "MC Time Residue, trig", 400,-100,300);
   TH1F* htResMC_interface = new TH1F("htResMC_interface", "MC Time Residue near boundary", 400,-100,300);
   TH1F* htResMC_av = new TH1F("htResMC_av", "MC Time Residue near AV", 400,-100,300);

   // use MPW straight light path
   TH1F* htRes = new TH1F("htRes_1p38486", "Time Residue, default = c/1.38486, trig", 400,-100,300);
   TH1F* htRes_rmOdd = new TH1F("htRes_rmOdd", "Time Residue, default = c/1.38486, remove odd PMTs", 1600,-100,300);
   TH1F* htRes_1p40 = new TH1F("htRes_1p40", "Time Residue, c/1.40, trig", 400,-100,300);
   TH1F* htRes_1p378 = new TH1F("htRes_1p378", "Time Residue, c/1.378, trig", 400,-100,300);
   // Check PMTs
   TH1F* hPMTid = new TH1F("hPMTid", "all PMT id", 10000, 0, 10000);
   TH1F* hPMTid_odd = new TH1F("hPMTid_odd", "PMT id in [-10,-5]ns ", 10000, 0, 10000);
   TH2F* hPMTthetaPhi = new TH2F("hPMTthetaPhi", "PMT distributions", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1); 
   TH2F* hPMTthetaPhi_odd = new TH2F("hPMTthetaPhi_odd", "PMT distributions in [-10,-5] ns", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1);

   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O fitted X, FECD cut", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O fitted Y, FECD cut", 2000,-9000,9000);  
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O fitted Z, FECD cut", 2000,-9000,9000);
   TH2F* hfitXY = new TH2F("hfitXY", "SNO+ H_{2}O fitted X vs Y, FECD cut", 2000,-9000,9000,2000,-9000,9000); 
   TH2F* hfitXZ = new TH2F("hfitXZ", "SNO+ H_{2}O fitted X vs Z, FECD cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ = new TH2F("hfitRZ", "SNO+ H_{2}O fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag = new TH1F("hPosMag","SNO+ H_{2}O fitted Mag()",1000,0,10000);
   TH2F* hfitRZ_odd = new TH2F("hfitRZ_odd", "SNO+ H_{2}O fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);


   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+ H_{2}O fitted X, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+ H_{2}O fitted Y, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+ H_{2}O fitted Z, FECD+trig cut", 2000,-9000,9000);
   TH2F* hfitXY_trig = new TH2F("hfitXY_trig", "SNO+ H_{2}O fitted X vs Y, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_trig = new TH2F("hfitXZ_trig", "SNO+ H_{2}O fitted X vs Z, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "SNO+ H_{2}O fitted R vs Z, FECD+trig cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_trig = new TH1F("hPosMag_trig","SNO+ H_{2}O fitted Mag(), FECD+trig cut",1000,0,10000);

   TH1F* hfitX_trigITR_beta14 = new TH1F("hfitX_trigITR_beta14", "SNO+ H_{2}O fitted X, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH1F* hfitY_trigITR_beta14  = new TH1F("hfitY_trigITR_beta14", "SNO+ H_{2}O fitted Y, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH1F* hfitZ_trigITR_beta14  = new TH1F("hfitZ_trigITR_beta14", "SNO+ H_{2}O fitted Z, FECD,trig,ITR,beta14", 2000,-9000,9000);
   TH2F* hfitXY_trigITR_beta14  = new TH2F("hfitXY_trigITR_beta14", "SNO+ H_{2}O fitted X vs Y, FECD,trig,ITR,beta14", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_trigITR_beta14  = new TH2F("hfitXZ_trigITR_beta14", "SNO+ H_{2}O fitted X vs Z, FECD,trig,ITR,beta14", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_trigITR_beta14  = new TH2F("hfitRZ_trigITR_beta14", "SNO+ H_{2}O fitted R vs Z, FECD,trig,ITR,beta14", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_trigITR_beta14 = new TH1F("hPosMag_trigITR_beta14","SNO+ H_{2}O fitted Mag(), FECD,trig,ITR,beta14",1000,0,10000);

   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 500,0,500);
   TH1F* hFitTime_odd = new TH1F("hFitTime_odd", "odd fitted time", 500,0,500);

   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "fitted time,FECD cuts, triggered", 800,0,800); 
   TH2F* hFlatVsZ = new TH2F("hFlatVsZ","fitted (x^2+y^2)/6000^2 vs Z",5000,0,5,2000,-9000,9000);
   TH1F* hfEnergy = new TH1F("hfEnergy","fitted Energy",100,0,100); 
   TH1F* hRSPfEnergy = new TH1F("hRSPfEnergy","RSPs fEnergy, 0<E<100",1000,0,100);
   TH1F* hRSPfNhits = new TH1F("hRSPfNhits","RSPs fNhits",100,0,100); 
   TH1F* hNhits = new TH1F("hNhits","Nhits",100,0,100); 
   TH1F* hNhitsFECD = new TH1F("hNhitsFECD_SNO_","Nhits FECD",100,0,100);
   TH2F* hNhitsVsMag = new TH2F("hNhitsVsMag","pos.Mag() vs Nhits",1000,0,10000,100,0,100);
   TH2F* hNhitVsFlat = new TH2F("hNhitsVsFlat","Nhits vs  (x^2+y^2)/6000^2",100,0,100,5000,0,5);
   TVector3 u_fit, pos_fit, pos_cor,  pos_fit_ItrBeta14, u_fit_pos_ItrBeta14;
   Double_t theta_e;

   vector<TH2F*> hCrateCardChannel;
   for(int i = 0;i<19;i++)
   {
     TH2F *htemp = new TH2F("htemp","card, channel", 16, 0, 15, 32,0,31);
     hCrateCardChannel.push_back((TH2F*)htemp->Clone(Form("hcrate%u",i)));
     delete htemp;
   }

   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   Double_t c_light = 299.792458;
   double n1 = 1.38486;
   double n2 = 1.40;
   double n3 = 1.378;

   Double_t rPSUP = 8390;
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
   unsigned int fecdID = 9207;
   TString fitName = "partialFitter";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
   {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     Int_t nevC = rDS.GetEVCount();
     const TVector3 pos_mc = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
     for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     double nhits = rev.GetNhits();
     hNhits->Fill(nhits);
     int trig = rev.GetTrigType();//get triggerType
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     countNoCleanTotal++;
     if(trig&(trigWord)){countNoCleanTrig++;}
     //RAT::DS::DataQCFlags dcflags = rev.GetDataCleaningFlags();
     //if( RAT::EventIsClean( rev, dcAnalysisWord ) )
     {
      countClean++;//count cleaned total events
      if(trig&(trigWord)){countTrig++;}
  
      for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)
      {
	   if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID){
		countFECDtotal++;
	     if(trig&(trigWord)){countTrigFECD++;}	
	   }
      }//count total FECD events

      if(rev.FitResultExists("partialFitter")){//find partialFitter exists
       try{
         RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
         if(rev.GetFitResult("partialFitter").GetValid())//Global Validity of the fitter !!! NOTE: this is different
         {
         //rev.GetFitResult("partialFitter").GetFOMNames()[0];
         //std::string FOMname ="partialFitter";
         //double fomValue =  rev.GetFitResult("partialFitter").GetFOM(FOMname);
         //hFOM->Fill(fomValue); hFOMvsNhit->Fill(rev.GetNhits(),fomValue);//cout<<rev.GetNhits()<<" "<<fomValue<<endl;

          pos_fit=fitVertex.GetPosition(); 
          //cout<<rev.GetClassifierNames().size()<<endl;//3 classifiers
          //cout<<rev.GetClassifierNames()[0]<<" "<<rev.GetClassifierNames()[1]<<" "<<rev.GetClassifierNames()[2]<<endl;
          //ITR:partialFitter QPDT:partialFitter isotropy:partialFitter

          //hFOM->Fill(fomValue); hFOMvsNhit->Fill(rev.GetNhits(),fomValue);//cout<<rev.GetNhits()<<" "<<fomValue<<endl;

          double tfit = fitVertex.GetTime();
          //if(fitVertex.ValidPosition() && fVertex1.ValidDirection())//NOTE: check fitPosition Valid !!
          { 
            double posX=(pos_fit.X()); double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
          
          for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
          {
	      if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID)
	      {
		   //std::cout<<"get fitted "<<calpmts.GetFECDCount()<<std::endl;
		   countSuccess++;//count events of fitValid && FECD 
		   //double knhit=rev.GetNhits(),kFOM = fomValue;
		   double kk = 7.97993; double mb = 19.0462;//5 MeV
		   hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
		   hFitTime->Fill(tfit);
		   hPosMag->Fill(pos_fit.Mag());
		   hfitXY->Fill(posX,posY);
		   hfitXZ->Fill(posX,posZ); 
		   hfitRZ->Fill(sqrt(posX*posX+posY*posY),posZ);
		   hFlatVsZ->Fill((posX*posX+posY*posY)/(6000*6000),posZ);
		   hPhiTheta->Fill(pos_fit.Phi(),pos_fit.CosTheta());
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
                        lightPath.CalcByPosition( pos_fit, pmtpos );
                        double distInInnerAV = lightPath.GetDistInInnerAV();
                        double distInAV = lightPath.GetDistInAV();
                        double distInWater = lightPath.GetDistInWater();
                        const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
                        htRes_lightPath->Fill( hitTime - transitTime - tfit);
                        //!!! z correction for cut
			TVector3 pos_fit_cor(pos_fit.X(), pos_fit.Y(), pos_fit.Z()-108);
			if(pos_fit_cor.Z()>4425 && pos_fit_cor.Z()<4450 ) htRes_interface->Fill( hitTime - transitTime - tfit);
                        if(pos_fit_cor.Mag()>5999 && pos_fit_cor.Mag()<6005 ) htRes_av->Fill( hitTime - transitTime - tfit);

			// Monte Carlo
                        lightPath.CalcByPosition( pos_mc, pmtpos );
                        double distInInnerAV = lightPath.GetDistInInnerAV();
                        double distInAV = lightPath.GetDistInAV();
                        double distInWater = lightPath.GetDistInWater();
                        const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
                        // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time
                        // hit times are relative to the trigger time, which will depend on event time and detector position so correct for that to line up events
                        // The 390ns corrects for the electronics delays and places the pulse in the middle of the window
                        double tResMC = hitTime - transitTime - 390 + rDS.GetMCEV(iev).GetGTTime();
			htResMC_lightPath->Fill(tResMC); 
                        //!!!  z correction for cut
                        TVector3 pos_mc_cor(pos_mc.X(), pos_mc.Y(), pos_mc.Z()-108);
			if(pos_mc_cor.Z()>4425 && pos_mc_cor.Z()<4450 ) htResMC_interface->Fill(tResMC);
                        if(pos_mc_cor.Mag()>5999 && pos_mc_cor.Mag()<6005 ) htResMC_av->Fill(tResMC);

			double tRes1 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n1);
                        double tRes2 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n2);
                        double tRes3 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n3);
			hPosXvsTime->Fill(hitTime,posX);
                        //if(trig&(2)
                        hPMTid->Fill(calpmts.GetPMT(ipmt).GetID());
			hPMTthetaPhi->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                        if(tRes1<-5 && tRes1>-10)
			{
        		  hPMTthetaPhi_odd->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                          hPMTid_odd->Fill(calpmts.GetPMT(ipmt).GetID());
                          hfitRZ_odd->Fill(sqrt(posX*posX+posY*posY),posZ);
                          hFitTime_odd->Fill(tfit);
                          int crate = calpmts.GetPMT(ipmt).GetCrate();
			  hCrateCardChannel[crate]->Fill(calpmts.GetPMT(ipmt).GetCard(),calpmts.GetPMT(ipmt).GetChannel());
			}
                       if(calpmts.GetPMT(ipmt).GetCrate()!=8 && calpmts.GetPMT(ipmt).GetCard()!=9)
		       { htRes_rmOdd->Fill(tRes1);}

			htRes->Fill(tRes1);
                        htRes_1p40->Fill(tRes2);
			htRes_1p378->Fill(tRes3);
		   }
			   
		     hTrig->Fill(0);
		     countTimeWindow++;//count events of fitValid && FECD && tFit Cut
		     hfitX_trig->Fill(posX);hfitY_trig->Fill(posY);hfitZ_trig->Fill(posZ);hfitXZ_trig->Fill(posX,posZ);
		     hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);hPosMag_trig->Fill(pos_fit.Mag());
		
		     //hTheta->Fill((pos_fit-sourcePos).Unit()*u_fit);
 		     double source2pos = (pos_fit-sourcePos).Mag();
		     //if(source2pos>1500) hThetaCut1->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1000) hThetaCut2->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>750) hThetaCut3->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1500) hThetaCut4->Fill((pos_fit-sourcePos).Unit()*u_fit);
		     //if(source2pos>1200) hThetaCut5->Fill((pos_fit-sourcePos).Unit()*u_fit);
 		     hFitTime_trig->Fill(tfit);
	    }//FECD==9188
          } //FECD loop
	}//fit position valid //!!!NOTE: this is necessary 
       }//if water fitter global valid //!!!NOTE: this is necessary
       }//try catch
       catch(exception& e)
       {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
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
  TString newfilename = "ResolMP_"+TString(filename); 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  

  for (size_t j=0;j<19;j++)
  {
    hCrateCardChannel[j]->Write();
  }

  hTrig->Write(); 
  hfitX->Write();hfitY->Write();hfitZ->Write();hfitXY->Write();hfitXZ->Write();hfitRZ->Write();hPosMag->Write();
  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hfitXZ_trig->Write();hfitRZ_trig->Write();hPosMag_trig->Write();

  hfitX_trigITR_beta14->Write();hfitY_trigITR_beta14->Write();hfitZ_trigITR_beta14->Write();hfitXY_trigITR_beta14->Write();
  hfitXZ_trigITR_beta14->Write();hfitRZ_trigITR_beta14->Write();hPosMag_trigITR_beta14->Write();
  hfitRZ_odd->Write();hFitTime_odd->Write();

  hNhitsFECD->Write();hPosXvsTime->Write();
  hFOM->Write();hFitTime->Write();hHitTime->Write();htRes->Write(); htRes_rmOdd->Write();
  htRes_lightPath->Write();htRes_1p40->Write();htRes_1p378->Write();
  htRes_interface->Write();htRes_av->Write();

  htResMC_lightPath->Write();htResMC_interface->Write();htResMC_av->Write();

  hPMTid->Write();hPMTid_odd->Write();hPMTthetaPhi->Write();hPMTthetaPhi_odd->Write();

  hFOMvsFitT->Write();hFOMvsFitX->Write();hFOMvsFitY->Write();hFOMvsFitZ->Write();hFOMvsNhit->Write();//hFitXYvsFOM->Write();
  hTheta->Write();hThetaCut1->Write();hThetaCut2->Write();hThetaCut3->Write();hThetaCut4->Write(); hThetaCut5->Write();
  fp->Close();
  
}
