//2017-6-20
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
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = MultiWater
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;const double ITR_PMT_Hi = 5.0; const double ITR_PMT_Lo = -2.5;
void analyMPWnew()
{  
   TString filename = "FitMP_output.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader("FitMP_output.root");
   TVector3 sourcePos; 
   sourcePos.SetXYZ(-186,256,-4397);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hTheta = new TH1F("hcosThetaSNOP_MultiWater", "SNO+ MPW (fitPos-sourcePos)*fitDirec, with corrected Pos and tResCut", 1000,-1,1);
   TH1F* hThetaCut1 = new TH1F("hcosThetaCut1SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.5 m", 1000,-1,1);
   TH1F* hThetaCut2 = new TH1F("hcosThetaCut2SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.0 m", 1000,-1,1);
   TH1F* hThetaCut3 = new TH1F("hcosThetaCut3SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.75m", 1000,-1,1);
   TH1F* hThetaCut4 = new TH1F("hcosThetaCut4SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>0.5 m", 1000,-1,1);
   TH1F* hThetaCut5 = new TH1F("hcosThetaCut5SNOP_MultiWater", "SNO+ H_{2}O MultiWater fitter (fitPos-sourcePos)*fitDirec, |fitPos-sourcePos|>1.2 m", 1000,-1,1);
   TH1F* hFOM = new TH1F("hFOM", "FOM MPW", 1000,0,1000); 
   TH2F* hFOMvsNhit = new TH2F("hFOMvsNhit","FOM vs Nhit",150,0,150,1500,0,1500);
   TH2F* hFOMvsFitT = new TH2F("hFOMvsFitT","FOM vs Fitted Time",300,0,300,1500,0,1500); 
   TH2F* hFOMvsFitX = new TH2F("hFOMvsFitX","FOM vs Fitted X",2000,-6000,6000,1500,0,1500);
   TH2F* hFOMvsFitY = new TH2F("hFOMvsFitY","FOM vs Fitted Y",2000,-6000,6000,1500,0,1500);
   TH2F* hFOMvsFitZ = new TH2F("hFOMvsFitZ","FOM vs Fitted Z",2000,-6000,6000,1500,0,1500);

   TH2F* hPhiTheta = new TH2F("hPhiTheta","theta vs phi",628,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);

   TH1F* hXcor = new TH1F("hXcor","Xcor",2000,-9000,9000);
   TH1F* hYcor = new TH1F("hYcor","Ycor",2000,-9000,9000);
   TH1F* hZcor = new TH1F("hZcor","Zcor",2000,-9000,9000);

   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O MultiWater fitter X", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O MultiWater fitter Y", 2000,-9000,9000);   
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O MultiWater fitter Z", 2000,-9000,9000);

   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);

   TVector3 u_fit,pos_fit, pos_cor;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   Double_t rPSUP = 8900;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   Double_t energy, wavelength ;
 
   int trigWord = 2;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   int nhitCut = 15;
   unsigned int fecdID = 9188;
   TString specFit = "MultiPathProcessor";
   //TString specFit = "fit10";
   ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   size_t evtNum = dsReader.GetEntryCount(); 
   for( size_t iEntry = 0; iEntry <evtNum; iEntry++ )
   {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     Int_t nevC = rDS.GetEVCount();
     for(Int_t iev=0;iev<nevC; iev++) {
     /// Get the Event information    
     const RAT::DS::EV& ev = rDS.GetEV(iev);
     if(ev.FitResultExists("multipath")) {
       try {
         RAT::DS::FitVertex fVertex = ev.GetFitResult("multipath").GetVertex(0);
         if(fVertex.ValidPosition()) {
           double tfit = fVertex.GetTime();
           hFitTime->Fill(tfit);
           pos_fit = fVertex.GetPosition(); //Note: here the position result is without drive correction!
           RAT::DS::FitVertex fVertex1 =ev.GetFitResult("multipathdirection").GetVertex(0);
           double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
           hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
           if(fVertex1.ValidDirection()) {
             try{
               TVector3 u_fit = fVertex1.GetDirection();
               TVector3 pos_cor = fVertex1.GetPosition(); //Note: here the position result is after drive correction!
               double Xcor = pos_cor.X();double Ycor = pos_cor.Y();double Zcor = pos_cor.Z();
               hXcor->Fill(Xcor);hYcor->Fill(Ycor);hZcor->Fill(Zcor);
             } catch(exception& e) {std::cout<<e.what()<<"fit direction is not valid"<<std::endl;}
           } // if valid direction
         }// if valid position
       } 
       catch(exception& e) {std::cout<<e.what()<<"fit position is not valid"<<std::endl;}
     }// fit exists
   } // evt
  }// loop evt
  TString newfilename = "ResolCorrected_"+filename; 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  
  hfitX->Write();hfitY->Write();hfitZ->Write();
  hXcor->Write();hYcor->Write();hZcor->Write();

  fp->Close();
}
