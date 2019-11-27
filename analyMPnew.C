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
#include "TF1.h" 
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = MultiPath
fit9 = FTP
*/
using namespace std ;
//gStyle->SetPalette(53);
void analyMPnew()
{  
   const char *filename = "FitMP_fillTop_test.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   bool trigCut = true;
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);

   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O MultiPath fitter X", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O MultiPath fitter Y", 2000,-9000,9000);   
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O MultiPath fitter Z", 2000,-9000,9000);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);
   TH1F* hDeltaX = new TH1F("hDeltaX", "MP fitted X - mc X, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaY = new TH1F("hDeltaY", "MP fitted Y - mc Y, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaZ = new TH1F("hDeltaZ", "MP fitted Z - mc Z, FECD cut", 2000,-9000,9000);
   TH2F* hfitRZ = new TH2F("hfitRZ", "MP fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+ H_{2}O MultiPath fitter X, trigCut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+ H_{2}O MultiPath fitter Y, trigCut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+ H_{2}O MultiPath fitter Z, trigCut", 2000,-9000,9000);
   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "all fitted time, trigCut", 800,0,800);
   TH1F* hDeltaX_trig = new TH1F("hDeltaX_trig", "MP fitted X - mc X, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaY_trig = new TH1F("hDeltaY_trig", "MP fitted Y - mc Y, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaZ_trig = new TH1F("hDeltaZ_trig", "MP fitted Z - mc Z, FECD cut, trigCut", 2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "MP fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);


   TH1F* hiterN= new TH1F("hiterN", "nstart", 200,0,200);

   TVector3 u_fit,pos_fit, pos_mc;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   Double_t rPSUP = 8390;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   Double_t energy, wavelength ;
 
   int trigWord = 2;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   int nhitCut = 15;
   unsigned int fecdID = 9207;
   size_t evtNum = dsReader.GetEntryCount(); 
   for( size_t iEntry = 0; iEntry <evtNum; iEntry++ )
   {
     if(iEntry%100 == 0) std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     Int_t nevC = rDS.GetEVCount();
     const RAT::DS::MC& rmc= rDS.GetMC();
     const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
     pos_mc =rmcparticle.GetPosition();

     if(nevC>1) nevC = 1;
     for(Int_t iev=0;iev<nevC; iev++) {
     /// Get the Event information    
     const RAT::DS::EV& ev = rDS.GetEV(iev);
     int trig = ev.GetTrigType();
     if(ev.FitResultExists("multipath")) {
       try{
         RAT::DS::FitVertex fVertex = ev.GetFitResult("multipath").GetVertex(0);
         std::string FOMpar1 = ev.GetFitResult("multipath").GetFOMNames()[2];
         double nstart  =  ev.GetFitResult("multipath").GetFOM("starts");

         if(ev.GetFitResult("multipath").GetValid()) {

//[0]--"goodfits" [1]--"multipath_scintwater" [2]--"starts"

           //hiterN->Fill(nstart);
           TVector3 pos_fit = fVertex.GetPosition();
           double tfit = fVertex.GetTime();
           hFitTime->Fill(tfit);
           double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
           double mcposX=(pos_mc.X());double mcposY=(pos_mc.Y());double mcposZ=(pos_mc.Z());
           hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
           hDeltaX->Fill(posX - mcposX); hDeltaY->Fill(posY - mcposY); hDeltaZ->Fill(posZ - mcposZ);
           hfitRZ->Fill(sqrt(posX*posX+posY*posY),posZ);
           if(trig == 95) {//95
             hfitX_trig->Fill(posX); hfitY_trig->Fill(posY); hfitZ_trig->Fill(posZ);
             hDeltaX_trig->Fill(posX - mcposX); hDeltaY_trig->Fill(posY - mcposY); hDeltaZ_trig->Fill(posZ - mcposZ);
             hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);
             hFitTime_trig->Fill(tfit);
           }
        }
       } catch(exception& e) {std::cout<<e.what()<<"fit position is not valid"<<std::endl;}
    }
   }
  }// loop evt

  TString newfilename = "ResolMPnew_" + TString(filename); 
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

  hfitX->Write();hfitY->Write();hfitZ->Write();hDeltaX->Write();hDeltaY->Write();hDeltaZ->Write();
  hfitRZ->Write();hFitTime->Write();
  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hDeltaX_trig->Write();hDeltaY_trig->Write();hDeltaZ_trig->Write();
  hfitRZ_trig->Write();hFitTime_trig->Write();
  hiterN->Write();
  fp->Close();
}
