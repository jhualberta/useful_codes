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
// This code is for solar nu directional study in scintillator

//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm

using namespace std ;
//gStyle->SetPalette(53);
void analyMPW()
{  
   const char *filename = "FitMP_N16_rat6168_snoplus_water_1e4evts.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   bool trigCut = true;//false;

   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   RAT::DU::LightPathCalculator lightPathMC = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path

   const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   const double c_light = 299.792458;
   double grV1 = c_light/1.38486;
   double grV2 = c_light/1.392;

   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hHitTime_modeCut = new TH1F("hHitTime_modeCut", "hit time, modeCut", 800,0,800);

   TH1F* hTheta = new TH1F("hcosTheta", "SNO+ MPW mcDirec*fitDirec", 1000,-1,1);
   TH1F* hfitX = new TH1F("hfitX", "SNO+ multipath X", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ multipath Y", 2000,-9000,9000);   
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ multipath Z", 2000,-9000,9000);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);
   TH1F* hDeltaX = new TH1F("hDeltaX", "multipath X - mc X, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaY = new TH1F("hDeltaY", "multipath Y - mc Y, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaZ = new TH1F("hDeltaZ", "multipath Z - mc Z, FECD cut", 2000,-9000,9000);
   TH2F* hfitRZ = new TH2F("hfitRZ", "multipath R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+ H_{2}O multipath X, trigCut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+ H_{2}O multipath Y, trigCut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+ H_{2}O multipath Z, trigCut", 2000,-9000,9000);
   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "all fitted time, trigCut", 800,0,800);
   TH1F* hDeltaX_trig = new TH1F("hDeltaX_trig", "multipath X - mc X, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaY_trig = new TH1F("hDeltaY_trig", "multipath Y - mc Y, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaZ_trig = new TH1F("hDeltaZ_trig", "multipath Z - mc Z, FECD cut, trigCut", 2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "multipath R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH1F* htRes_trig = new TH1F("htRes_trig","time residual",400,-100,300);
   TH1F* htResMC = new TH1F("htResMC","time residual, mc",400,-100,300);
   TH1F* htResMC_straight1 = new TH1F("htResMC_straight1", "n = 1.378", 400,-100,300);
   TH1F* htResMC_straight2 = new TH1F("htResMC_straight2", "n = 1.38486", 400,-100,300);
   TH1F* htResMC_straight3 = new TH1F("htResMC_straight3", "n = 1.392", 400,-100,300);

   TH1F* htRes_straight1 = new TH1F("htRes_straight1","n = 1.378",400,-100,300);
   TH1F* htRes_straight2 = new TH1F("htRes_straight2","n = 1.38486",400,-100,300);
   TH1F* htRes_straight3 = new TH1F("htRes_straight3","n = 1.392",400,-100,300);

   TH1F* hcosTheta = new TH1F("hcosThetaPMT","angular distribution of pmt",200,-1,1);
   TH2F* hPMTphiTheta = new TH2F("hPMTphiTheta","phi theta distribution of pmt",628, -TMath::Pi(), TMath::Pi(), 200,-1,1);
   
   TH1F* hcosTheta_cut = new TH1F("hcosThetaPMT_cut","angular distribution of pmt",200,-1,1);
   TH2F* hPMTphiTheta_cut = new TH2F("hPMTphiTheta_cut","phi theta distribution of pmt",628, -TMath::Pi(), TMath::Pi(), 200,-1,1);
   TH1F* hcostTheta_fitCut = new TH1F("hcosThetaPMT_fitCut","angular distribution of pmt",200,-1,1);
   TH2F* htResVsCosTheta = new TH2F("htResVsCosTheta","tRes vs cosTheta",400,-100,300,200,-1,1);

   TVector3 u_mc, u_fit,pos_fit, pos_mc;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   Double_t rPSUP = 8900;
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
     pos_mc = rmcparticle.GetPosition();
     u_mc = rmcparticle.GetMomentum().Unit();
     if(nevC>1) nevC = 1;
     for(Int_t iev=0;iev<nevC; iev++) {
     /// Get the Event information    
     const RAT::DS::EV& ev = rDS.GetEV(iev);
     int trig = ev.GetTrigType();
     const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
 
     if(ev.FitResultExists("multipath")) {
       try{
         RAT::DS::FitVertex fVertex = ev.GetFitResult("multipath").GetVertex(0);
         if(ev.GetFitResult("multipath").GetValid()) 
	 {
           TVector3 pos_fit = fVertex.GetPosition();
           RAT::DS::FitVertex fVertex1 =ev.GetFitResult("multipathdirection").GetVertex(0);
           try{
             u_fit = fVertex1.GetDirection();
           }
           catch(exception& e)
           {std::cout<<e.what()<<" problems in dir fit"<<std::endl;}

           double tfit = fVertex.GetTime();
           hFitTime->Fill(tfit);
           double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
           double mcposX=(pos_mc.X());double mcposY=(pos_mc.Y());double mcposZ=(pos_mc.Z());
           hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
           hDeltaX->Fill(posX - mcposX); hDeltaY->Fill(posY - mcposY); hDeltaZ->Fill(posZ - mcposZ);
           hfitRZ->Fill(sqrt(posX*posX+posY*posY),posZ);
           //if(trig == 31 ) //TeLoaded trigType == 31; partialScint == 95
           {
             hfitX_trig->Fill(posX); hfitY_trig->Fill(posY); hfitZ_trig->Fill(posZ);
             hDeltaX_trig->Fill(posX - mcposX); hDeltaY_trig->Fill(posY - mcposY); hDeltaZ_trig->Fill(posZ - mcposZ);
             hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);
             hFitTime_trig->Fill(tfit);
             // hTheta->Fill(u_mc*u_fit);

             /// mode time cut, afterwards
             TH1F *hmode = new TH1F("hmode","",800,0,800);
             for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
             {
               double pmtTime =(calpmts.GetPMT(ipmt)).GetTime();
               hmode->Fill(pmtTime);
             }
             int binmax = hmode->GetMaximumBin();
             double modeTime = hmode->GetXaxis()->GetBinCenter(binmax);
              
             delete hmode;
             for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
             {
               TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
               double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
               hHitTime->Fill(hitTime);
               if( hitTime<modeTime+100 && hitTime>modeTime-50 ) hHitTime_modeCut->Fill(hitTime);
               lightPath.CalcByPosition( pos_fit, pmtpos );
               double distInInnerAV = lightPath.GetDistInInnerAV();
               double distInAV = lightPath.GetDistInAV();
               double distInWater = lightPath.GetDistInWater();
               const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
               double tRes = hitTime - transitTime - tfit; 
               double tRes_straight1 = hitTime - (pmtpos-pos_fit).Mag()/grVelocity - tfit;
               double tRes_straight2 = hitTime - (pmtpos-pos_fit).Mag()/grV1 - tfit;
               double tRes_straight3 = hitTime - (pmtpos-pos_fit).Mag()/grV2 - tfit;

               lightPathMC.CalcByPosition( pos_mc, pmtpos );
               double distInInnerAV_mc = lightPathMC.GetDistInInnerAV();
               double distInAV_mc = lightPathMC.GetDistInAV();
               double distInWater_mc = lightPathMC.GetDistInWater();
               const double transitTime2 = groupVelocity.CalcByDistance( distInInnerAV_mc, distInAV_mc, distInWater_mc );
               double tResMC = hitTime - transitTime2 - 390 + rDS.GetMCEV(iev).GetGTTime();
               double tResMC_straight1 = hitTime - (pmtpos-pos_mc).Mag()/grVelocity - 390 + rDS.GetMCEV(iev).GetGTTime();
               double tResMC_straight2 = hitTime - (pmtpos-pos_mc).Mag()/grV1 - 390 + rDS.GetMCEV(iev).GetGTTime();
               double tResMC_straight3 = hitTime - (pmtpos-pos_mc).Mag()/grV2 - 390 + rDS.GetMCEV(iev).GetGTTime();

               htResMC->Fill(tResMC);
               htResMC_straight1->Fill(tResMC_straight1);
               htResMC_straight2->Fill(tResMC_straight2);
               htResMC_straight3->Fill(tResMC_straight3);

	       htRes_trig->Fill(tRes);
               htRes_straight1->Fill(tRes_straight1);
               htRes_straight2->Fill(tRes_straight2);
               htRes_straight3->Fill(tRes_straight3);

               TVector3 Xdiff = (pmtpos - pos_mc).Unit();
               hcosTheta->Fill(Xdiff*u_mc);
               hPMTphiTheta->Fill(pmtpos.Phi(), pmtpos.CosTheta());
               //if(hitTime>215 && hitTime<219)
               if( hitTime<modeTime+100 && hitTime>modeTime-50 )
               {

                 if(tResMC<-1 && tResMC>-5.5)
                 hcosTheta_cut->Fill(((pmtpos - pos_mc).Unit())*u_mc); 

                 if(-6<tRes && tRes<-1) {
                   hPMTphiTheta_cut->Fill(pmtpos.Phi(), pmtpos.CosTheta());
                   hcostTheta_fitCut->Fill(((pmtpos - pos_fit).Unit())*u_mc);
                 }
               } // modeCut
              }//for trig event, loop pmts and play with different tRes windows
            } // trigger word
           } // fitValid 
       } catch(exception& e) {std::cout<<e.what()<<"fit position is not valid"<<std::endl;}
     } // fit exists
     } // if exist triggered events
  }// loop evt

  TString newfilename = "GetPMT_lot_" + TString(filename); 
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
  hfitRZ->Write();hTheta->Write();hcosTheta->Write();hPMTphiTheta->Write();
  hHitTime->Write();hHitTime_modeCut->Write();hFitTime->Write();hFitTime_trig->Write();
  htRes_trig->Write();htResMC->Write();
  htResMC_straight1->Write();
  htResMC_straight2->Write();
  htResMC_straight3->Write();

  htRes_straight1->Write();
  htRes_straight2->Write();
  htRes_straight3->Write();

  hcosTheta_cut->Write();hcostTheta_fitCut->Write();hPMTphiTheta_cut->Write();
  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hDeltaX_trig->Write();hDeltaY_trig->Write();hDeltaZ_trig->Write();
  hfitRZ_trig->Write();
  fp->Close();
}
