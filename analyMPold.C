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
fit10 = MPW
fit9 = FTP
*/
using namespace std ;
//gStyle->SetPalette(53);
void analyMPold()
{ 
  // change gaus fit range
   const char *filename = "FitMP_MC_realTrig106904_2p2MeVgamma_IsoCenter_1e4evt.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   bool trigCut = true;
   double driveCor_p0 = 0.995765; // drive parameters in MPW; see Water-Phase-Unidoc; x_correct = x_fit*p0+u_fit*p1.   ,
   double driveCor_p1 = -63.826;//-38.8275;
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TH1F* hTrig = new TH1F("HTrig","",100,0,100);
   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hNhit = new TH1F("hNhit","total nhits",100,0,100);   
   TH2F* hNhitsVsPosMag = new TH2F("hNhitsVsPosR","nhits vs posMag",9000,0,9000, 300,0,3000);
   TH1F* hNhitSelect = new TH1F("hNhitSelect","total nhits",100,0,100);

   TH1F* hFOM = new TH1F("hFOM", "FOM MPW", 1000,0,1000);
   TH1F* hFOM_nhitScaled = new TH1F("hFOM_nhitScaled", "FOM scaled to nhit", 200,0,20);
   TH1F* hFOM_nhitSelect = new TH1F("hFOM_nhitSelect", "FOM scaled to selected nhit", 200,0,20);
   TH2F* hFOMvsNhit = new TH2F("hFOMvsNhit","FOM vs Nhit",100,0,100,1000,0,1000);

   TH1F* hfitX = new TH1F("hfitX", "SNO+ H_{2}O MPW fitter X", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+ H_{2}O MPW fitter Y", 2000,-9000,9000);   
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+ H_{2}O MPW fitter Z", 2000,-9000,9000);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 800,0,800);

   /// only for MC comparisons
   TH1F* hDeltaX = new TH1F("hDeltaX", "MPW fitted X - mc X, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaY = new TH1F("hDeltaY", "MPW fitted Y - mc Y, FECD cut", 2000,-9000,9000);
   TH1F* hDeltaZ = new TH1F("hDeltaZ", "MPW fitted Z - mc Z, FECD cut", 2000,-9000,9000);

   TH2F* hDeltaXvsNhits = new TH2F("hDeltaXsNhits", "MPW fitted X - mc X vs nhits", 2000,-9000,9000,100,0,100);
   TH2F* hDeltaYvsNhits = new TH2F("hDeltaYsNhits", "MPW fitted Y - mc Y vs nhits", 2000,-9000,9000,100,0,100);
   TH2F* hDeltaZvsNhits = new TH2F("hDeltaZsNhits", "MPW fitted Z - mc Z vs nhits", 2000,-9000,9000,100,0,100);

   TH2F* hfitRZ = new TH2F("hfitRZ", "MPW fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+ H_{2}O MPW fitter X, trigCut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+ H_{2}O MPW fitter Y, trigCut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+ H_{2}O MPW fitter Z, trigCut", 2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "MPW fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "all fitted time, trigCut", 800,0,800);

   /// only for MC comparisons
   TH1F* hDeltaX_trig = new TH1F("hDeltaX_trig", "MPW fitted X - mc X, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaY_trig = new TH1F("hDeltaY_trig", "MPW fitted Y - mc Y, FECD cut, trigCut", 2000,-9000,9000);
   TH1F* hDeltaZ_trig = new TH1F("hDeltaZ_trig", "MPW fitted Z - mc Z, FECD cut, trigCut", 2000,-9000,9000);
   TH2F* hDeltaXvsFomPos_trig = new TH2F("hDeltaXvsFomPos_trig", "MPW fitted X - mc X vs Fom of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,500);
   TH2F* hDeltaYvsFomPos_trig = new TH2F("hDeltaYvsFomPos_trig", "MPW fitted Y - mc Y vs Fom of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,500);
   TH2F* hDeltaZvsFomPos_trig = new TH2F("hDeltaZvsFomPos_trig", "MPW fitted Z - mc Z vs Fom of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,500);

   TH2F* hDeltaXvsFomPosScaled_trig = new TH2F("hDeltaXvsFomPosScaled_trig", "MPW fitted X - mc X vs Fom scaled of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,100);
   TH2F* hDeltaYvsFomPosScaled_trig = new TH2F("hDeltaYvsFomPosScaled_trig", "MPW fitted Y - mc Y vs Fom scaled of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,100);
   TH2F* hDeltaZvsFomPosScaled_trig = new TH2F("hDeltaZvsFomPosScaled_trig", "MPW fitted Z - mc Z vs Fom scaled of pos, FECD cut, trigCut", 2000,-9000,9000,100,0,100);

   TH1F* hEnergy = new TH1F("hEnergy", "SNO+ H_{2}O MPW fitter energy", 100,0,10);   

   TH1F* hiterN= new TH1F("hiterN", "nstart", 200,0,200);
   TH1F* hcosTheta = new TH1F("hcosTheta","u_mc*u_fit",200,-1,1);
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
   int nhitCut = 0;
   unsigned int fecdID = 9207;
   size_t evtNum = dsReader.GetEntryCount(); 
   for( size_t iEntry = 0; iEntry <evtNum; iEntry++ )
   {
     // if(iEntry%100 == 0) std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     Int_t nevC = rDS.GetEVCount();
     const RAT::DS::MC& rmc= rDS.GetMC();
     const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
     pos_mc =rmcparticle.GetPosition();
     TVector3 u_mc = rmcparticle.GetMomentum();
     u_mc = u_mc.Unit();
     std::string fitName = "multipath";
     for(Int_t iev=0;iev<nevC; iev++) {
     /// Get the Event information    
     const RAT::DS::EV& ev = rDS.GetEV(iev);
     int trig = ev.GetTrigType();
     double nhits = ev.GetNhitsCleaned();
     hNhit->Fill(nhits);
     if( nhits>nhitCut ) {
     if(ev.FitResultExists(fitName)) {
       try{
         RAT::DS::FitVertex fVertex = ev.GetFitResult(fitName).GetVertex(0);
         double fomValue  =  ev.GetFitResult(fitName).GetFOM("multipath_waterposition");
         double fomSelectedNhits = ev.GetFitResult(fitName).GetFOM("multipath_SelectedNHit_waterposition");
         double fomScale =  fomValue/fomSelectedNhits;//ev.GetFitResult(fitName).GetFOM("multipath_LogLSelectedNHit_waterposition");
         vector<string> fomname = ev.GetFitResult(fitName).GetFOMNames();
        // for(vector<string>::iterator i = fomname.begin(); i != fomname.end(); ++i) std::cout<<*i<<" ";
         //cout<<endl;
         /* fom pars
           Directiongoodfits Directionmultipath_waterdirection Directionstarts EnergyGtest EnergyMedianDev EnergyMedianDevHit EnergyMedianProb EnergyMedianProbHit EnergyRSP EnergyUtest EnergynCerPhotons EnergypromptHits Positiongoodfits Positionmultipath_waterposition Positionstarts
         */
       if(ev.GetFitResult(fitName).GetValid()) {
        if( fomSelectedNhits>=4 ) 
        {
         hNhitSelect->Fill(fomSelectedNhits);
         hFOM->Fill(fomValue); hFOMvsNhit->Fill(nhits,fomValue);
         hFOM_nhitScaled->Fill(fomScale);
         hFOM_nhitSelect->Fill(fomValue/fomSelectedNhits);

         TVector3 pos_fit = fVertex.GetPosition();
         double tfit = fVertex.GetTime();
         RAT::DS::FitVertex fVertex1 = ev.GetFitResult("multipathdirection").GetVertex(0);

         if(pos_fit.Mag()<8390 && fVertex.ValidPosition() && fVertex1.ValidDirection()) 
         {
           if(fomScale>8)
           {
           u_fit = fVertex1.GetDirection();
           pos_fit = driveCor_p0*pos_fit+driveCor_p1*u_fit;// drive correction afterwards
           double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
           double mcposX=(pos_mc.X());double mcposY=(pos_mc.Y());double mcposZ=(pos_mc.Z());
           hfitX->Fill(posX); hfitY->Fill(posY); hfitZ->Fill(posZ);
           hDeltaX->Fill(posX - mcposX); hDeltaY->Fill(posY - mcposY); hDeltaZ->Fill(posZ - mcposZ);
           hfitRZ->Fill(sqrt(posX*posX+posY*posY),posZ);
           hFitTime->Fill(tfit);
           //if(trig == 2 ) 
           {
             hfitX_trig->Fill(posX); hfitY_trig->Fill(posY); hfitZ_trig->Fill(posZ);
             hDeltaX_trig->Fill(posX - mcposX); hDeltaY_trig->Fill(posY - mcposY); hDeltaZ_trig->Fill(posZ - mcposZ);
             hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);
             hFitTime_trig->Fill(tfit);
             hcosTheta->Fill(u_mc*u_fit);
             hDeltaXvsFomPos_trig->Fill(posX-mcposX, fomValue);
             hDeltaYvsFomPos_trig->Fill(posY-mcposY, fomValue);
             hDeltaZvsFomPos_trig->Fill(posZ-mcposZ, fomValue);
             hDeltaXvsFomPosScaled_trig->Fill(posX-mcposX, fomScale);//fomValue/nhits);
             hDeltaYvsFomPosScaled_trig->Fill(posY-mcposY, fomScale);//fomValue/nhits);
             hDeltaZvsFomPosScaled_trig->Fill(posZ-mcposZ, fomScale);//fomValue/nhits);

             hDeltaXvsNhits->Fill(posX-mcposX,nhits);
             hDeltaYvsNhits->Fill(posY-mcposY,nhits);
             hDeltaZvsNhits->Fill(posZ-mcposZ,nhits);
           }
           }
        } // if pos Valid
      } // NhitSelected>4
     } // global Valid
     } catch(exception& e) {std::cout<<e.what()<<"fit position is not valid"<<std::endl;}
    }
   }
  } // nhits
  }// loop evt

  TString newfilename = "ResolMPWoldFomScale_" + TString(filename); 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  
  double fitRange = 1000;
  TF1 *gx = new TF1("gx","gaus",hfitX_trig->GetMean()-fitRange,hfitX_trig->GetMean()+fitRange);
  TF1 *gy = new TF1("gy","gaus",hfitY_trig->GetMean()-fitRange,hfitY_trig->GetMean()+fitRange);
  TF1 *gz = new TF1("gz","gaus",hfitZ_trig->GetMean()-fitRange,hfitZ_trig->GetMean()+fitRange);
  gx->SetLineColor(kRed);gy->SetLineColor(kRed);gz->SetLineColor(kRed);
  if(!trigCut) {
    hDeltaX_trig->Fit(gx);
    hDeltaY_trig->Fit(gy);
    hDeltaZ_trig->Fit(gz);
  }
  else {
    hDeltaX->Fit(gx,"R");
    hDeltaY->Fit(gy,"R");
    hDeltaZ->Fit(gz,"R");
  }
  //cout<<gx->GetParameters(0)<<endl;
  hNhit->Write(); hNhitSelect->Write(); hFOM_nhitSelect->Write();
  hfitX->Write();hfitY->Write();hfitZ->Write();hDeltaX->Write();hDeltaY->Write();hDeltaZ->Write();
  hfitRZ->Write();hFitTime->Write();
  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hDeltaX_trig->Write();hDeltaY_trig->Write();hDeltaZ_trig->Write();
  hfitRZ_trig->Write();hFitTime_trig->Write();
  //hiterN->Write();
  hFOM->Write(); hFOMvsNhit->Write();
  hFOM_nhitScaled->Write();
  hcosTheta->Write();
  hEnergy->Write();
  hDeltaXvsFomPos_trig->Write();
  hDeltaYvsFomPos_trig->Write();
  hDeltaZvsFomPos_trig->Write();
  hDeltaXvsNhits->Write();hDeltaYvsNhits->Write();hDeltaZvsNhits->Write();
  hDeltaXvsFomPosScaled_trig->Write();
  hDeltaYvsFomPosScaled_trig->Write();
  hDeltaZvsFomPosScaled_trig->Write();
 
  fp->Close();
}
