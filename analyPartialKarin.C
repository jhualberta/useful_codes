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
//!!! for SNO+ 16N data
//Unit: mm

using namespace std ;
void analyPartialKarin()
{ 
//   TVector3 srcPos(-1120.8, 1041.4, 6172.5);    	
   TVector3 srcPos(-1120.8,1041.4,6108.0);
   const char* filename = "FitPartial_SNOP_0000252308_000.root";
   //Analysis10_r0000250425_s000_p001.root
   double driveP0 = 0.995765, driveP1 = -63.826;
   double distCut = 500;
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
//   const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();

   TH1F* hTrig = new TH1F("HTrig","",100,0,100);

   TH1F* hLogL = new TH1F("hLogL", "LogL", 1000,0,1000); 
   TH2F* hLogLvsNhit = new TH2F("hLogLvsNhit","LogL vs Nhit",1000,0,1000,1000,0,10000);
   TH2F* hLogLvsFitT = new TH2F("hLogLvsFitT","LogL vs Fitted Time",300,0,300,1000,0,10000); 
   TH2F* hLogLvsFitX = new TH2F("hLogLvsFitX","LogL vs Fitted X",2000,-6000,6000,1000,0,10000);
   TH2F* hLogLvsFitY = new TH2F("hLogLvsFitY","LogL vs Fitted Y",2000,-6000,6000,1000,0,10000);
   TH2F* hLogLvsFitZ = new TH2F("hLogLvsFitZ","LogL vs Fitted Z",2000,-6000,6000,1000,0,10000);

   TH1F* hITR = new TH1F("hITR", "itr", 100,0,1);
   TH1F* hITR_scint = new TH1F("hITR_scint", "itr in scint", 100,0,1);
   TH1F* hITR_water = new TH1F("hITR_water", "itr in water", 100,0,1);
   TH1F* hITR_exwater = new TH1F("hITR_exwater", "itr in exwater", 100,0,1);
   TH1F* htResITRlt0p3 = new TH1F("htResITRlt0p3", "Time Residue,itr<0.3", 400,-100,300);
   TH1F* htResITRht0p3 = new TH1F("htResITRht0p3", "Time Residue,itr>0.3", 400,-100,300);
   TH1F* htResITR0p3 = new TH1F("htResITR0p3", "Time Residue, |itr-0.3|<0.1", 400,-100,300);
   TH2F* hfitRhoZ_ITRlt0p3 = new TH2F("hfitRhoZ_ITRlt0p3", "SNO+  fitted Rho vs Z, itr<0.3", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRZ_ITRlt0p3 = new TH2F("hfitRZ_ITRlt0p3", "SNO+  fitted R vs Z, itr<0.3", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRhoZ_ITRht0p3 = new TH2F("hfitRhoZ_ITRht0p3", "SNO+  fitted Rho vs Z, itr>0.3", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRZ_ITRht0p3 = new TH2F("hfitRZ_ITRht0p3", "SNO+  fitted R vs Z, itr>0.3", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRhoZ_ITR0p3 = new TH2F("hfitRhoZ_ITR0p3", "SNO+  fitted Rho vs Z, |itr-0.3|<0.1", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRZ_ITR0p3 = new TH2F("hfitRZ_ITR0p3", "SNO+  fitted R vs Z, |itr-0.3|<0.1", 1000, 0,9000,2000,-9000,9000);

   TH1F* hscaleLogL = new TH1F("hscaleLogL", "scaledLogL", 1000,0,1000);
   TH2F* hscaleLogLvsNhit = new TH2F("hscaleLogLvsNhit","scaledLogL vs Nhit",1000,0,1000,100,0,100);
   TH2F* hscaleLogLvsFitT = new TH2F("hscaleLogLvsFitT","scaledLogL vs Fitted Time",300,0,300,100,0,100);
   TH2F* hscaleLogLvsFitX = new TH2F("hscaleLogLvsFitX","scaledLogL vs Fitted X",2000,-6000,6000,100,0,100);
   TH2F* hscaleLogLvsFitY = new TH2F("hscaleLogLvsFitY","scaledLogL vs Fitted Y",2000,-6000,6000,100,0,100);
   TH2F* hscaleLogLvsFitZ = new TH2F("hscaleLogLvsFitZ","scaledLogL vs Fitted Z",2000,-6000,6000,100,0,100);

   TH1F* hLogL_leonCut = new TH1F("hLogL_leonCut", "LogL, leonCut", 1000,0,1000);
   TH2F* hLogLvsNhit_leonCut = new TH2F("hLogLvsNhit_leonCut","LogL vs Nhit",1000,0,1000,10000,0,10000);
   TH2F* hLogLvsFitT_leonCut = new TH2F("hLogLvsFitT_leonCut","LogL vs Fitted Time, leonCut",300,0,300,1000,0,10000);
   TH2F* hLogLvsFitX_leonCut = new TH2F("hLogLvsFitX_leonCut","LogL vs Fitted X, leonCut",2000,-6000,6000,1000,0,10000);
   TH2F* hLogLvsFitY_leonCut = new TH2F("hLogLvsFitY_leonCut","LogL vs Fitted Y, leonCut",2000,-6000,6000,1000,0,10000);
   TH2F* hLogLvsFitZ_leonCut = new TH2F("hLogLvsFitZ_leonCut","LogL vs Fitted Z, leonCut",2000,-6000,6000,1000,0,10000);

   TH1F* hscaleLogL_leonCut = new TH1F("hscaleLogL_leonCut", "scaledLogL", 100,0,100);
   TH2F* hscaleLogLvsNhit_leonCut = new TH2F("scalehLogLvsNhit","LogL vs Nhit",1000,0,1000,100,0,100);
   TH2F* hscaleLogLvsFitT_leonCut = new TH2F("hscaleLogLvsFitT_leonCut","scaledLogL vs Fitted Time, leonCut",300,0,300,100,0,100);
   TH2F* hscaleLogLvsFitX_leonCut = new TH2F("hscaleLogLvsFitX_leonCut","scaledLogL vs Fitted X, leonCut",2000,-6000,6000,100,0,100);
   TH2F* hscaleLogLvsFitY_leonCut = new TH2F("hscaleLogLvsFitY_leonCut","scaledLogL vs Fitted Y, leonCut",2000,-6000,6000,100,0,100);
   TH2F* hscaleLogLvsFitZ_leonCut = new TH2F("hscaleLogLvsFitZ_leonCut","scaledLogL vs Fitted Z, leonCut",2000,-6000,6000,100,0,100);

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
   TH1F* htRes = new TH1F("htRes", "Time Residue, using lightPath, trig", 400,-100,300);
   TH1F* htRes_nhit100 = new TH1F("htRes_nhit100", "Time Residue, using lightPath, trig", 400,-100,300);
   TH1F* htRes_nhit300 = new TH1F("htRes_nhit300", "Time Residue, using lightPath, trig", 400,-100,300);
   TH2F* htResVsPMTz = new TH2F("htResVsPMTz", "Time Residue vs PMTz, using lightPath", 400,-100,300, 2000,-9000,9000);
   TH2F* htResVsPMTz_nhit100 = new TH2F("htResVsPMTz_nhit100", "Time Residue vs PMTz, using lightPath, nhit>100", 400,-100,300, 2000,-9000,9000);

   TH1F* htResInScint = new TH1F("htResInScint", "Time Residue, using lightPath, in scint", 400,-100,300);
   TH2F* htResInScintVsPMTz = new TH2F("htResInScintVsPMTz", "Time Residue vs PMTz, using lightPath, in scint", 400,-100,300, 2000,-9000,9000);
   TH1F* htResInScint_nhit100 = new TH1F("htResInScint_nhit100", "Time Residue, using lightPath, in scint, nhit>100", 400,-100,300);
   TH2F* htResInScintVsPMTz_nhit100 = new TH2F("htResInScintVsPMTz_nhit100", "Time Residue vs PMTz, using lightPath, in scint, nhit>100", 400,-100,300, 2000,-9000,9000);
   TH1F* htResInScint_nhit300 = new TH1F("htResInScint_nhit300", "Time Residue, using lightPath, in scint, nhit>300", 400,-100,300);

   TH1F* htResInWater = new TH1F("htResInWater", "Time Residue, using lightPath, in water", 400,-100,300);
   TH2F* htResInWaterVsPMTz = new TH2F("htResInWaterVsPMTz", "Time Residue vs PMTz, using lightPath, in water", 400,-100,300, 2000,-9000,9000);
   TH1F* htResInWater_nhit100 = new TH1F("htResInWater_nhit100", "Time Residue, using lightPath, in water, nhit<100", 400,-100,300);
   TH2F* htResInWaterVsPMTz_nhit100 = new TH2F("htResInWaterVsPMTz_nhit100", "Time Residue vs PMTz, using lightPath, in water, nhit>100", 400,-100,300, 2000,-9000,9000);
   TH1F* htResInWater_nhit300 = new TH1F("htResInWater_nhit300", "Time Residue, using lightPath, in water, nhit<300", 400,-100,300);

   TH2F* htResVsX = new TH2F("htResVsX","tRes Vs X",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsX_nhit100 = new TH2F("htResVsX_nhit100","tRes Vs X, nhit100",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsX_nhit300 = new TH2F("htResVsX_nhit300","tRes Vs X, nhit300",400, -100, 300,2000,-9000,9000);

   TH2F* htResVsY = new TH2F("htResVsY","tRes Vs Y",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsY_nhit100 = new TH2F("htResVsY_nhit100","tRes Vs Y, nhit100",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsY_nhit300 = new TH2F("htResVsY_nhit300","tRes Vs Y, nhit300",400, -100, 300,2000,-9000,9000);

   TH2F* htResVsZ = new TH2F("htResVsZ","tRes Vs Z",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsZ_nhit100 = new TH2F("htResVsZ_nhit100","tRes Vs Z, nhit100",400, -100, 300,2000,-9000,9000);
   TH2F* htResVsZ_nhit300 = new TH2F("htResVsZ_nhit300","tRes Vs Z, nhit300",400, -100, 300,2000,-9000,9000);

   TH2F* hNhitsVsX = new TH2F("hNhitsVsX","nhits vs x", 2000, 0, 2000, 2000, -9000, 9000);
   TH2F* hNhitsVsY = new TH2F("hNhitsVsY","nhits vs y", 2000, 0, 2000, 2000, -9000, 9000);
   TH2F* hNhitsVsZ = new TH2F("hNhitsVsZ","nhits vs z", 2000, 0, 2000, 2000, -9000, 9000);
   TH2F* hNhitsVsR = new TH2F("hNhitsVsR","nhits vs r", 2000, 0, 2000, 1000, 0, 9000);
  
   // use MPW straight light path
   TH1F* htRes_1p38486 = new TH1F("htRes_1p38486", "Time Residue, default = c/1.38486, trig", 400,-100,300);
   //TH1F* htRes_rmOdd = new TH1F("htRes_rmOdd", "Time Residue, default = c/1.38486, remove odd PMTs", 1600,-100,300);
   TH1F* htRes_1p41 = new TH1F("htRes_1p41", "Time Residue, c/1.40, trig", 400,-100,300);
   TH1F* htRes_1p63 = new TH1F("htRes_1p63", "Time Residue, c/1.378, trig", 400,-100,300);
   TH1F* hSrcToEvt = new TH1F("hSrcToEvt","|Xsrc-Xevt|",1000,0,9000);
   // Check PMTs, distanceCut>0.5 m
   TH1F* hcosTheta = new TH1F("hcosTheta","srcToPos*posToPMT",200,-1,1);
   TH1F* hcosTheta_rmOdd = new TH1F("hcosTheta_rmOdd","srcToPos*posToPMT, !(5900<posZcor<6100)",200,-1,1);
   TH1F* hcosTheta_rmOdd1 = new TH1F("hcosTheta_rmOdd1","srcToPos*posToPMT, !(5900<posZcor<6100 || 5090<posZcor<5110)",200,-1,1);
   TH2F* hcosThetaVsNhits = new TH2F("hcosThetaVsNhits", "Angle between Directions;Cos(#theta_{PMT}); Nhits, fecd == 9188", 200, -1.0, 1.0, 200, 0, 2000);
   TH1F* hcosTheta_prompt = new TH1F("hcosTheta_prompt","srcToPos*posToPMT",200,-1,1);
   TH2F* hcosThetaVsNhits_prompt = new TH2F("hcosThetaVsNhits_prompt", "Angle between Directions;Cos(#theta_{PMT}); Nhits, fecd == 9188", 200, -1.0, 1.0, 200, 0, 2000);
   TH2F* hcosThetaVsTres = new TH2F("hcosThetaVsTres", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188", 200, -1.0, 1.0, 400, -100, 300);
   TH2F* hcosThetaVsTres_nhit100 = new TH2F("hcosThetaVsTres_nhit100", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188, nhits>100", 200, -1.0, 1.0, 400, -100, 300);
   TH2F* hcosThetaVsTres_nhitLess100 = new TH2F("hcosThetaVsTres_nhitLess100", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188, nhits<100", 200, -1.0, 1.0, 400, -100, 300);


   TH1F* hcosTheta_nhit100 = new TH1F("hcosTheta_nhit100","srcToPos*posToPMT",200,-1,1);
   TH1F* hcosTheta_prompt_nhit100 = new TH1F("hcosTheta_prompt_nhit100","srcToPos*posToPMT, nhit>100",200,-1,1);

   TH1F* hcosTheta_nhitLess100 = new TH1F("hcosTheta_nhitLess100","srcToPos*posToPMT, Nhit<100",200,-1,1);
   TH1F* hcosTheta_prompt_nhitLess100 = new TH1F("hcosTheta_prompt_nhitLess100","srcToPos*posToPMT, nhit<100",200,-1,1);

   TH1F* hPMTid = new TH1F("hPMTid", "all PMT id", 10000, 0, 10000);
   TH1F* hPMTid_odd = new TH1F("hPMTid_odd", "PMT id in [-10,-5]ns ", 10000, 0, 10000);
   TH2F* hPMTthetaPhi = new TH2F("hPMTthetaPhi", "PMT distributions", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1); 
   TH2F* hPMTthetaPhi_odd = new TH2F("hPMTthetaPhi_odd", "PMT distributions in [-10,-5] ns", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1);
   int fNBinsPhi = 45, fNBinsTheta = 45;
   TH2F* hPMTthetaPhi_sinu = new TH2F("hPMTthetaPhi_sinu", "PMT distributions, sinusoidal ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi()); 
   TH2F* hPMTthetaPhi_sinu_nhit100 = new TH2F("hPMTthetaPhi_sinu_nhit100", "PMT distributions, sinusoidal, nhit>100 ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());
   TH2F* hPMTthetaPhi_sinu_nhit300 = new TH2F("hPMTthetaPhi_sinu_nhit300", "PMT distributions, sinusoidal, nhit>300 ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());
   TH2F* hPMTthetaPhi_sinu_leonCut = new TH2F("hPMTthetaPhi_sinu_leonCut", "PMT distributions, sinusoidal,leonCut ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi()); 
   TH2F* hPMTthetaPhi_prompt = new TH2F("hPMTthetaPhi_prompt", "PMT distributions, prompt cut", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1);
   TH2F* hPMTthetaPhi_sinu_prompt = new TH2F("hPMTthetaPhi_sinu_prompt", "PMT distributions, sinusoidal, prompt cut ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());

   TH2F* hPMTthetaPhi_prompt_nhit100 = new TH2F("hPMTthetaPhi_prompt_nhit100", "PMT distributions, prompt cut, nhits>100", 628, -TMath::Pi(), TMath::Pi(), 200,-1,1);
   TH2F* hPMTthetaPhi_sinu_prompt_nhit100 = new TH2F("hPMTthetaPhi_sinu_prompt_nhit100", "PMT distributions, sinusoidal, prompt cut, nhits>100 ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());

   TH1F* hfitX = new TH1F("hfitX", "SNO+  fitted X, FECD cut", 2000,-9000,9000);   
   TH1F* hfitY = new TH1F("hfitY", "SNO+  fitted Y, FECD cut", 2000,-9000,9000);  
   TH1F* hfitZ = new TH1F("hfitZ", "SNO+  fitted Z, FECD cut", 2000,-9000,9000);
   TH2F* hfitXY = new TH2F("hfitXY", "SNO+  fitted X vs Y, FECD cut", 2000,-9000,9000,2000,-9000,9000); 
   TH2F* hfitXZ = new TH2F("hfitXZ", "SNO+  fitted X vs Z, FECD cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRhoZ = new TH2F("hfitRhoZ", "SNO+  fitted Rho vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRZ = new TH2F("hfitRZ", "SNO+  fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH2F* hfitXY_nhit100 = new TH2F("hfitXY_nhit100", "SNO+  fitted X vs Y, FECD cut, nhit>100", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_nhit100 = new TH2F("hfitXZ_nhit100", "SNO+  fitted X vs Z, FECD cut, nhit>100", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRhoZ_nhit100 = new TH2F("hfitRhoZ_nhit100", "SNO+  fitted Rho vs Z, FECD cut, nhit>100", 1000, 0,9000,2000,-9000,9000);
   TH2F* hfitRZ_nhit100 = new TH2F("hfitRZ_nhit100", "SNO+  fitted R vs Z, FECD cut, nhit>100", 1000, 0,9000,2000,-9000,9000);

   TH1F* hPosMag = new TH1F("hPosMag","SNO+  fitted Mag()",1000,0,10000);
   TH2F* hfitRZ_odd = new TH2F("hfitRZ_odd", "SNO+  fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);

   TH1F* hfitX_trig = new TH1F("hfitX_trig", "SNO+  fitted X, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitY_trig = new TH1F("hfitY_trig", "SNO+  fitted Y, FECD+trig cut", 2000,-9000,9000);
   TH1F* hfitZ_trig = new TH1F("hfitZ_trig", "SNO+  fitted Z, FECD+trig cut", 2000,-9000,9000);
   TH2F* hfitXY_trig = new TH2F("hfitXY_trig", "SNO+  fitted X vs Y, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitXZ_trig = new TH2F("hfitXZ_trig", "SNO+  fitted X vs Z, FECD+trig cut", 2000,-9000,9000,2000,-9000,9000);
   TH2F* hfitRZ_trig = new TH2F("hfitRZ_trig", "SNO+  fitted R vs Z, FECD+trig cut", 1000, 0,9000,2000,-9000,9000);
   TH1F* hPosMag_trig = new TH1F("hPosMag_trig","SNO+  fitted Mag(), FECD+trig cut",1000,0,10000);

   TH1F* hHitTime = new TH1F("hHitTime", "all hit time", 800,0,800);
   TH1F* hFitTime = new TH1F("hFitTime", "all fitted time", 500,0,500);
   TH1F* hFitTime_odd = new TH1F("hFitTime_odd", "odd fitted time", 500,0,500);

   TH1F* hFitTime_trig = new TH1F("hFitTime_trig", "fitted time,FECD cuts, triggered", 800,0,800); 
   TH2F* hFlatVsZ = new TH2F("hFlatVsZ","fitted (x^2+y^2)/6000^2 vs Z",5000,0,5,2000,-9000,9000);
   TH1F* hNhits = new TH1F("hNhits","Nhits",100,0,100); 
   TH1F* hNhitsFECD = new TH1F("hNhitsFECD","Nhits FECD",100,0,100);
   TH2F* hNhitsVsMag = new TH2F("hNhitsVsMag","pos.Mag() vs Nhits",1000,0,10000,100,0,100);
   TH2F* hNhitVsFlat = new TH2F("hNhitsVsFlat","Nhits vs  (x^2+y^2)/6000^2",100,0,100,5000,0,5);

   TH1F* hwaterlevel = new TH1F("hwaterlevel", "water level", 2000,4000,6000);
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
   double n2 = 1.41;
   double n3 = 1.63;

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
   int nhitCut = 200;//15;
   unsigned int fecdID = 9188;
   TString fitName = "partialFitter";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
   {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     Int_t nevC = rDS.GetEVCount();
    for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     int trig = rev.GetTrigType();//get triggerType
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     countNoCleanTotal++;
     if(trig&(trigWord)){countNoCleanTrig++;}
     //RAT::DS::DataQCFlags dcflags = rev.GetDataCleaningFlags(); //!! dataClean not used for calibration source
     //if( RAT::EventIsClean( rev, dcAnalysisWord ) )
     {
      countClean++;//count cleaned total events
      if(trig&(trigWord)){countTrig++;}
      double nhits = rev.GetNhitsCleaned();//NOTE: use nhits_cleaned
      hNhits->Fill(nhits);
      for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)
      {
	   if(calpmts.GetFECDPMT(ipmt).GetID()==fecdID){
		countFECDtotal++;
	     if(trig&(trigWord)){countTrigFECD++;}
	   }
      }//count total FECD events
      if(nhits>10) {
      if(rev.FitResultExists("partialFitter")){//find partialFitter exists
       try{
         RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
	 //if(rev.GetFitResult("partialFitter").GetValid())//Global Validity of the fitter !!! NOTE: this is different
         {
          //double fomValue =  rev.GetFitResult("partialFitter").GetFOM(FOMname);
          pos_fit=fitVertex.GetPosition(); 
          //cout<<rev.GetClassifierNames().size()<<endl;//3 classifiers
          //cout<<rev.GetClassifierNames()[0]<<" "<<rev.GetClassifierNames()[1]<<" "<<rev.GetClassifierNames()[2]<<endl;
          //ITR:partialFitter QPDT:partialFitter isotropy:partialFitter

	  double scaleLogL = rev.GetFitResult("partialFitter").GetFOM("Positionmultipath_scintwater")/rev.GetFitResult("partialFitter").GetFOM("Positionmultipath_SelectedNHit_scintwater");
          double LogL = rev.GetFitResult("partialFitter").GetFOM("Positionmultipath_scintwater");
          double itr = rev.GetClassifierResult("ITR:partialFitter").GetClassification("ITR");
         // RAT::DS::FitVertex fVertex1 = rev.GetFitResult("partialFitterdirection").GetVertex(0);
         // if(fVertex1.ValidDirection()) {
         //    try {
         //         u_fit = fVertex1.GetDirection();
         //         pos_cor = fVertex1.GetPosition(); //Note: here the position result is after drive correction!
         //         pos_fit = pos_cor;
	 //    } catch(exception& e) {std::cout<<e.what()<<"fit direction is not valid"<<std::endl;}
         // } // if valid direction

          double tfit = fitVertex.GetTime();
          if(fitVertex.ValidPosition())// && fVertex1.ValidDirection())//NOTE: check fitPosition Valid !!
          { 
            double posX=(pos_fit.X()); double posY=(pos_fit.Y());double posZ=(pos_fit.Z());
            hITR->Fill(itr);
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
		   hfitXZ->Fill(posX,posZ-108);
		   hfitRZ->Fill(pos_fit.Mag(),posZ-108);
		   hfitRhoZ->Fill(sqrt(posX*posX+posY*posY),posZ-108);
                   if(nhits>nhitCut) 
	           {		   
                     hfitXY_nhit100->Fill(posX,posY);
                     hfitXZ_nhit100->Fill(posX,posZ-108);
                     hfitRZ_nhit100->Fill(pos_fit.Mag(),posZ-108);
                     hfitRhoZ_nhit100->Fill(sqrt(posX*posX+posY*posY),posZ-108);
                   }

                   TVector3 pos_cor(pos_fit.X(),pos_fit.Y(),pos_fit.Z()-108);
                   if(itr<0.3) {
                     hfitRhoZ_ITRlt0p3->Fill(sqrt(posX*posX+posY*posY),posZ-108);
		     hfitRZ_ITRlt0p3->Fill(sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108)),posZ-108);
                   }
                   if(itr>0.3) {
                     hfitRhoZ_ITRht0p3->Fill(sqrt(posX*posX+posY*posY),posZ-108);
                     hfitRZ_ITRht0p3->Fill(sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108)),posZ-108);
                   }
                   if(abs(itr-0.3)<0.1) {
                     hfitRhoZ_ITR0p3->Fill(sqrt(posX*posX+posY*posY),posZ-108);
                     hfitRZ_ITR0p3->Fill(sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108)), posZ-108);
                   }
                   hSrcToEvt->Fill( (srcPos-pos_cor).Mag() );
		   hNhitsVsX->Fill(nhits,pos_fit.X());hNhitsVsY->Fill(nhits,pos_fit.Y());hNhitsVsZ->Fill(nhits,pos_fit.Z());
                   hNhitsVsR->Fill(nhits,pos_fit.Mag());
		   hFlatVsZ->Fill((posX*posX+posY*posY)/(6000*6000),posZ);
		   hPhiTheta->Fill(pos_fit.Phi(),pos_fit.CosTheta());
                   hLogL->Fill(LogL); hLogLvsNhit->Fill(rev.GetNhits(),LogL);
                   hLogLvsFitT->Fill(tfit,LogL);hLogLvsFitX->Fill(posX,LogL);hLogLvsFitY->Fill(posY,LogL);hLogLvsFitZ->Fill(posZ,LogL);
                   hscaleLogL->Fill(scaleLogL); hscaleLogLvsNhit->Fill(rev.GetNhits(),scaleLogL);
                   hscaleLogLvsFitT->Fill(tfit,scaleLogL);hscaleLogLvsFitX->Fill(posX,scaleLogL);hscaleLogLvsFitY->Fill(posY,scaleLogL);hscaleLogLvsFitZ->Fill(posZ,scaleLogL);
                   if( sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))<6000 && posZ-108>5100) hITR_scint->Fill(itr);
		   if( sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))>6005 || posZ-108<5100) hITR_water->Fill(itr);
		   if( sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))>6005) hITR_exwater->Fill(itr);
 		   ///Leon cuts
		   if( sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))<6000 && posZ-108>5100 && sqrt(posX*posX+posY*posY)<1000 && sqrt(posX*posX+posY*posY)>400 )
                   {
		     hLogL_leonCut->Fill(LogL); hLogLvsNhit_leonCut->Fill(rev.GetNhits(),LogL);
                     hLogLvsFitT_leonCut->Fill(tfit,LogL);hLogLvsFitX_leonCut->Fill(posX,LogL);hLogLvsFitY_leonCut->Fill(posY,LogL);hLogLvsFitZ_leonCut->Fill(posZ,LogL);
                     hscaleLogL_leonCut->Fill(scaleLogL); hscaleLogLvsNhit_leonCut->Fill(rev.GetNhits(),scaleLogL);
                     hscaleLogLvsFitT_leonCut->Fill(tfit,scaleLogL);hscaleLogLvsFitX_leonCut->Fill(posX,scaleLogL);hscaleLogLvsFitY_leonCut->Fill(posY,scaleLogL);
		     hscaleLogLvsFitZ_leonCut->Fill(posZ,scaleLogL);
		   }
		   //  std::cout<<"fitted time  "<<tfit<<std::endl;
		   Double_t radius = pos_fit.Mag();
		   //Double_t fEnergy = fitVertex.GetEnergy();
		   hNhitsFECD->Fill(rev.GetNhits());
		   hNhitVsFlat->Fill(rev.GetNhits(),(posX*posX+posY*posY)/(6000*6000));
		   //hfEnergy->Fill(fEnergy);
		   //double temp1fit = pow(u_fit*pos_fit,2)-pos_fit.Mag2()+rPSUP*rPSUP;
		   //double dPSUPfit=-u_fit*pos_fit+sqrt(temp1fit);
		   //hMfit->Fill(dPSUPfit); //distance
                   bool condition = (pos_fit.Mag()<6000 && pos_fit.Z()>4435);
                   for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++) 
                   {
			TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
			double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
			hHitTime->Fill(hitTime);

                        lightPath.CalcByPositionPartial( pos_fit, pmtpos );
                        double distInInnerAV = lightPath.GetDistInInnerAV();
                        double distInAV = lightPath.GetDistInAV();
                        double distInWater = lightPath.GetDistInWater();
                        double distInUpperTarget = lightPath.GetDistInUpperTarget();
                        double distInLowerTarget = lightPath.GetDistInLowerTarget();
                        const double transitTime = effectiveVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );

                        //lightPath.CalcByPosition( pos_fit, pmtpos );
			//double distInInnerAV = lightPath.GetDistInInnerAV();
                        //double distInAV = lightPath.GetDistInAV();
                        //double distInWater = lightPath.GetDistInWater();
                        //const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
                        //const double transitTime = effectiveVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater );
			//const double transitTime = groupVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );

			double tRes = ( hitTime - transitTime - tfit );
			htRes->Fill( tRes );
                        double cosTheta = (pmtpos-pos_cor).Unit()*(pos_cor-srcPos).Unit();
			if( (pos_cor-srcPos).Mag()>distCut) 
			{ 
		          hcosTheta->Fill( cosTheta );
			  if(!(pos_cor.Z()>5900 && pos_cor.Z()<6100)) hcosTheta_rmOdd->Fill( cosTheta );
			  if( !( (pos_cor.Z()>5900 && pos_cor.Z()<6100) || (pos_cor.Z()>5090 && pos_cor.Z()<5110)) ) hcosTheta_rmOdd1->Fill( cosTheta );
			  hcosThetaVsNhits->Fill( cosTheta, nhits);
                          if(nhits>nhitCut) hcosTheta_nhit100->Fill( cosTheta );
			  if(nhits<nhitCut) hcosTheta_nhitLess100->Fill( cosTheta );
                          hcosThetaVsTres->Fill(cosTheta,tRes);
			  if(nhits>nhitCut) hcosThetaVsTres_nhit100->Fill(cosTheta,tRes);
			  if(nhits<nhitCut) hcosThetaVsTres_nhitLess100->Fill(cosTheta,tRes);
			  if(tRes>-5 && tRes<-1) 
		          {
		            hcosTheta_prompt->Fill( cosTheta );
		            hcosThetaVsNhits_prompt->Fill( cosTheta, nhits);
                            if(nhits>nhitCut) hcosTheta_prompt_nhit100->Fill( cosTheta );
                            if(nhits<nhitCut) hcosTheta_prompt_nhitLess100->Fill( cosTheta );
			  }

			}
			if(itr<0.3) htResITRlt0p3->Fill(tRes);
			if(itr>0.3) { htResITRht0p3->Fill(tRes);}
			if(abs(itr-0.3)<0.1) { htResITR0p3->Fill(tRes);}
			double tRes1 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n1);
                        double tRes2 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n2);
                        double tRes3 = (calpmts.GetPMT(ipmt)).GetTime()-fitVertex.GetTime()-(pmtpos-pos_fit).Mag()/(c_light/n3);
			hPosXvsTime->Fill(hitTime,posX);
                        //if(trig&(2)
                        hPMTid->Fill(calpmts.GetPMT(ipmt).GetID());
			hPMTthetaPhi->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                        hPMTthetaPhi_sinu->Fill(pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta());
			if(nhits>nhitCut) hPMTthetaPhi_sinu_nhit100->Fill(pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta());
                        if(nhits>300) hPMTthetaPhi_sinu_nhit300->Fill(pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta());
			//if(tRes1<-5 && tRes1>-10)
			//
			//
			if(tRes>-5 && tRes<-1) 
			{
		          hPMTthetaPhi_prompt->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                          hPMTthetaPhi_sinu_prompt->Fill(pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta());	
                          if(nhits>100) 
			  {
			    hPMTthetaPhi_prompt_nhit100->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                            hPMTthetaPhi_sinu_prompt_nhit100->Fill(pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta());
			  }
			} 	
                        if( pmtpos.Phi() * sin(pmtpos.Theta())<-1.45 && pmtpos.Phi() * sin(pmtpos.Theta())>-1.65 && pmtpos.Theta()>1.6 && pmtpos.Theta()<1.68 )
			{
        		  hPMTthetaPhi_odd->Fill(pmtpos.Phi(),pmtpos.CosTheta());
                          hPMTid_odd->Fill(calpmts.GetPMT(ipmt).GetID());
                          hfitRZ_odd->Fill(sqrt(posX*posX+posY*posY),posZ);
                          hFitTime_odd->Fill(tfit);
                          int crate = calpmts.GetPMT(ipmt).GetCrate();
			  hCrateCardChannel[crate]->Fill(calpmts.GetPMT(ipmt).GetCard(),calpmts.GetPMT(ipmt).GetChannel());
			}
			htRes->Fill(tRes);
                        htResVsPMTz->Fill(tRes,pmtpos.Z());
			if(nhits>100) htResVsPMTz_nhit100->Fill(tRes,pmtpos.Z());
			if(condition) 
			{ htResInScint->Fill(tRes);
                         htResInScintVsPMTz->Fill(tRes,pmtpos.Z()); 
      		         if(nhits>100) { 
		            htResInScint_nhit100->Fill(tRes);
                            htResInScintVsPMTz_nhit100->Fill(tRes,pmtpos.Z());
			  }  

			  if(nhits>300) htResInScint_nhit300->Fill(tRes);
                        }
                        else{
                          if(! (pos_fit.Z()>6000 && sqrt(pos_fit.X()*pos_fit.X()+pos_fit.Y()*pos_fit.Y())<780))
			  {
                          htResInWater->Fill(tRes);htResInWaterVsPMTz->Fill(tRes,pmtpos.Z());
                          if(nhits<100) { 
			    htResInWater_nhit100->Fill(tRes);
			    htResInWaterVsPMTz_nhit100->Fill(tRes,pmtpos.Z());
			  }
                          if(nhits<300) htResInWater_nhit300->Fill(tRes);
			  }
			}
			if(nhits>100) htRes_nhit100->Fill(tRes);
                        if(nhits>300) htRes_nhit300->Fill(tRes);

                        htResVsX->Fill(tRes,pos_fit.X());
                        if(nhits>100) htResVsX_nhit100->Fill(tRes, pos_fit.X());
                        if(nhits>300) htResVsX_nhit300->Fill(tRes, pos_fit.X());

                        htResVsY->Fill(tRes,pos_fit.Y());
                        if(nhits>100) htResVsY_nhit100->Fill(tRes, pos_fit.Y());
                        if(nhits>300) htResVsY_nhit300->Fill(tRes, pos_fit.Y());

			htResVsZ->Fill(tRes,pos_fit.Z());
                        if(nhits>100) htResVsZ_nhit100->Fill(tRes, pos_fit.Z());
                        if(nhits>300) htResVsZ_nhit300->Fill(tRes, pos_fit.Z());
			htRes_1p38486->Fill(tRes1);
     			htRes_1p41->Fill(tRes2);
			htRes_1p63->Fill(tRes3);
		   }
			   
		     hTrig->Fill(0);
		     countTimeWindow++;//count events of fitValid && FECD && tFit Cut
		     hfitX_trig->Fill(posX);hfitY_trig->Fill(posY);hfitZ_trig->Fill(posZ);hfitXZ_trig->Fill(posX,posZ);
		     hfitRZ_trig->Fill(sqrt(posX*posX+posY*posY),posZ);hPosMag_trig->Fill(pos_fit.Mag());
                     //hTheta->Fill((pos_fit-sourcePos).Unit()*u_fit);
 		     hFitTime_trig->Fill(tfit);
	    }//FECD==9188
          } //FECD loop
	}//fit position valid //!!!NOTE: this is necessary 
       }//if water fitter global valid //!!!NOTE: this is necessary
       }//try catch
       catch(exception& e)
       {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
     }//if water fitter exists
    } //nhit cuts
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
  TString newfilename = "ResolTest3_"+TString(filename); 
  TFile *fp = new TFile(newfilename,"recreate");
  fp->cd();  
  hTrig->Write(); 
  hfitX->Write();hfitY->Write();hfitZ->Write();
  hfitXY->Write();hfitXZ->Write();hfitRhoZ->Write();hfitRZ->Write();hPosMag->Write();
  hfitXY_nhit100->Write();hfitXZ_nhit100->Write();hfitRhoZ_nhit100->Write();hfitRZ_nhit100->Write();

  hfitX_trig->Write();hfitY_trig->Write();hfitZ_trig->Write();hfitXZ_trig->Write();hfitRZ_trig->Write();hPosMag_trig->Write();
  hITR->Write();hITR_water->Write();hITR_scint->Write();hITR_exwater->Write();
  htResITRlt0p3->Write();htResITRht0p3->Write();htResITR0p3->Write();
  hfitRhoZ_ITRlt0p3->Write();hfitRZ_ITRlt0p3->Write();
  hfitRhoZ_ITRht0p3->Write();hfitRZ_ITRht0p3->Write();
  hfitRhoZ_ITR0p3->Write();hfitRZ_ITR0p3->Write();

  hfitRZ_odd->Write();hFitTime_odd->Write();
  hNhitsVsX->Write();hNhitsVsY->Write();hNhitsVsZ->Write();hNhitsVsR->Write();
  hSrcToEvt->Write(); 
  hPMTid->Write();hPMTid_odd->Write();hPMTthetaPhi->Write();hPMTthetaPhi_sinu->Write();hPMTthetaPhi_sinu_nhit100->Write();hPMTthetaPhi_sinu_nhit300->Write();
  hPMTthetaPhi_sinu_leonCut->Write();

  hPMTthetaPhi_prompt->Write();
  hPMTthetaPhi_sinu_prompt->Write();
  hcosTheta->Write();
  hcosTheta_rmOdd->Write();
  hcosTheta_rmOdd1->Write();

  hcosTheta_prompt->Write();
  hcosThetaVsNhits->Write();
  hcosThetaVsNhits_prompt->Write();

  hPMTthetaPhi_prompt_nhit100->Write();
  hPMTthetaPhi_sinu_prompt_nhit100->Write();
  hcosTheta_nhit100->Write();
  hcosTheta_prompt_nhit100->Write();

  hcosTheta_nhitLess100->Write();
  hcosTheta_prompt_nhitLess100->Write();
  hcosThetaVsTres->Write();hcosThetaVsTres_nhit100->Write(); hcosThetaVsTres_nhitLess100->Write();

  hNhitsFECD->Write();hPosXvsTime->Write();
  hFitTime->Write();hHitTime->Write();
  htRes->Write(); htRes_nhit100->Write(); htRes_nhit300->Write();
  htResInScint->Write(); htResInScint_nhit100->Write(); htResInScint_nhit300->Write();
  htResInScintVsPMTz->Write();htResInScintVsPMTz_nhit100->Write();
  htResInWater->Write(); htResInWater_nhit100->Write(); htResInWater_nhit300->Write();
  htResVsPMTz->Write();htResVsPMTz_nhit100->Write(); htResInWaterVsPMTz->Write();htResInWaterVsPMTz_nhit100->Write();

  htResVsX->Write();htResVsX_nhit100->Write();htResVsX_nhit300->Write();
  htResVsY->Write();htResVsY_nhit100->Write();htResVsY_nhit300->Write();
  htResVsZ->Write();htResVsZ_nhit100->Write();htResVsZ_nhit300->Write();

  hLogL->Write(); hLogLvsNhit->Write();
  hLogLvsFitT->Write();hLogLvsFitX->Write();hLogLvsFitY->Write();hLogLvsFitZ->Write();
  hscaleLogL->Write(); hscaleLogLvsNhit->Write();
  hscaleLogLvsFitT->Write();hscaleLogLvsFitX->Write();hscaleLogLvsFitY->Write();hscaleLogLvsFitZ->Write();
  hLogL_leonCut->Write(); hLogLvsNhit_leonCut->Write();
  hLogLvsFitT_leonCut->Write();hLogLvsFitX_leonCut->Write();hLogLvsFitY_leonCut->Write();hLogLvsFitZ_leonCut->Write();
  hscaleLogL_leonCut->Write(); hscaleLogLvsNhit_leonCut->Write();
  hscaleLogLvsFitT_leonCut->Write();hscaleLogLvsFitX_leonCut->Write();hscaleLogLvsFitY_leonCut->Write();
  hscaleLogLvsFitZ_leonCut->Write();

  fp->Close();
  
}
