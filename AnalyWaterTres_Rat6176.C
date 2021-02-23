// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TVector3.h>
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <vector>
//Modified from AnalysisKarinN16.C, get tRes only!
using namespace std;
void AnalyWaterTres_Rat6176()
{
  const char* filename = "Analysis_r0000200658_s000_p004.root";
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  // cannot be used in 6.17.3
  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();

//  RAT::DU::LightPathCalculator lightPathMC = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  // TVector3 srcPos(-1120.8,1041.4,6172.5);// For N16 r251748
  ifstream in0, in1;
  ofstream out,outputbadfiles;
 //!!! RSP correctiion
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();

  //ostringstream oss;
  TTree *tree= new TTree("T","solar");
  double tfit=0; 
  /// timeRes
  std::vector<double> *gCosPMT = new std::vector<double>;
  std::vector<double> *gtRes = new std::vector<double>;
  std::vector<double> *gCosPMTfit = new std::vector<double>;

  std::vector<UInt_t> vfecd;

//  std::vector<double> *gtResMC = new std::vector<double>;

  double grVelocity = 2.17554021555098529e+02;
  double posx, posy, posz, posRad, energyOrigin, energy, time, dirx, diry, dirz, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso, thij;
  double posxcor, posycor, poszcor, posRadCor;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t nhits, day, sec, runNumber, subRun, eventGTID, fecdID, triggerWord;
  double nsecs;
  double Gtest, Utest, posFOM, posFOM2, scaleLogL, dirFOM, dirFOM2, dirscaleLogL;

  /// TVector3 objects  
  //tree->Branch("gFitPosition",&gFitPosition);
  //tree->Branch("gTruePosition",&gTruePosition);
  //tree->Branch("gTrueDirection",&gTrueDirection);
  //tree->Branch("gSolarDirection",&gSolarDirection);
  /// double values
  tree->Branch("gCosPMT", &gCosPMT);
  tree->Branch("gCosPMTfit", &gCosPMTfit);
  /// vector containers 
  // tree->Branch("gPMTTrueTime",&gPMTTrueTime);
  tree->Branch("gtRes",&gtRes);
//  tree->Branch("cosPMTMC", &cosPMTMC);
//  tree->Branch("gtResMC",&gtResMC);

  std::vector<UInt_t>* pvfecd = &vfecd; 
  tree->Branch("posx",&posx,"posx/D");
  tree->Branch("posy",&posy,"posy/D");
  tree->Branch("posz",&posz,"posz/D");
  tree->Branch("posTheta",&posTheta,"posTheta/D");
  tree->Branch("posPhi",&posPhi,"posPhi/D");
  tree->Branch("time",&time,"time/D");
  //tree->Branch("timeRes",&vtRes);
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("dirx",&dirx,"dirx/D");
  tree->Branch("diry",&diry,"diry/D");
  tree->Branch("dirz",&dirz,"dirz/D");
  tree->Branch("nhits",&nhits,"nhits/i");
  tree->Branch("energyOrigin",&energyOrigin,"energyOrigin/D");
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("subRun", &subRun, "subRun/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("triggerWord",&triggerWord,"triggerWord/I");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&iso,"iso/D");
  tree->Branch("thetaij",&thij,"thij/D");
  tree->Branch("Gtest",&Gtest,"Gtest/D");
  tree->Branch("Utest",&Utest,"Utest/D");
  tree->Branch("posFOM",&posFOM,"posFOM/D");
  tree->Branch("posFOM2",&posFOM2,"posFOM2/D");
  tree->Branch("scaleLogL",&scaleLogL,"scaleLogL/D");
  tree->Branch("dirFOM",&dirFOM,"dirFOM/D");
  tree->Branch("dirFOM2",&dirFOM2,"dirFOM2/D");
  tree->Branch("dirscaleLogL",&dirscaleLogL,"dirscaleLogL/D"); 
  tree->Branch("prompt_cut",&prompt_cut,"prompt_cut/D");

  //char filename[1500];

  Int_t i=0;
  Int_t TriggerType=0;
  Int_t nentries=0;
  ULong64_t p2_mask = 0xfb7ffff97ffe;//0xFB0000017FFE for 6.5.3, 0xfb7ffff97ffe for 6.17.6

  cout<<"!!!!! filename "<<filename<<endl;
   // Load the RAT file
  RAT::DU::DSReader dsReader(filename);
  // To later get the run info
  const RAT::DS::Run& run = dsReader.GetRun();
  size_t numEvents = dsReader.GetEntryCount();
  cout<<numEvents<<endl; 
  // Loop over all of the events
  for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
  {
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
    const RAT::DS::Run& run = dsReader.GetRun();
    gtRes->clear();
    gCosPMT->clear();
    gCosPMTfit->clear();
    //MC only
    //const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event
    //gtResMC->clear();
    //cosPMTMC->clear();

    size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
    ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
    if (iEntry%10000 == 0) {
      cout<<"processed "<<iEntry<<endl;
    }
    
    // Looping over triggered events in each ds
    if (iEVCount) {
     for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
     {
          // Get the event
          const RAT::DS::EV& ev = rDS.GetEV(iEV);
          const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
          runNumber = run.GetRunID();
          RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
          ULong64_t dcApplied;
          ULong64_t dcFlagged;
	  Int_t fPass = ev.GetDataCleaningFlags().GetLatestPass();
          if(fPass >= 0)
          {
            dcApplied = ev.GetDataCleaningFlags().GetApplied(fPass).GetULong64_t(0);
            dcFlagged = ev.GetDataCleaningFlags().GetFlags(fPass).GetULong64_t(0);
          }
          else continue;

	  //if(RAT::EventIsClean( ev, dcAnalysisWord )
          if( ((dcApplied & p2_mask) & dcFlagged) == (dcApplied & p2_mask) )
	  {
	     //cout<<" iEntry "<<iEntry<< endl;
	     TVector3 pos_fit;
	     bool directionFlag = false;
	     TVector3 u_fit;
	     const string fitName = "waterFitter";
             pvfecd->clear();
             for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
             {
               fecdID = calpmts.GetFECDPMT(ipmt).GetID();pvfecd->push_back(fecdID);
             }
//             // MC
//	     //TVector3 pos_mc =mc.GetMCParticle(0).GetPosition();
//	     //TVector3 u_mc =  mc.GetMCParticle(0).GetMomentum().Unit();
             bool checkFitExist = ev.FitResultExists(fitName);
	     nhits= ev.GetNhitsCleaned();
             if(nhits>15) // nhit cut for cosTheta calculation!!
             {
	           //if(((TriggerType&2)==2))// || ((TriggerType&18)==18))
	           {
                const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
                  if(checkFitExist)
                  {  // It needs to exist
                    try {//!!! fVertex somehow broken here 
                     RAT::DS::FitVertex fVertex = ev.GetFitResult("waterFitter").GetVertex(0);
		     if( fVertex.ValidPosition() )
		     {
                      eventGTID = ev.GetGTID();
	       	      pos_fit = fVertex.GetPosition(); 
          
                      RAT::DS::FitVertex fVertex1 =ev.GetFitResult("waterFitter").GetVertex(0);
                      if(fVertex1.ValidDirection()) {
                       try{
                         u_fit = fVertex1.GetDirection();
                       }
                       catch(exception& e)
                       {std::cout<<e.what()<<" problems in dir fit"<<std::endl;}
                        //NOTE: drive correction here
		       //pos_fit = driveP0*pos_fit + driveP1*u_fit;
		      time = fVertex.GetTime();
		      //!!! z - 108
                      posx = pos_fit.X();posy = pos_fit.Y();posz = pos_fit.Z();
                      posFOM = ev.GetFitResult(fitName).GetFOM("PositionLogL");
                      posFOM2 = ev.GetFitResult(fitName).GetFOM("PositionSelectedNHit");
		      scaleLogL = posFOM/posFOM2;
                      dirFOM = ev.GetFitResult(fitName).GetFOM("DirectionLogL");
                      dirFOM2 = ev.GetFitResult(fitName).GetFOM("DirectionSelectedNHit");
                      dirscaleLogL = dirFOM/dirFOM2;
                      double radiusOfEvent = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
                      posTheta = pos_fit.CosTheta(); posPhi = pos_fit.Phi();
                      posRad = radiusOfEvent;
		      //radiusOfEventCor = pos_cor.Mag();
		      //radiusCor->Fill(radiusOfEventCor);
                      dirx = u_fit.X();diry = u_fit.Y();dirz = u_fit.Z();		
                      if( !ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) continue;
		      RAT::DS::FitResult fResult = ev.GetFitResult(fitName);
		      if(!fResult.GetValid()) continue;
		      //cout<< " Fit is Valid "<<endl ;
                      energyOrigin = fResult.GetVertex(0).GetEnergy();
                      //energy = eCorr->CorrectEnergyRSP(energyOrigin,2);
		      energy = eCorr.CorrectEnergyRSP(energyOrigin,2);
                      Gtest =  fResult.GetFOM("EnergyGtest");
                      Utest =  fResult.GetFOM("EnergyUtest");

		      if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
                      if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
                      itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );

                      if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
                      if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
                      iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );

                      if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
                      if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
                      thij = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "thetaij" );

		      // Now look the direction to the Sun
                      cosThetaToSun = -2.;
                      RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
                      day = eventTime.GetDays() ;
                      sec = eventTime.GetSeconds() ;
                      nsecs = eventTime.GetNanoSeconds() ;
                      TVector3 directionSun = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
                      sunDirX = directionSun.X();sunDirY = directionSun.Y();sunDirZ = directionSun.Z();
                      TVector3 eventDir; eventDir.SetXYZ(dirx,diry,dirz);
                      cosThetaToSun = eventDir*(-1.0*directionSun);
		      for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                      {
                         TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                         double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                         /// reserved for MPW calculation
                         lightPath.CalcByPosition( pos_fit, pmtpos );
                         double distInInnerAV = lightPath.GetDistInInnerAV();
                         double distInAV = lightPath.GetDistInAV();
                         double distInWater = lightPath.GetDistInWater();
                         const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater );
                         //double tRes = (calpmts.GetPMT(ipmt)).GetTime()-time-(pmtpos-pos_fit).Mag()/grVelocity;
                         double tRes = hitTime - transitTime - time;
                         gtRes->push_back(tRes);

                         double cosThetaTrue = (pmtpos - pos_fit).Unit()*-1*directionSun;
                         gCosPMT->push_back(cosThetaTrue);
                         double cosThetaFit = (pmtpos - pos_fit).Unit()*u_fit;
                         gCosPMTfit->push_back(cosThetaFit);
                         /// for MC
                         //lightPathMC.CalcByPosition( pos_mc, pmtpos );
                         //double distInInnerAV2 = lightPathMC.GetDistInInnerAV();
                         //double distInAV2 = lightPathMC.GetDistInAV();
                         //double distInWater2 = lightPathMC.GetDistInWater();
                         //const double transitTime2 = groupVelocity.CalcByDistance( distInInnerAV2, distInAV2, distInWater2 ); // Assumes a 400nm photon
                         //  double tResMC = hitTime - transitTime2 - 390 + rDS.GetMCEV(iEV).GetGTTime(); //has problem
                         //gtResMC->push_back(tResMC);
                         //double cosThetaMC = (pmtpos - pos_mc).Unit()*-1*directionSun;
                         //cosPMTMC->push_back(cosThetaMC);
		      } //loop PMT
		      triggerWord = TriggerType; // !!! Fill trigTyp
                      tree->Fill();
                  }//fit Valid Direction
                }//fit Valid Position
             }//try catch
             catch(exception& e)
             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
            } //fitName exist
	    } //triggered
	   }//nhit cut>15
	}//data clean 
     }// event loop
    }// triggered event>1
   } // End ds event loop  
  // Give everything an axis label, as should always be

  delete gCosPMT;
  delete gCosPMTfit;
  delete gtRes;

//  delete cosPMTMC;
//  delete gtResMC;

//  cout<<"OK4?"<<endl;

  // Write the histograms to fileName
  TString oldfilename(filename);
//  cout<<"OK3? "<<oldfilename<<endl;

  TString newfilename = "ExtractTres_"+oldfilename;
  cout<<newfilename<<endl;
  TFile *file=new TFile(newfilename,"RECREATE");
  file->cd();
  tree->Write();
  file->Close(); 
}
