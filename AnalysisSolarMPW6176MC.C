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
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <vector>

void AnalysisSolarMPW6176MC()
{
  //const double driveP0 = 0.995765, driveP1 = -63.826;//new param. for Jeff
 const char* filename = "WaterMP6176_WaterSolar_NueRun_r206392_s0_p0.root";
  const double c_light = 299.792458;
  double waterRI = 1.38486;
  double grVelocity = c_light/waterRI;//2.17554021555098529e+02;
  double nhitCut = 25;
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  ostringstream oss;
//  in0.open("fList.dat");
  //!!!! energy correction for 6.17.x
  //RAT::DU::ReconCorrector *eCorr = RAT::DU::ReconCorrector::Get();
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  // Define all of the histograms and ntuples that are to be made  
  //MPW no drive correction, with drive correction and ITR+beta14 cuts
  TH1D *hITR = new TH1D("hITR","trig+FECD, ITR",100,0,1.0);
  TH1D *hBeta14 = new TH1D("hBeta14","trig+FECD, beta14",200,-1.0,4.0);
  TH1D *hThetaij = new TH1D("hThetaij","trig+FECD, thetaij",200,-1.0,4.0);
  TH2D *hITR_vs_Beta14 = new TH2D("hITR_vs_Beta14","trig+FECD, itr vs beta14",100,0,1.0,200,-1.0,4.0);
  TH1D *hTrig = new TH1D("hTrig","trigger type",50,0,50);  
 
  TH1D* radius = new TH1D( "radius", "Radius of all triggered events", 100, 0., 15000.);
  TH1D* radiusCor = new TH1D( "radiusCor", "corrected radius of all triggered events", 100, 0., 15000.);
  
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5);
  TH2D* nHit_toSun = new TH2D("nHit_toSun", "nHit vs direction to the Sun of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1.);
  TH2D* nHit_toSunAV = new TH2D("nHit_toSunAV", "nHit vs direction to the Sun of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5, 100, -1., 1.); 
  TH3D* nHit_toSun_radius = new TH3D("nHit_toSun_radius", "nHit vs direction to the Sun vs radius of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1., 100., 0., 15000.);

//  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events, without drive correction", "runNumber:eventGTID:nHit:radius:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir:itr:iso");
  
  TTree *tree= new TTree("T","Simulation");
  
  double posx, posy, posz, posRad, energyOrigin, energy, time, dirx, diry, dirz, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso, thij;
  double mcposx, mcposy, mcposz, mcdirx, mcdiry, mcdirz, energyMC;
  double posxcor, posycor, poszcor, posRadCor;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t nhits, day, sec, runNumber, subRun, eventGTID, fecdID, triggerWord;
  double nsecs;
  double Gtest, Utest;

  //std::vector<double> vtRes;
  std::vector<UInt_t> vfecd;   

  //std::vector<double>* pvtRes = &vtRes;
  std::vector<UInt_t>* pvfecd = &vfecd; 
  tree->Branch("posx",&posx,"posx/D");
  tree->Branch("posy",&posy,"posy/D");
  tree->Branch("posz",&posz,"posz/D");
//!!!!! always use drive-corrected values
//  tree->Branch("posxcor",&posxcor,"posxcor/D");
//  tree->Branch("posycor",&posycor,"posycor/D");
//  tree->Branch("poszcor",&poszcor,"poszcor/D");
//  tree->Branch("posRadCor",&posRadCor,"posRadCor/D");
  
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
  //for mc
  tree->Branch("energyMC",&energyMC,"energyMC/D");

  tree->Branch("mcposx",&mcposx,"mcposx/D");
  tree->Branch("mcposy",&mcposy,"mcposy/D");
  tree->Branch("mcposz",&mcposz,"mcposz/D");

  tree->Branch("mcdirx",&mcdirx,"mcdirx/D");
  tree->Branch("mcdiry",&mcdiry,"mcdiry/D");
  tree->Branch("mcdirz",&mcdirz,"mcdirz/D");

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
  tree->Branch("prompt_cut",&prompt_cut,"prompt_cut/D");
  //char filename[1500];
  Int_t countValid = 0; 
  Int_t i=0;
  UInt_t TriggerType;
  Int_t nentries=0;
  ULong64_t bitword;
  cout<<" filename "<<filename<<endl;
   //TFile temp(filename);
   //try{ temp.IsZombie();}
   //catch(exception& e)
   //{std::cout<<e.what()<<"file is corrupted"<<std::endl;continue;}

    // Load the RAT file
    RAT::DU::DSReader dsReader(filename);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();

    cout<<numEvents<<endl; 
    // Loop over all of the events
    for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
    {
       //posx=-10000, posy=-10000, posz=-10000, posRad=-10000, time=99999, ITR=-2, beta14=-10, dirx=-2.0, diry=-2.0, dirz=-2.0, nhits=0, sunDirX=0.0,sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10, prompt_cut = -1;
       //posTheta = -10;  posPhi = -10;
       //for drive correction, only for position !!!
       //posxcor=-10000, posycor=-10000, poszcor=-10000, posRadCor=-10000;
       //day=0, sec=0, runNumber=0, eventGTID=0;
       //nsecs=0 ;	
       const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
       // const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
       
       size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
       
       //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

       if(iEntry%10000 == 0) {cout<<iEntry<<endl;}
	
    	// Looping over triggered events in each ds
	if (iEVCount) {
    	 for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
         {
	     // Get the event
	     const RAT::DS::EV& ev = rDS.GetEV(iEV);
	     const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
	     runNumber = run.GetRunID();
	     subRun = run.GetSubRunID();
	     eventGTID = ev.GetGTID();
             const RAT::DS::MC& rmc= rDS.GetMC();

             const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
             energyMC = rmcparticle.GetKineticEnergy();
	     TVector3 posmc = rmcparticle.GetPosition();
	     TVector3 dirmc = rmcparticle.GetMomentum().Unit();
             mcposx = posmc.X(); mcposy = posmc.Y(); mcposz = posmc.Z();
             mcdirx = dirmc.X(); mcdiry = dirmc.Y(); mcdirz = dirmc.Z();

	     //RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
             //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
             {	  
		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhitsCleaned();
		  nHit->Fill(nHitOfEvent);
		  nhits= nHitOfEvent;
		  TriggerType = ev.GetTrigType();
                  //pvtRes->clear();
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1., radiusOfEventCor = -1.;
		  TVector3 pos_fit, pos_cor;
		  TVector3 u_fit;
		  const string fitName = "waterFitterMP";
                  pvfecd->clear();
                  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
                  {
                    fecdID = calpmts.GetFECDPMT(ipmt).GetID();pvfecd->push_back(fecdID);
                  }
                  bool checkFitExist = ev.FitResultExists(fitName);
		  //if(((TriggerType&2)==2) || ((TriggerType&18)==18) )
		  {
		    if(checkFitExist) 
		    {  // It needs to exist
                     try {//!!! fVertex somehow broken here 
                      RAT::DS::FitVertex fVertex = ev.GetFitResult("waterFitterMP").GetVertex(0);
                      //std::cout<<"eventGTID "<<eventGTID<<" "<<iEntry<<std::endl;
                      if( !ev.GetFitResult("waterFitterMP").GetVertex(0).ContainsPosition() ) continue;
                      if(!fVertex.ValidPosition()) continue;
                      pos_fit = fVertex.GetPosition();
		      if( !ev.GetFitResult("waterFitterMP").GetVertex(0).ContainsDirection() ) continue;
                      if(!fVertex.ValidDirection()) continue;
                      u_fit = fVertex.GetDirection();
                       //NOTE: drive correction here
		      //pos_fit = driveP0*pos_fit + driveP1*u_fit;
                      time = fVertex.GetTime();
		      //!!! z - 108
                      posx = pos_fit.X();posy = pos_fit.Y();posz = pos_fit.Z();
                      radiusOfEvent = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
		      posTheta = pos_fit.CosTheta(); posPhi = pos_fit.Phi();
                      posRad = radiusOfEvent;
                      radius->Fill(radiusOfEvent);
		      //radiusOfEventCor = pos_cor.Mag();
		      //radiusCor->Fill(radiusOfEventCor);
                      dirx = u_fit.X();diry = u_fit.Y();dirz = u_fit.Z();		
                      if( !ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) continue;
		      if(!fVertex.ValidEnergy()) continue;
                      countValid++;
		      //cout<< " Fit is Valid "<<endl ;
		      energyOrigin = fVertex.GetEnergy();
                      //energy = eCorr->CorrectEnergyRSP(energyOrigin,2);
		      energy = eCorr.CorrectEnergyRSP(energyOrigin,2);
                      Gtest =  ev.GetFitResult(fitName).GetFOM("EnergyGtest");
                      Utest =  ev.GetFitResult(fitName).GetFOM("EnergyUtest");
                      //posxcor = pos_cor.X();posycor = pos_cor.Y(); poszcor = pos_cor.Z();
                      //!!!! z-108
		      //posRadCor = sqrt(pos_cor.Perp()**2+(pos_cor.Z()-108)**2);
                      //!!! apply ITR and iso cuts here, use pos_cor 
                      //!!!!!!!!!!!!!!!
		      //!!!!!NOTE the difference use pos_cor in the following calculation
                      double countPMT_ITR = 0;
                      vector<TVector3> pmtDir;
                      const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
                      double countSelectPMT = 0;
                      //TH1F *hmode = new TH1F("hmode","",800,0,800);
                      //for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                      //{
                      //   double pmtTime =(calpmts.GetPMT(ipmt)).GetTime();
                      //   hmode->Fill(pmtTime);
                      //}
                      //int binmax = hmode->GetMaximumBin();
                      //double modeTime = hmode->GetXaxis()->GetBinCenter(binmax);
                      //delete hmode;

		      //for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                      //{
                      //  TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                      //  double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                      //  if (hitTime>modeTime+modeCut_Lo && hitTime<modeTime+modeCut_Hi) //modeCut
		      //  {	 
		      //     double tRes = hitTime - time - (pmtpos-pos_fit).Mag()/grVelocity;
                      //     pvtRes->push_back(tRes);
                      //     if(tRes<straightTresCut_Hi  && tRes> straightTresCut_Lo) { //tRes cut
                      //       countSelectPMT++;// after mode and tRes cuts
                      //       pmtDir.push_back( (pmtpos -pos_fit).Unit() );
		      //       //for beta14 cuts
                      //       if(tRes<ITR_PMT_Hi && tRes>ITR_PMT_Lo) countPMT_ITR++;
	 	      //       //const TVector3 pmtDir1 = (pmtpos - pos_fit).Unit();
		      //       //do promt time cut
		      //       prompt_cut = abs((pos_fit-pmtpos).Mag() - c_light/waterRI*(hitTime-fVertex.GetTime())); 
		      //     }
		      //  }
		      //  //cout<<prompt_cut<<" "<<(pos_fit-pmtpos).Mag()<<" "<<hitTime-fVertex.GetTime()<<" "<<c_light/1.40*(hitTime-fVertex.GetTime())<<endl;
		      //}//for PMT loop
                      //struct isoVal isoresult = IsoClassifier(pmtDir);
                      //iso = isoresult.beta14;
		      //thij = isoresult.thetaij;
                      //double ITR = countPMT_ITR/countSelectPMT;
                      //itr = ITR;
                       if( !ev.ClassifierResultExists("ITR:waterFitterMP") ) continue;
                       if ( !ev.GetClassifierResult( "ITR:waterFitterMP" ).GetValid() ) continue;
                       itr = ev.GetClassifierResult( "ITR:waterFitterMP" ).GetClassification( "ITR" );

                       if( !ev.ClassifierResultExists("isotropy:waterFitterMP") ) continue;
                       if ( !ev.GetClassifierResult( "isotropy:waterFitterMP" ).GetValid() ) continue;
                       iso = ev.GetClassifierResult( "isotropy:waterFitterMP" ).GetClassification( "snobeta14" );

		       if( !ev.ClassifierResultExists("isotropy:waterFitterMP") ) continue;
                       if ( !ev.GetClassifierResult( "isotropy:waterFitterMP" ).GetValid() ) continue;
                       thij = ev.GetClassifierResult( "isotropy:waterFitterMP" ).GetClassification( "thetaij" );

		       hITR->Fill(itr);hBeta14->Fill(iso);hITR_vs_Beta14->Fill(itr,iso); hThetaij->Fill(thij);
                      // Fill the nHit histogram if the event was inside the AV
                      if (radiusOfEventCor >= 0. && radiusOfEventCor <= 6000.) {
                         nHitAV->Fill(nHitOfEvent);
                      }
                      
                      //------------------------------------  
		       cosThetaToSun = -2.;  
		      // Now look the direction to the Sun
                       RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
		       day = eventTime.GetDays() ;
		       sec = eventTime.GetSeconds() ;
		       nsecs = eventTime.GetNanoSeconds() ;
		       TVector3 directionSun = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
		       sunDirX = directionSun.X();sunDirY = directionSun.Y();sunDirZ = directionSun.Z();
		       TVector3 eventDir; eventDir.SetXYZ(dirx,diry,dirz);
		       cosThetaToSun = eventDir*(-1.0*directionSun);
		       // Fill the historgrams
		       nHit_toSun->Fill(nHitOfEvent, cosThetaToSun);
		       if(radiusOfEventCor >= 0. && radiusOfEventCor <= 6000.) {
		         nHit_toSunAV->Fill(nHitOfEvent, cosThetaToSun);   
		       } 
                      //------------------------------------
			if(radiusOfEvent >= 0.) {
			  nHit_toSun_radius->Fill(nHitOfEvent, cosThetaToSun, radiusOfEventCor);   
			}
                       triggerWord = TriggerType; // !!! Fill trigTyp 
		       // cout<<dirx<<" "<<diry<<" "<<dirz<<" "<<cosThetaToSun<<" "<<dirx*-sunDirX+diry*-sunDirY+dirz*-sunDirZ<<endl;
                       tree->Fill();
		       // Now store the event information in an nTuple of those events with:
		       // - high nHit (>=25)
		       // - that have all of position, direction, and solar direction information
		        
		       // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
		      
		      // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
		      //if (nHitOfEvent >= nhitCut && radiusOfEvent >=0.  && cosThetaToSun >= -1.) {
		      //  candidateEventInfo->Fill(run.GetRunID(), ev.GetGTID(), nHitOfEvent, radiusOfEventCor, cosThetaToSun, pos_cor.X(), pos_cor.Y(), pos_cor.Z(), u_fit.X(), u_fit.Y(), u_fit.Z(), itr, iso);
		      //}
             }//try catch
             catch(exception& e)
             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
              }//fit exist
	    } // trig cuts
		   hTrig->Fill(TriggerType);
	   } //data cleaning
	 }//loop triggered events
    }//triggered events
   }//loop ds events 
  // Give everything an axis label, as should always be
  radius->SetXTitle( "Radius [mm]" );radius->SetYTitle( "Counts" );
  radiusCor->SetXTitle( "RadiusCor [mm]" );radiusCor->SetYTitle( "Counts" );
  nHit->SetXTitle( "nHit" );nHit->SetYTitle( "Counts" );  
  nHitAV->SetXTitle( "nHit" );nHitAV->SetYTitle( "Counts" );  
  nHit_toSun->SetXTitle( "nHit" );nHit_toSun->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun->SetZTitle( "Counts" );
  nHit_toSunAV->SetXTitle( "nHit" );nHit_toSunAV->SetYTitle( "cos(#theta) to Sun" );nHit_toSunAV->SetZTitle( "Counts" );
  //nHit_toSun_radius->SetXTitle( "nHit" );nHit_toSun_radius->SetYTitle( "cos(#theta) to Sun" );nHit_toSun_radius->SetZTitle( "Radius [mm]" );
  
  cout<<countValid<<endl;
  // Write the histograms to fileName
  TString newname = "MPWsolarTest_";
  TString fileName(filename);
  TString runID(fileName(fileName.Index("r00")+5,6));
//  TString processname = newname+runID+"_test.root";
  TString processname = newname+filename;

  TFile *file=new TFile(processname,"RECREATE");
  file->cd();
  tree->Write();
  hITR->Write();hBeta14->Write();hITR_vs_Beta14->Write();hTrig->Write();hThetaij->Write();
  radius->Write();nHit->Write();nHitAV->Write();nHit_toSun->Write();nHit_toSunAV->Write();//nHit_toSun_radius->Write();
  // candidateEventInfo->Write();
}
