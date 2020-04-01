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
//HERE FOR MC, NOT DATA!!
struct isoVal
{
   double beta14;
   double thetaij;
};//store ISO classifier results

struct isoVal IsoClassifier(vector<TVector3>);//customized ISO classifier

void searchSolarJeff1a()
{
  const double ITRval = 0.55;const double ITR_PMT_Hi = 5.0; const double ITR_PMT_Lo = -2.5;
  const double driveP0 = 0.995765, driveP1 = -63.82600;//new param. for Jeff

  const double c_light = 299.792458;
  double waterRI = 1.38486;
  double grVelocity = c_light/waterRI;//2.17554021555098529e+02;
  double nhitCut = 25;
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  ostringstream oss;
  in0.open("test1.dat");

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

  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events, without drive correction", "runNumber:eventGTID:nHit:radius:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir:itr:iso");
  
  TTree *tree= new TTree("T","Simulation");
  
  double posX, posY, posZ, posRad, time, dirX, dirY, dirZ, nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso, thij;
  double posXcor, posYcor, posZcor, posRadCor;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t day, sec, runNumber, eventGTID, fecdID;
  double nsecs;
  std::vector<double> vtRes;
  std::vector<UInt_t> vfecd;   

  std::vector<double>* pvtRes = &vtRes;
  std::vector<UInt_t>* pvfecd = &vfecd; 
  tree->Branch("posX",&posX,"posX/D");
  tree->Branch("posY",&posY,"posY/D");
  tree->Branch("posZ",&posZ,"posZ/D");

  tree->Branch("posXcor",&posXcor,"posXcor/D");
  tree->Branch("posYcor",&posYcor,"posYcor/D");
  tree->Branch("posZcor",&posZcor,"posZcor/D");
  tree->Branch("posRadCor",&posRadCor,"posRadCor/D");
  
  tree->Branch("posTheta",&posTheta,"posTheta/D");
  tree->Branch("posPhi",&posPhi,"posPhi/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("timeRes",&vtRes);
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("dirX",&dirX,"dirX/D");
  tree->Branch("dirY",&dirY,"dirY/D");
  tree->Branch("dirZ",&dirZ,"dirZ/D");
  tree->Branch("nhits",&nhits,"nhits/D");

  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");

  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&iso,"iso/D");
  tree->Branch("theta_ij",&thij,"thij/D");
  tree->Branch("prompt_cut",&prompt_cut,"prompt_cut/D");
  char filenames[1500];
  Int_t countValid = 0; 
  Int_t i=0;
  Int_t TriggerType;
  Int_t nentries=0;
  ULong64_t bitword;
  while(in0>>filenames){
   cout<<" filenames "<<filenames<<endl;
   //TFile temp(filenames);
   //try{ temp.IsZombie();}
   //catch(exception& e)
   //{std::cout<<e.what()<<"file is corrupted"<<std::endl;continue;}

    // Load the RAT file
    RAT::DU::DSReader dsReader(filenames);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();

    cout<<numEvents<<endl; 
    // Loop over all of the events
    for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
    {
       //posX=-10000, posY=-10000, posZ=-10000, posRad=-10000, time=99999, ITR=-2, beta14=-10, dirX=-2.0, dirY=-2.0, dirZ=-2.0, nhits=0, sunDirX=0.0,sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10, prompt_cut = -1;
       //posTheta = -10;  posPhi = -10;
       //for drive correction, only for position !!!
       //posXcor=-10000, posYcor=-10000, posZcor=-10000, posRadCor=-10000;
       //day=0, sec=0, runNumber=0, eventGTID=0;
       //nsecs=0 ;	
       const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
       // const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
       
       size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
       
       ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

       if(iEntry%10000 == 0) {cout<<iEntry<<endl;}
	
    	// Looping over triggered events in each ds
	if (iEVCount) {
    	 for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
         {
	     // Get the event
	     const RAT::DS::EV& ev = rDS.GetEV(iEV);
	     const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
	     runNumber = run.GetRunID();
	     eventGTID = ev.GetGTID();
	     int fPass = ev.GetDataCleaningFlags().GetLatestPass();
	     int dcFlagged = ev.GetDataCleaningFlags().GetFlags(fPass).GetULong64_t(0);
             int p1_mask = 0xFB0000017FFE;
             // p2_mask = 0xfb3ffff97ffe;
             if ( (dcFlagged & p1_mask ) != p1_mask )  continue;
             else
	     //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
             {	  
		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhitsCleaned();
		  nHit->Fill(nHitOfEvent);
		  nhits= nHitOfEvent;
		  TriggerType = ev.GetTrigType() ;
                  pvtRes->clear();
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1., radiusOfEventCor = -1.;
		  TVector3 pos_fit, pos_cor;
		  TVector3 u_fit;
		  const string fitName = "alberta";
                  pvfecd->clear();
                  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
                  {
                    fecdID = calpmts.GetFECDPMT(ipmt).GetID();pvfecd->push_back(fecdID);
                  }
                  bool checkFitExist = ev.FitResultExists(fitName);
		  if(((TriggerType&2)==2) || ((TriggerType&18)==18) )
		  {
		    if(checkFitExist) 
		    {  // It needs to exist
                     try{//!!! fVertex somehow broken here 
                     RAT::DS::FitVertex fVertex = ev.GetFitResult("alberta").GetVertex(0);
                     
		     if(fVertex.ValidPosition()) 
		     {  
                      //std::cout<<"eventGTID "<<eventGTID<<" "<<iEntry<<std::endl;
                      pos_fit = fVertex.GetPosition(); 
          
                      RAT::DS::FitVertex fVertex1 =ev.GetFitResult("albertadirection").GetVertex(0);
                      if(fVertex1.ValidDirection()) {
                       try{
                         u_fit = fVertex1.GetDirection();
                       }
                       catch(exception& e)
                       {std::cout<<e.what()<<" problems in dir fit"<<std::endl;}
                        //NOTE: drive correction here
		       pos_cor = driveP0*pos_fit + driveP1*u_fit;
                       countValid++;
                       time = fVertex.GetTime();
                       radiusOfEvent = pos_fit.Mag();
                       posX = pos_fit.X();posY = pos_fit.Y();posZ = pos_fit.Z();
                       posTheta = pos_fit.CosTheta(); posPhi = pos_fit.Phi();
                       posRad = radiusOfEvent;
                       radius->Fill(radiusOfEvent);
		       radiusOfEventCor = pos_cor.Mag();
		       radiusCor->Fill(radiusOfEventCor);
                       dirX = u_fit.X();dirY = u_fit.Y();dirZ = u_fit.Z();		

                       posXcor = pos_cor.X();posYcor = pos_cor.Y(); posZcor = pos_cor.Z();
		       posRadCor = pos_cor.Mag();
                       //!!! apply ITR and iso cuts here, use pos_cor 
                       //!!!!!!!!!!!!!!!
		       //!!!!!NOTE the difference use pos_cor in the following calculation
                       double countPMT_ITR = 0;
                       vector<TVector3> pmtDir;
                       const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
                       for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                       {
                         TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                         pmtDir.push_back( (pmtpos -pos_cor).Unit() );
                         double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                         double tRes = (calpmts.GetPMT(ipmt)).GetTime()-fVertex.GetTime()-(pmtpos-pos_cor).Mag()/grVelocity;
                         pvtRes->push_back(tRes);
                         if(tRes<ITR_PMT_Hi && tRes>ITR_PMT_Lo) countPMT_ITR++;
                         //for beta14 cuts
                         const TVector3 pmtDir1 = (pmtpos - pos_cor).Unit();
                           
		         //do promt time cut
		         prompt_cut = abs((pos_cor-pmtpos).Mag() - c_light/waterRI*(hitTime-fVertex.GetTime())); 
		         //cout<<prompt_cut<<" "<<(pos_fit-pmtpos).Mag()<<" "<<hitTime-fVertex.GetTime()<<" "<<c_light/1.40*(hitTime-fVertex.GetTime())<<endl;
		       }//for PMT loop
                       struct isoVal isoresult = IsoClassifier(pmtDir);
                       iso = isoresult.beta14;
		       thij = isoresult.thetaij;
                       double ITR = countPMT_ITR/double(calpmts.GetCount()); 
                       itr = ITR;	
                       hITR->Fill(ITR);hBeta14->Fill(isoresult.beta14);hITR_vs_Beta14->Fill(ITR,isoresult.beta14); hThetaij->Fill(isoresult.thetaij);
 
                       tree->Fill();
                       }//fit Valid Direction
                      }//fit Valid Position
       
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
			TVector3 eventDir; eventDir.SetXYZ(dirX,dirY,dirZ);
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
		     // Now store the event information in an nTuple of those events with:
		     // - high nHit (>=25)
		     // - that have all of position, direction, and solar direction information
		      
		     // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
		      
		  // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
		  if (nHitOfEvent >= nhitCut && radiusOfEvent >=0.  && cosThetaToSun >= -1.) {
		    candidateEventInfo->Fill(run.GetRunID(), ev.GetGTID(), nHitOfEvent, radiusOfEventCor, cosThetaToSun, pos_cor.X(), pos_cor.Y(), pos_cor.Z(), u_fit.X(), u_fit.Y(), u_fit.Z(), itr, iso);
		  }
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
  }// end of the while loop for reading files
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
  TString newname = "JeffSolar_";
  TString fileName(filenames);
  TString runID(fileName(fileName.Index("r00")+5,6));
  TString processname = newname+runID+"_solar.root";
  TFile *file=new TFile(processname,"RECREATE");
  file->cd();
  tree->Write();
  hITR->Write();hBeta14->Write();hITR_vs_Beta14->Write();hTrig->Write();hThetaij->Write();
  radius->Write();nHit->Write();nHitAV->Write();nHit_toSun->Write();nHit_toSunAV->Write();//nHit_toSun_radius->Write();
  candidateEventInfo->Write();
}

struct isoVal IsoClassifier(vector<TVector3> pmtDir)
{
    double thetaij = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0, p4 = 0.0, tij = 0.0;
    struct isoVal result;
    result.beta14 = 0;
    result.thetaij = 0;
    if(pmtDir.size()==0) return result;
    
    for( size_t iPMT = 0; iPMT < pmtDir.size(); ++iPMT )
    {
      const TVector3 pmtDir1 = pmtDir[iPMT];
      for( size_t iPMT2 = iPMT + 1; iPMT2 < pmtDir.size(); ++iPMT2 )
        {
          const TVector3 pmtDir2 = pmtDir[iPMT2];
          thetaij = pmtDir[iPMT].Angle(pmtDir[iPMT2]);
          tij += thetaij;
          const double cosTheta12 = pmtDir1.Dot( pmtDir2 );
          p1 += cosTheta12; // px is Legendre Polynomial order x
          const double cosThetaSquared = cosTheta12 * cosTheta12;
          p2 += ( 3.0 * cosThetaSquared - 1.0 ) / 2.0;
          p3 += ( 5.0 * cosThetaSquared - 3.0 ) * cosTheta12 / 2.0;
          p4 += ( 35.0 * cosThetaSquared * cosThetaSquared - 30.0 * cosThetaSquared + 3.0 ) / 8.0;
        }
    }
    p1 = 2.0 * p1 / static_cast<double>( pmtDir.size() * ( pmtDir.size() - 1 ) );
    p2 = 2.0 * p2 / static_cast<double>( pmtDir.size() * ( pmtDir.size() - 1 ) );
    p3 = 2.0 * p3 / static_cast<double>( pmtDir.size() * ( pmtDir.size() - 1 ) );
    p4 = 2.0 * p4 / static_cast<double>( pmtDir.size() * ( pmtDir.size() - 1 ) );
    tij = 2.0 * tij / static_cast<double>( pmtDir.size() * ( pmtDir.size() - 1 ) );
    const double beta14_temp = p1 + 4 * p4;
    result.beta14 = beta14_temp;
    result.thetaij = tij;
    return result;
}
