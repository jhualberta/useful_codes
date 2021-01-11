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

void AnalysisSolarMPW655(const char* filename)
{
  const double ITRval = 0.55;const double ITR_PMT_Hi = 5.0; const double ITR_PMT_Lo = -2.5;
  const double straightTresCut_Lo = -10; const double straightTresCut_Hi = 120;
  const double modeCut_Lo = -50.5; const double modeCut_Hi = 100.5;

//  const double driveP0 = 0.9868, driveP1 = -78.417;  //6.17.x
  const double driveP0 = 0.995765, driveP1 = -63.826;//new param. for Jeff

  const double c_light = 299.792458;
  double waterRI = 1.38486;
  double grVelocity = c_light/waterRI;//2.17554021555098529e+02;
  double nhitCut = 25;
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  ostringstream oss;
//  in0.open("fList.dat");

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
  
  double posx, posy, posz, posRad, time, dirx, diry, dirz, nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso, thij;
  double posxcor, posycor, poszcor, posRadCor;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t day, sec, runNumber, subRun, eventGTID, fecdID, triggerWord;
  double nsecs;
  std::vector<double> vtRes;
  std::vector<UInt_t> vfecd;   

  std::vector<double>* pvtRes = &vtRes;
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
  tree->Branch("timeRes",&vtRes);
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("dirx",&dirx,"dirx/D");
  tree->Branch("diry",&diry,"diry/D");
  tree->Branch("dirz",&dirz,"dirz/D");
  tree->Branch("nhits",&nhits,"nhits/D");

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
  tree->Branch("theta_ij",&thij,"thij/D");
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
	     subRun = run.GetSubRunID();
	     eventGTID = ev.GetGTID();
	     RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
             if( RAT::EventIsClean( ev, dcAnalysisWord ) )
             {	  
		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhitsCleaned();
		  nHit->Fill(nHitOfEvent);
		  nhits= nHitOfEvent;
		  TriggerType = ev.GetTrigType();
                  pvtRes->clear();
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1., radiusOfEventCor = -1.;
		  TVector3 pos_fit, pos_cor;
		  TVector3 u_fit;
		  const string fitName = "multipath";
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
                     RAT::DS::FitVertex fVertex = ev.GetFitResult("multipath").GetVertex(0);
                     
		     if( fVertex.ValidPosition() )
		     {
                      //std::cout<<"eventGTID "<<eventGTID<<" "<<iEntry<<std::endl;
                      pos_fit = fVertex.GetPosition(); 
          
                      RAT::DS::FitVertex fVertex1 =ev.GetFitResult("multipathdirection").GetVertex(0);
                      if(fVertex1.ValidDirection()) {
                       try{
                         u_fit = fVertex1.GetDirection();
                       }
                       catch(exception& e)
                       {std::cout<<e.what()<<" problems in dir fit"<<std::endl;}
                        //NOTE: drive correction here
		       pos_fit = driveP0*pos_fit + driveP1*u_fit;
                       countValid++;
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
                         if (hitTime>modeTime+modeCut_Lo && hitTime<modeTime+modeCut_Hi) //modeCut
			 {	 
			    double tRes = hitTime - time - (pmtpos-pos_fit).Mag()/grVelocity;
                            pvtRes->push_back(tRes);
                            if(tRes<straightTresCut_Hi  && tRes> straightTresCut_Lo) { //tRes cut
                              countSelectPMT++;// after mode and tRes cuts
                              pmtDir.push_back( (pmtpos -pos_fit).Unit() );
			      //for beta14 cuts
                              if(tRes<ITR_PMT_Hi && tRes>ITR_PMT_Lo) countPMT_ITR++;
	 		      //const TVector3 pmtDir1 = (pmtpos - pos_fit).Unit();
		              //do promt time cut
		              prompt_cut = abs((pos_fit-pmtpos).Mag() - c_light/waterRI*(hitTime-fVertex.GetTime())); 
			    }
			 }
			 //cout<<prompt_cut<<" "<<(pos_fit-pmtpos).Mag()<<" "<<hitTime-fVertex.GetTime()<<" "<<c_light/1.40*(hitTime-fVertex.GetTime())<<endl;
		       }//for PMT loop
                       struct isoVal isoresult = IsoClassifier(pmtDir);
                       iso = isoresult.beta14;
		       thij = isoresult.thetaij;
                       double ITR = countPMT_ITR/countSelectPMT;
                       itr = ITR;
                       hITR->Fill(ITR);hBeta14->Fill(isoresult.beta14);hITR_vs_Beta14->Fill(ITR,isoresult.beta14); hThetaij->Fill(isoresult.thetaij);
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
	       }//fit Valid Direction
              }//fit Valid Position
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
