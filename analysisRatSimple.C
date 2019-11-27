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

void analysisRatSimple()
{
  double grVelocity = 2.17554021555098529e+02 ;
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  
  ostringstream oss;
  const int eventNumber = 13;
  UInt_t interestGTID[eventNumber]={7369689,11108354,5079885,14629491,9750261,10122239,5741456,5248454,6406449,7882360,5954038,15767175,14404357};
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  in0.open("processFile.dat");
  
  // Define all of the histograms and ntuples that are to be made  
  
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5); 
  
  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events", "runNumber:eventGTID:nHit:radius_mm:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir");
  TNtuple* calPMTinfo = new TNtuple("calPMTinfo", "Information for calPMTs", "runNumber:eventGTID:calPMTtime:calPMTid:calPMTposX:calPMTposY:calPMTposZ");

  TTree *tree= new TTree("T","SolarNu");
  
  double posX, posY, posZ, posRad, time, energy, ITR, beta14, dirX, dirY, dirZ, Nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, itr, iso ;

  UInt_t day, sec, runNumber, eventGTID;
 
  //double calPMTtime;
  //int calPMTid;
  //double calPMTposX, calPMTposY, calPMTposZ;
  double nsecs ;
  
  tree->Branch("posX",&posX,"posX/D");
  tree->Branch("posY",&posY,"posY/D");
  tree->Branch("posZ",&posZ,"posZ/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("iso",&iso,"iso/D");
  tree->Branch("dirX",&dirX,"dirX/D");
  tree->Branch("dirY",&dirY,"dirY/D");
  tree->Branch("dirZ",&dirZ,"dirZ/D");
  tree->Branch("Nhits",&Nhits,"Nhits/D");
  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  char filenames[1500];
  
  Int_t i=0;
  Int_t TriggerType ;
  Int_t nentries=0;
  ULong64_t bitword ;
  
  while(in0>>filenames){
    cout<<" filenames "<<filenames<<endl ;
    // Load the RAT file
    RAT::DU::DSReader dsReader(filenames);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();
    
    // Loop over all of the events
    for(size_t iEntry = iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
   {
	posX=-10000, posY=-10000, posZ=-10000, posRad=-10000, time=99999, energy=0, ITR=-2, beta14=-10, dirX=-2.0, dirY=-2.0, dirZ=-2.0, Nhits=0, sunDirX=0.0, sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10 ;
	
    day=0, sec=0, runNumber=0, eventGTID=0;
    nsecs=0 ;	
	const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
   	// const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
	
	size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
	
	ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

	if (iEntry%5000 == 0) {
    	  cout<<iEntry<<endl;    
	}
	
    	// Looping over triggered events in each ds
	if (iEVCount) {
      for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
	  {
	      
	      // Get the event
	      const RAT::DS::EV& ev = rDS.GetEV(iEV);
	      const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
	      runNumber = run.GetRunID();
	      eventGTID = ev.GetGTID();
	      RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();

	    bool solarGTID = 0;
	    for(int i =0; i< eventNumber;i++)
	    {
		  if(eventGTID==interestGTID[i])
		  {solarGTID = 1; break;}	
		}
	    
        if(solarGTID)
		{  
		  //if bitword==dcflags.GetFlags( dcflags.GetLatestPass() ).ToString())
		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhitsCleaned();
		  //nHit->Fill(nHitOfEvent);
		  Nhits= nHitOfEvent ;
		  
		  TriggerType = ev.GetTrigType() ;
		  //cout<<" TriggerType "<<TriggerType<<endl ;
		  
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1.;
		  TVector3 pos_fit;
		  bool directionFlag = false;
		  TVector3 u_fit;
		  const string fitName = "waterFitter";
		  if(((TriggerType&2)==2) || ((TriggerType&18)==18))
		  {
		    if( ev.FitResultExists(fitName) ) {  // It needs to exist 
			RAT::DS::FitResult fResult = ev.GetFitResult(fitName);  
				
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) { // Needs to contain position
			  if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) { // Needs a valid position
			    
			       // Get fit result
			    if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			      RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
			      pos_fit = fVertex.GetPosition();
			      time = fVertex.GetTime();
			      radiusOfEvent = pos_fit.Mag();
			      posX = pos_fit.X();
			      posY = pos_fit.Y();
			      posZ = pos_fit.Z();
			      posRad = radiusOfEvent ;
			    }
			  } // end of loop if valid position
			}// end of loop if contains position
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsDirection() ) { // Needs to contain direction
			  if( ev.GetFitResult(fitName).GetVertex(0).ValidDirection() ) { // Needs a valid direction
			    if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			      RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
			      u_fit = fVertex.GetDirection();
			      dirX = u_fit.X();
			      dirY = u_fit.Y();
			      dirZ = u_fit.Z();
			      directionFlag = true;
			    }
			  }
			}
			
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) {
			  RAT::DS::FitResult fResult = ev.GetFitResult(fitName);
			  if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			    energy= fResult.GetVertex(0).GetEnergy(); }
			}
				
			if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
			//try{ ( ev.GetClassifierResult( "ITR:waterFitter" ).GetValid()== false ) ; continue; }
			//catch(RAT::DS::ClassifierResult::NoClassifierResultError& e){cout << "No 'ITR' classifier invoked in the macro." << endl; return;}
			//try { itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );}
			//catch(RAT::DS::ClassifierResult::NoClassificationError& e){ cout<< " ITR classifier has failed "<< endl ; return ;}
			if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
			itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );		    
			
			if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
			if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
			iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );
			//cout <<" itr "<<itr<<" iso "<<iso<<endl;
			
		      } // end of fit results
		      
		      // Fill the nHit histogram if the event was inside the AV
		      if (radiusOfEvent >= 0. && radiusOfEvent <= 6000.) {
				nHitAV->Fill(nHitOfEvent);
		      }
		      
		      cosThetaToSun = -2.;  
		      // Now look the direction to the Sun
		      //if (directionFlag) 
		      {
			RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
			day = eventTime.GetDays() ;
			sec = eventTime.GetSeconds() ;
			nsecs = eventTime.GetNanoSeconds() ;
			TVector3 directionSun = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
			sunDirX = directionSun.X();
			sunDirY = directionSun.Y();
			sunDirZ = directionSun.Z();
			TVector3 eventDir; eventDir.SetXYZ(dirX,dirY, dirZ) ;
			
			cosThetaToSun = eventDir*(-1.0*directionSun);
			
		      } // end of directionFlag
		      
		      // Now store the event information in an nTuple of those events with:
		      // - high nHit (>=40)
		      // - that have all of position, direction, and solar direction information
		      // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
		      
		      // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
		    candidateEventInfo->Fill(run.GetRunID(), ev.GetGTID(), nHitOfEvent, radiusOfEvent, cosThetaToSun, pos_fit.X(), pos_fit.Y(), pos_fit.Z(), u_fit.X(), u_fit.Y(), u_fit.Z());  
		  
		  for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
          {
            TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
			RAT::DS::FitResult fResult = ev.GetFitResult(fitName);
            RAT::DS::FitVertex fVertex = fResult.GetVertex(0); 
            double calPMTposX = pmtpos.X();
            double calPMTposY = pmtpos.Y();
            double calPMTposZ = pmtpos.Z();
            double hitTime =(calpmts.GetPMT(ipmt)).GetTime();   
            UInt_t calPMTid = calpmts.GetPMT(ipmt).GetID(); 
            double tRes = (calpmts.GetPMT(ipmt)).GetTime()-fVertex.GetTime()-(pmtpos-pos_fit).Mag()/grVelocity;
         
		    calPMTinfo->Fill(run.GetRunID(), ev.GetGTID(), hitTime, calPMTid ,calPMTposX,calPMTposY,calPMTposZ);
		  }

		  tree->Fill(); 
		  } //trigger cuts
		} // end of if loop for triggered word  
	    } // End loop over triggered events
	  }
	
      } // End ds loop  
  } // end of the while loop
  // Give everything an axis label, as should always be
  // Write the histograms to fileName
  TFile *file=new TFile("FilteredSolar_RatFit.root","RECREATE");
  file->cd();
  tree->Write();
  candidateEventInfo->Write();
  calPMTinfo->Write();
}
