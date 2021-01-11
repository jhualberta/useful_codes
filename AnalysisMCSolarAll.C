// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
//#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>
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

void AnalysisMCSolarAll(char *fname)
{
  
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  
  ostringstream oss;

  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();

  //!!! RSP correctiion
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
//  in0.open("f0.dat");

  /// timeRes
  std::vector<double> *cosPMT = new std::vector<double>;
  std::vector<double> *gtRes = new std::vector<double>;
  
  // Define all of the histograms and ntuples that are to be made  
  TH1D* radius = new TH1D( "radius", "Radius of all triggered events", 100, 0., 15000.);
  
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5); 
  
  TH2D* nHit_toSun = new TH2D("nHit_toSun", "nHit vs direction to the Sun of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1.);
  TH2D* nHit_toSunAV = new TH2D("nHit_toSunAV", "nHit vs direction to the Sun of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5, 100, -1., 1.);
  
  TH3D* nHit_toSun_radius = new TH3D("nHit_toSun_radius", "nHit vs direction to the Sun vs radius of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1., 100., 0., 15000.);
  
  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events", "runNumber:eventGTID:nHit:radius_mm:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir");
  
  TTree *tree= new TTree("T","SolarNu");
  
  double posX, posY, posZ, posRad, time, energy, energyCalib, ITR, beta14, dirX, dirY, dirZ, Nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, itr, iso, thetaij;
  
  UInt_t day, sec, runNumber, eventGTID;
  double nsecs;
  
  tree->Branch("posX",&posX,"posX/D");
  tree->Branch("posY",&posY,"posY/D");
  tree->Branch("posZ",&posZ,"posZ/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("energyCalib",&energyCalib,"energyCalib/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("iso",&iso,"iso/D");
  tree->Branch("thetaij",&thetaij,"thetaij/D");
  tree->Branch("dirX",&dirX,"dirX/D");
  tree->Branch("dirY",&dirY,"dirY/D");
  tree->Branch("dirZ",&dirZ,"dirZ/D");
  tree->Branch("Nhits",&Nhits,"Nhits/D");
  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("gtRes",&gtRes);
  tree->Branch("cosPMT", &cosPMT);
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");  
  //char filenames[1500];

  Int_t i=0;
  Int_t TriggerType ;
  Int_t nentries=0;
  ULong64_t bitword ;
  
  //while(in0>>filenames)
  {
    cout<<" filenames "<<fname<<endl;
    // Load the RAT file
    RAT::DU::DSReader dsReader( fname );
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();
    
    // Loop over all of the events
    for(size_t iEntry = iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
      {
	posX=-10000, posY=-10000, posZ=-10000, posRad=-10000, time=99999, energy=0, energyCalib=0,ITR=-2, beta14=-10, dirX=-2.0, dirY=-2.0, dirZ=-2.0, Nhits=0, sunDirX=0.0, sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10 ;
	
        day=0, sec=0, runNumber=0, eventGTID=0;
        nsecs=0 ;	
	const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
   	// const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
	
	size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
	
	// ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

	if (iEntry%15000 == 0) {
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
	      //RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
              //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
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
		  TVector3 directionOfEvent;
		  const string fitName = "waterFitter";
		  //if(  ((TriggerType&2)==2) || ((TriggerType&18)==18)   )
		    {
		      if( ev.FitResultExists(fitName) ) {  // It needs to exist 
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() && ev.GetFitResult(fitName).GetVertex(0).ContainsDirection() &&  ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) { // global valid
			  if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() && ev.GetFitResult(fitName).GetVertex(0).ValidDirection() && ev.GetFitResult(fitName).GetVertex(0).ValidEnergy() ) { // global fit valid 
			    
			    RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
			    if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			      RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
			      pos_fit = fVertex.GetPosition();
			      time = fVertex.GetTime();
			      double tfit = time;
			      radiusOfEvent = pos_fit.Mag();
			      posX = pos_fit.X();
			      posY = pos_fit.Y();
			      posZ = pos_fit.Z();
			      posRad = radiusOfEvent;
			      //radius->Fill(radiusOfEvent);
                              //cout<< " Fit is Valid "<<endl ;
			      directionOfEvent = fVertex.GetDirection();
			      dirX = directionOfEvent.X();
			      dirY = directionOfEvent.Y();
			      dirZ = directionOfEvent.Z();
			      directionFlag = true;
			      energy= fResult.GetVertex(0).GetEnergy();
                              energy = eCorr.CorrectEnergyRSP(energy,2); // Correct the reconstructed energy
			      double rho = sqrt(posX*posX+posY*posY);
                              energyCalib = eCalib.CalibrateEnergyRSP(true,energy,rho,posZ);  // Calibrate the reconstructed energy (after correction)
			      //cout<<posX<<" "<<posY<<" "<<posZ<<endl;

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

			      for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                              {
                                  TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                                  double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
          			                  /// reserved for MPW calculation
                                  //double tRes = (calpmts.GetPMT(ipmt)).GetTime()-tfit-(pmtpos-pos_fit).Mag()/grVelocity;
          			  lightPath.CalcByPosition( pos_fit, pmtpos );
                                  double distInInnerAV = lightPath.GetDistInInnerAV();
                                  double distInAV = lightPath.GetDistInAV();
                                  double distInWater = lightPath.GetDistInWater();
                                  const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater );
                                  double tRes = hitTime - transitTime - tfit;
          			  gtRes->push_back(tRes);
          
          			  double cosTheta = (pmtpos - pos_fit).Unit()*-1*directionSun;
                                  cosPMT->push_back(cosTheta);
			     }
			    }
			} // global valid
			}// global contained
			
			if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
			//try{ ( ev.GetClassifierResult( "ITR:waterFitter" ).GetValid()== false ) ; continue; }
			//catch(RAT::DS::ClassifierResult::NoClassifierResultError& e){cout << "No 'ITR' classifier invoked in the macro." << endl; return;}
			//try { itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );}
			//catch(RAT::DS::ClassifierResult::NoClassificationError& e){ cout<< " ITR classifier has failed "<< endl ; return ;}
			if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
			itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );		    
			
			if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
			//try{ (ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid()== false ) ; continue;}
			//catch(RAT::DS::ClassifierResult::NoClassifierResultError& e){cout << "No 'Beta14' classifier invoked in the macro." << endl; return;}
			//try { iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );}
			//catch(RAT::DS::ClassifierResult::NoClassificationError& e){ cout<< " Beta14 classifier has failed "<< endl ; return ;}
			if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
			iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );
			thetaij = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "thetaij" );

			//cout <<" itr "<<itr<<" iso "<<iso<<endl;
			
		      } // end of fit results
		      
		      // Fill the nHit histogram if the event was inside the AV
		      if (radiusOfEvent >= 0. && radiusOfEvent <= 6000.) {
			nHitAV->Fill(nHitOfEvent);
		      }
		      
		      cosThetaToSun = -2.;  
		      // Now look the direction to the Sun
		      //if (directionFlag) 
		      
		      // Now store the event information in an nTuple of those events with:
		      // - high nHit (>=40)
		      // - that have all of position, direction, and solar direction information
		      
		      // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
		      
		      // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
		  if (nHitOfEvent >= 40 && radiusOfEvent >=0. && directionFlag && cosThetaToSun >= -1.) {
		    candidateEventInfo->Fill(run.GetRunID(), ev.GetGTID(), nHitOfEvent, radiusOfEvent, cosThetaToSun, pos_fit.X(), pos_fit.Y(), pos_fit.Z(), directionOfEvent.X(), directionOfEvent.Y(), directionOfEvent.Z());  
		  }
		  
		  tree->Fill(); 
		    } // trigger
		} // end of if loop for triggered word  
	    } // End loop over triggered events
	}
	
      } // End ds loop  
  } // end of the while loop
  // Give everything an axis label, as should always be
  radius->SetXTitle( "Radius [mm]" );
  radius->SetYTitle( "Counts" );
  nHit->SetXTitle( "nHit" );
  nHit->SetYTitle( "Counts" );  
  nHitAV->SetXTitle( "nHit" );
  nHitAV->SetYTitle( "Counts" );  
  nHit_toSun->SetXTitle( "nHit" );
  nHit_toSun->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun->SetZTitle( "Counts" );
  nHit_toSunAV->SetXTitle( "nHit" );
  nHit_toSunAV->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSunAV->SetZTitle( "Counts" );
  nHit_toSun_radius->SetXTitle( "nHit" );
  nHit_toSun_radius->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun_radius->SetZTitle( "Radius [mm]" );
  
  // Write the histograms to fileName
  TString sss(fname);
  TString sfile = "FilteredSolar_"+sss;
  TFile *file=new TFile(sfile,"RECREATE");
  file->cd();
  tree->Write();
  radius->Write();
  nHit->Write();
  nHitAV->Write();
  nHit_toSun->Write();
  nHit_toSunAV->Write();
  nHit_toSun_radius->Write();
  candidateEventInfo->Write();
  
}
