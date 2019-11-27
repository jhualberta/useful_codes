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
void AnalysisKarinN16(const char* filenames)
{
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   // cannot be used in 6.17.3
   // RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   // const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   // const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();
   // RAT::DU::LightPathCalculator lightPathMC = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path

   // TVector3 srcPos(-1120.8,1041.4,6172.5);// For N16 r251748
   TVector3 srcPos(-5861,-2524,5109);//-1120.8, 1041.4, 6108.0);

   ifstream in0, in1;
   ofstream out,outputbadfiles;
   ostringstream oss;
   //ostringstream oss;
//   in0.open("files1.dat");
   TTree *tree= new TTree("gOutput","16N");

   TVector3 gFitPosition;
   TVector3 gTruePosition;
   TVector3 gTrueDirection;
   // TVector3 gSolarDirection;

   UInt_t evtNum = 0;
   UInt_t trig = 0;
   UInt_t eventGTID = 0;
   double gTrueEnergy = 0;
   double gFitEnergy = 0;
   double gFitTime = 0;
   double posX = 0, posY = 0, posZ = 0, tfit = 0;
   double gFecdTime = 0;
   double gNhitClean = 0; //use as energy
   double itr = 0, theta_ij = 0, beta14 = 0;
   double distSrcToFitPos = 0;
   // std::vector<double>* vtRes = new std::vector<double>;//timeRes saved in vector, for one event tRes is multi-values 
   std::vector<UInt_t> *gFECD;// = new std::vector<UInt_t>;//fecdID saved in vector, for one event fecd can be multi-values

   /// timeRes
   std::vector<double> *cosPMT;
   std::vector<double> *gPMTTrueTime;// = new std::vector<double>;
   std::vector<double> *gtRes;// = new std::vector<double>;

   std::vector<double> *gQhs;// = new std::vector<double>;
   std::vector<double> *gQhl;// = new std::vector<double>;

   std::vector<double> *gPMTX;// = new std::vector<double>;
   std::vector<double> *gPMTY;// = new std::vector<double>;
   std::vector<double> *gPMTZ;// = new std::vector<double>;

  double grVelocity = 2.17554021555098529e+02;

  tree->Branch("evtNum", &evtNum);
  tree->Branch("trig", &trig);
  /// TVector3 objects  
  //tree->Branch("gFitPosition",&gFitPosition);
  tree->Branch("posX",&posX);
  tree->Branch("posY",&posY);
  tree->Branch("posZ",&posZ);
  tree->Branch("tfit",&tfit);

  //tree->Branch("gTruePosition",&gTruePosition);
  //tree->Branch("gTrueDirection",&gTrueDirection);
  //tree->Branch("gSolarDirection",&gSolarDirection);
  /// double values
  tree->Branch("cosPMT", &cosPMT);
  tree->Branch("gFitTime",&gFitTime);
  tree->Branch("gFecdTime",&gFecdTime);
  tree->Branch("gTrueEnergy",&gTrueEnergy);
  tree->Branch("gFitEnergy",&gFitEnergy);
  tree->Branch("gNhitClean",&gNhitClean);
  tree->Branch("itr",&itr);
  tree->Branch("distSrcToFitPos", &distSrcToFitPos);
  /// vector containers 
  // tree->Branch("gPMTTrueTime",&gPMTTrueTime);
  tree->Branch("gtRes",&gtRes);
  tree->Branch("gQhs",&gQhs);
  tree->Branch("gQhl",&gQhl);

  tree->Branch("gFECD",&gFECD);
  tree->Branch("gPMTX",&gPMTX);
  tree->Branch("gPMTY",&gPMTY);
  tree->Branch("gPMTZ",&gPMTZ);

  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&beta14,"beta14/D");
  tree->Branch("theta_ij",&theta_ij,"theta_ij/D");

  //char filenames[1500];
  Int_t i=0;
  Int_t TriggerType;
  Int_t nentries=0;
  ULong64_t bitword;
//  while(in0>>filenames)
   {
   cout<<" filenames "<<filenames<<endl;
    // Load the RAT file
    RAT::DU::DSReader dsReader(filenames);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();
    cout<<numEvents<<endl; 
    // Loop over all of the events
    for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
    {
	const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
   	// const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event
        gtRes->clear();gFECD->clear();
        gQhs->clear();gQhl->clear();
        gPMTX->clear();gPMTY->clear();gPMTZ->clear();
        cosPMT->clear();
	size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
	//ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
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
	      //runNumber = run.GetRunID();
	      eventGTID = ev.GetGTID();
	      //RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
              //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
              {
		  //cout<< " iEntry "<<iEntry<<" GTID "<<eventGTID<<" Applied bits "<< dcflags.GetApplied( dcflags.GetLatestPass() ).ToString() <<" bits flagged "<<dcflags.GetFlags( dcflags.GetLatestPass() ).ToString()<<endl ;
		  //cout<<"PASS 1"<<endl;
		  //if bitword==dcflags.GetFlags( dcflags.GetLatestPass() ).ToString())
		  // Fill the nHit histogram, all events should have nHit
		  gNhitClean = ev.GetNhitsCleaned();
		  TriggerType = ev.GetTrigType();
		  //cout<<" TriggerType "<<TriggerType<<endl;
		  // Now get the radius, if applicable
		  TVector3 pos_fit;
		  bool directionFlag = false;
		  TVector3 u_fit;
		  const string fitName = "waterFitter";
                  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
                  {
                    int fecdID = calpmts.GetFECDPMT(ipmt).GetID();gFECD->push_back(fecdID);
                  }

		  //if(((TriggerType&2)==2))// || ((TriggerType&18)==18))
		  {
                    const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
		    if( ev.FitResultExists(fitName) ) {  // It needs to exist 
		      if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) { // Needs to contain position
                        if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) { // Needs a valid position
                           RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
                           if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
                             RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
                             pos_fit = fVertex.GetPosition();
                             tfit = fVertex.GetTime();
                             posX = pos_fit.X();
                             posY = pos_fit.Y();
                             posZ = pos_fit.Z();
			     TVector3 srcToFitPos(posX-srcPos.X(), posY-srcPos.Y(), posZ-srcPos.Z());
                             distSrcToFitPos = srcToFitPos.Mag();
			     for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                              {
                                 TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                                 double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                                 double tRes = (calpmts.GetPMT(ipmt)).GetTime()-tfit-(pmtpos-pos_fit).Mag()/grVelocity;
				 gtRes->push_back(tRes);
                                 gPMTX->push_back(pmtpos.X());gPMTY->push_back(pmtpos.Y());gPMTZ->push_back(pmtpos.Z());
                                 double qhs = calpmts.GetPMT(ipmt).GetQHS();
                                 double qhl = calpmts.GetPMT(ipmt).GetQHL();
                                 gQhs->push_back(qhs);
                                 gQhl->push_back(qhl);
                                 double cosTheta = (pmtpos - pos_fit).Unit()*srcToFitPos.Unit();
                                 cosPMT->push_back(cosTheta);
                              }
			   }
                         } // end of loop if valid position
                       }// end of loop if contains position
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) {
			  RAT::DS::FitResult fResult = ev.GetFitResult(fitName);
			  if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			    gFitEnergy = fResult.GetVertex(0).GetEnergy(); }
			}

			if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
			if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
			itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );
			
			if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
			if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
			beta14 = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );
                        theta_ij = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "thetaij" );
			//cout <<" itr "<<itr<<" iso "<<iso<<endl;
			
		        } // end of fit results
		      
                        tree->Fill(); 
                        // gtRes->clear();gFECD->clear();
                        // gQhs->clear();gQhl->clear();
                        // gPMTX->clear();gPMTY->clear();gPMTZ->clear();
		    } //trigger cuts
		} // end of if loop for triggered word  
	    } // End loop over triggered events
	} // if triggered 
      } // End ds loop  
  } // end of the while loop
  // Give everything an axis label, as should always be
  
  // Write the histograms to fileName
  TString oldfilename(filenames);
  TString newfilename = "Tree_"+oldfilename;
  TFile *file=new TFile(newfilename,"RECREATE");
  file->cd();
  tree->Write();
  file->Close(); 
}
