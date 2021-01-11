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
void AnalyWaterTres()
{
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  // cannot be used in 6.17.3
  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();

  RAT::DU::LightPathCalculator lightPathMC = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  // TVector3 srcPos(-1120.8,1041.4,6172.5);// For N16 r251748
  ifstream in0, in1;
  ofstream out,outputbadfiles;
  //ostringstream oss;
  TTree *tree= new TTree("T","solar");
  double tfit=0; 
  /// timeRes
  std::vector<double> *cosPMT = new std::vector<double>;
  std::vector<double> *gtRes = new std::vector<double>;

  std::vector<double> *cosPMTMC = new std::vector<double>;
  std::vector<double> *gtResMC = new std::vector<double>;

  UInt_t eventGTID = 0;
  UInt_t runNumber = 0;
  double grVelocity = 2.17554021555098529e+02;

  /// TVector3 objects  
  //tree->Branch("gFitPosition",&gFitPosition);
  //tree->Branch("gTruePosition",&gTruePosition);
  //tree->Branch("gTrueDirection",&gTrueDirection);
  //tree->Branch("gSolarDirection",&gSolarDirection);
  /// double values
  tree->Branch("cosPMT", &cosPMT);
  /// vector containers 
  // tree->Branch("gPMTTrueTime",&gPMTTrueTime);
  tree->Branch("gtRes",&gtRes);

  tree->Branch("cosPMTMC", &cosPMTMC);
  tree->Branch("gtResMC",&gtResMC);
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("runNumber",&runNumber,"runNumber/i");

  //char filenames[1500];

  Int_t i=0;
  Int_t TriggerType;
  Int_t nentries=0;
  ULong64_t bitword;

  const char* filenames = "WaterSolar_NueRun_r206392_s0_p0.root";

  cout<<"!!!!! filenames "<<filenames<<endl;
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
    const RAT::DS::Run& run = dsReader.GetRun();

    gtRes->clear();
    cosPMT->clear();

    //MC only
    const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event
    gtResMC->clear();
    cosPMTMC->clear();

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
          //RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
          //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
          {
	     //cout<<" iEntry "<<iEntry<< endl;
	     TVector3 pos_fit;
	     bool directionFlag = false;
	     TVector3 u_fit;
	     const string fitName = "waterFitter";
             // MC
	     TVector3 pos_mc =mc.GetMCParticle(0).GetPosition();
	     TVector3 u_mc =  mc.GetMCParticle(0).GetMomentum().Unit();
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
                      //solar information
                      RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
                      UInt_t day = eventTime.GetDays();
                      UInt_t sec = eventTime.GetSeconds();
                      UInt_t nsecs = eventTime.GetNanoSeconds();
                      TVector3 directionSun = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
                      runNumber = run.GetRunID();
                      eventGTID = ev.GetGTID();

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

                        //lightPathMC.CalcByPosition( pos_mc, pmtpos );
                        //double distInInnerAV2 = lightPathMC.GetDistInInnerAV();
                        //double distInAV2 = lightPathMC.GetDistInAV();
                        //double distInWater2 = lightPathMC.GetDistInWater();
                        //const double transitTime2 = groupVelocity.CalcByDistance( distInInnerAV2, distInAV2, distInWater2 ); // Assumes a 400nm photon
                        //  double tResMC = hitTime - transitTime2 - 390 + rDS.GetMCEV(iEV).GetGTTime(); //has problem
                        //gtResMC->push_back(tResMC);

                        double cosThetaMC = (pmtpos - pos_mc).Unit()*-1*directionSun;
                        cosPMTMC->push_back(cosThetaMC);
			}
		     } //fit valid
                 } // end of loop if valid position
               }// end of loop if contains position
			
	} // end of fit results
	//cout<<"OK?"<<endl;
        tree->Fill();
// gtRes->clear();gFECD->clear();
// gQhs->clear();gQhl->clear();
// gPMTX->clear();gPMTY->clear();gPMTZ->clear();
		    } //trigger cuts
		   } // end of if loop for triggered word  
	    } // End loop over triggered events
	   } // if triggered 
   } // End ds loop  
  // Give everything an axis label, as should always be

  delete cosPMT;
  delete gtRes;

  delete cosPMTMC;
  delete gtResMC;

  cout<<"OK4?"<<endl;

  // Write the histograms to fileName
  TString oldfilename(filenames);
  cout<<"OK3? "<<oldfilename<<endl;

  TString newfilename = "ExtractTres_"+oldfilename;
  cout<<newfilename<<endl;
  TFile *file=new TFile(newfilename,"RECREATE");
  file->cd();
  tree->Write();
  file->Close(); 
}
