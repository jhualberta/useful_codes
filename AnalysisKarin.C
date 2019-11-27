// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)
// !!!!!!!!!!! This is for data !!!!!!!!!!!!!!!!!!!!!!!!!
// Jie Hu Oct 10, 2019
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
using namespace::std;
void AnalysisKarin()
{

   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   //const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();

   RAT::DU::LightPathCalculator lightPathMC = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path

   // TVector3 srcPos(-1120.8,1041.4,6172.5);// For N16 r251748
   TVector3 srcPos(-1120.8, 1041.4, 6108.0);

   ifstream in0, in1;
   //ostringstream oss;
   in0.open("test1.dat");
   TTree *tree= new TTree("gOutput","16N");
  
   TVector3 gFitPosition;
   TVector3 gTruePosition;
   TVector3 gTrueDirection;
   TVector3 gSolarDirection;
 
   UInt_t evtNum = 0;
   UInt_t trig = 0;
   double gTrueEnergy = 0;
   double gFitEnergy = 0;
   double gFitTime = 0;
   
   double gFecdTime = 0;
   double gTimeResidual = 0;
   double gNhitClean = 0; //use as energy
   double itr = 0;
   double distSrcToFitPos = 0;  
   // std::vector<double>* vtRes = new std::vector<double>;//timeRes saved in vector, for one event tRes is multi-values 
   std::vector<UInt_t> *gFECD;// = new std::vector<UInt_t>;//fecdID saved in vector, for one event fecd can be multi-values
  
   /// timeRes
   double cosTheta = 0;
   std::vector<double> *gPMTTrueTime;// = new std::vector<double>;
   std::vector<double> *gPMTRecoTime;// = new std::vector<double>;
  
   std::vector<double> *gQhs;// = new std::vector<double>;
   std::vector<double> *gQhl;// = new std::vector<double>;
  
   std::vector<double> *gPMTX;// = new std::vector<double>;
   std::vector<double> *gPMTY;// = new std::vector<double>;
   std::vector<double> *gPMTZ;// = new std::vector<double>;
   // pointers
//   std::vector<double>* pvtRes = &gPMTRecoTime;
//   std::vector<UInt_t>* pvfecd = &gFECD;
//   std::vector<double>* pPMTX = &gPMTX;
//   std::vector<double>* pPMTY = &gPMTY;
//   std::vector<double>* pPMTZ = &gPMTZ;

   tree->Branch("evtNum", &evtNum); 
   tree->Branch("trig", &trig);
   /// TVector3 objects  
   tree->Branch("gFitPosition",&gFitPosition);
   //tree->Branch("gTruePosition",&gTruePosition);
   //tree->Branch("gTrueDirection",&gTrueDirection);
   tree->Branch("gSolarDirection",&gSolarDirection);
   /// double values
   tree->Branch("cosTheta", &cosTheta);
   tree->Branch("gFitTime",&gFitTime);
   tree->Branch("gFecdTime",&gFecdTime);
   tree->Branch("gTrueEnergy",&gTrueEnergy);
   tree->Branch("gFitEnergy",&gFitEnergy);
   tree->Branch("gNhitClean",&gNhitClean);
   tree->Branch("itr",&itr);
   tree->Branch("distSrcToFitPos", &distSrcToFitPos); 
   /// vector containers 
   // tree->Branch("gPMTTrueTime",&gPMTTrueTime);
   tree->Branch("gPMTRecoTime",&gPMTRecoTime);
   tree->Branch("gQhs",&gQhs);
   tree->Branch("gQhl",&gQhl);

   tree->Branch("gFECD",&gFECD);
   tree->Branch("gPMTX",&gPMTX);
   tree->Branch("gPMTY",&gPMTY);
   tree->Branch("gPMTZ",&gPMTZ);
  
   //tree->Branch("gTimeResidual",&gTimeResidual,"gTimeResidual/D");
  
   char filenames[1500];
   
   Int_t i=0;
   Int_t TriggerType;
   Int_t nentries=0;
   ULong64_t bitword;
   
   while(in0>>filenames) {
    cout<<" filenames "<<filenames<<endl ;
    // Load the RAT file
    RAT::DU::DSReader dsReader( filenames);
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    size_t numEvents = dsReader.GetEntryCount();
    // Loop over all of the events
    for(size_t iEntry = iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
      //const RAT::DS::MC& rmc= rDS.GetMC();
      //const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
      gPMTRecoTime->clear();gFECD->clear();
      gQhs->clear();gQhl->clear();
      gPMTX->clear();gPMTY->clear();gPMTZ->clear();
    
      size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
      // ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
      if (iEntry%15000 == 0) {cout<<iEntry<<endl;}
      // Looping over triggered events in each ds
      for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
      {
        // Get the event
        const RAT::DS::EV& ev = rDS.GetEV(iEV);
        const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();

        // RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
        // if( RAT::EventIsClean( ev, dcAnalysisWord ) )
        {
          //if bitword==dcflags.GetFlags( dcflags.GetLatestPass() ).ToString())
          // Fill the nHit histogram, all events should have nHit
          Float_t nHitOfEvent = ev.GetNhitsCleaned();
          TriggerType = ev.GetTrigType() ;
          //cout<<" TriggerType "<<TriggerType<<endl ;
          // Now get the radius, if applicable
          Double_t radiusOfEvent = -1.;
          TVector3 positionOfEvent;
          bool directionFlag = false;
          TVector3 directionOfEvent;
          const string fitName = "partialFitter";
          ///!!! only save info for the fitted events
          //if(((TriggerType&2)==2) || ((TriggerType&18)==18))
          {
            if( ev.FitResultExists(fitName) ) {  // It needs to exist 
              if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) { // Needs to contain position
                if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) { // Needs a valid position
                  RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
                  RAT::DS::FitVertex fVertex = fResult.GetVertex(0);//!! change from GetValid to ValidPosition
                  if(fVertex.ValidPosition()) {
                    /// RAT::DS::FitVertex fVertex = fResult.GetVertex(0);// Get first fit vertex
                    evtNum = iEntry;
                    trig = ev.GetTrigType(); 
		    gFitPosition = fVertex.GetPosition();
                    gFitEnergy = 0;//fVertex.GetEnergy();
	            gFitTime = fVertex.GetTime();
                    gNhitClean = ev.GetNhitsCleaned();
                    ///MC
                    //gTruePosition = rmcparticle.GetPosition(); 
                    //gTrueEnergy = rmcparticle.GetMomentum().Mag();
                    //gTrueDirection = rmcparticle.GetMomentum().Unit();
                    /// classifiers
                    if( !ev.ClassifierResultExists("ITR:partialFitter") ) continue;
                    itr = ev.GetClassifierResult( "ITR:partialFitter" ).GetClassification( "ITR" );

                    /// PMT info
                    for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts, for 16N
                    {
                       UInt_t fecd = calpmts.GetFECDPMT(ipmt).GetID();
                       gFECD->push_back(fecd);
		       // if(fecd == 9188) 
		       gFecdTime= calpmts.GetFECDPMT(ipmt).GetTime();
		    }
		    TVector3 srcToFitPos(gFitPosition.X() - srcPos.X(), gFitPosition.Y() - srcPos.Y(), (gFitPosition.Z()-108) - srcPos.Z());
                    distSrcToFitPos = srcToFitPos.Mag();
		    for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                    {
                       //RAT::DS::PMTCal pmtCalSelect = rev.GetCalPMTs().GetPMT(ipmt);//NOTE: Javi's new selector, don't mix with calPMTs class!!!
                       //if( pmtCalStat.GetHitStatus(pmtCalSelect) == 0 )//Javi: PMT selector 
                       TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                       gPMTX->push_back(pmtpos.X()); gPMTY->push_back(pmtpos.Y()); gPMTZ->push_back(pmtpos.Z()); 
	               double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
     
     	               lightPath.CalcByPositionPartial( gFitPosition, pmtpos );
                       double distInInnerAV = lightPath.GetDistInInnerAV();
                       double distInAV = lightPath.GetDistInAV();
                       double distInWater = lightPath.GetDistInWater();
                       double distInUpperTarget = lightPath.GetDistInUpperTarget();
                       double distInLowerTarget = lightPath.GetDistInLowerTarget();
                       const double transitTime = effectiveVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );
                       double tRes = hitTime - transitTime - gFitTime;
                       double tRes_fecd = hitTime - transitTime - gFecdTime;
		       gPMTRecoTime->push_back(tRes);
                       double qhs = calpmts.GetPMT(ipmt).GetQHS();
		       double qhl = calpmts.GetPMT(ipmt).GetQHL();
                       gQhs->push_back(qhs);
                       gQhl->push_back(qhl);
		       cosTheta = (pmtpos - gFitPosition).Unit()*srcToFitPos.Unit();
	     	  /// !!! for MC
                       //lightPathMC.CalcByPositionPartial( gTruePosition, pmtpos );
                       //double distMCInInnerAV = lightPathMC.GetDistInInnerAV();
                       //double distMCInAV = lightPathMC.GetDistInAV();
                       //double distMCInWater = lightPathMC.GetDistInWater();
                       //double distMCInUpperTarget = lightPathMC.GetDistInUpperTarget();
                       //double distMCInLowerTarget = lightPathMC.GetDistInLowerTarget(); 
                       //const double transitTimeMC = effectiveVelocity.CalcByDistance( distMCInUpperTarget, distMCInAV, distMCInWater+distInLowerTarget );
                       //double tResMC = hitTime - transitTimeMC - 390 + rDS.GetMCEV(iEV).GetGTTime();
                       //gPMTTrueTime->push_back(tResMC);
	            }

                    /// solar information
                    RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
                    //day = eventTime.GetDays() ;
                    //sec = eventTime.GetSeconds() ;
                    //nsecs = eventTime.GetNanoSeconds() ;
                    TVector3 sunDir = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
                    gSolarDirection = sunDir;
               }
             } // end of loop if valid position
           }// end of loop if contains position
         // Now store the event information in an nTuple of those events with:
         // - high nHit (>=40)
         // - that have all of position, direction, and solar direction information
         
         // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
         
         // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
         tree->Fill(); 
        } //fit exist 
        } // trigWord
       } // dataClean 
      } // loop triggered events
     } // loop all events 
  } // end of the while loop for files
  // Give everything an axis label, as should always be
  
  // Write the histograms to fileName
  TString oldfilename(filenames);
  TString newfilename = "Tree_"+oldfilename;
  TFile *file=new TFile(newfilename,"RECREATE");
  file->cd();
  tree->Write();
  
}
