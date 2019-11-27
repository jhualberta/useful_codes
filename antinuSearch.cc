
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DataCleaningUtility.hh>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TNtuple.h>

#include <bitset>
#include <iostream>
#include <sstream>
#include <stdlib.h>

// To remove soon
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>

using namespace std;

// ********** All search parameter values defined here **********

int verbose = 0; // Set to 0 for no screen output
double FVcutValueNeutron = 6000.; // mm
double FVcutValuePositron = 6000.; // mm
double timeCutValue = 1000.; // us, can be less than of equal to classifier value, but not greater than it
double distanceCutValue = 2000.; // mm, not currently being applied to anything
const string fitName = "waterFitter"; // "alberta";// "waterFitter"; //"fit10"; // Name of the fitter from which to retrieve event data
bool requireGlobalValidFit = true; // Does the fit need to be globally valid
bool requireValidPosition = true; // Does the fit need to have its position valid flag, obsolete theoretically, always needed
int minNHitForFit = 6; // Any event with nHit below this cannot have a successful fit
bool keepAllToPSUP = true; // keep all events that reconstruct within an 8000. mm FV, rather than the previous FV cut, this restriction is only applied to the first event
int minNHitExtWater = 20.; // Force a minimum nHit onto the events in the external water, to keep the data set small

const double driveP0 = 0.995765, driveP1 = -63.826; // Drive corrections for 'alberta' fitter
const double c_light = 299.792458;
double waterRI = 1.38486;
double grVelocity = c_light/waterRI;//2.17554021555098529e+02;
// **************************************************************

// Code to convert to bit strings, used for converting the trigger word
std::vector< int > get_bits( unsigned long x ) {
    std::string chars( std::bitset< sizeof(long) * CHAR_BIT >( x )
        .to_string( char(0), char(1) ) );
    return std::vector< int >( chars.begin(), chars.end() );
}

struct isoVal
{
   double beta14;
   double thetaij;
};//store ISO classifier results

struct isoVal IsoClassifier(vector<TVector3>);//customized ISO classifier

// Function that returns triggered time of an event is nanoseconds
Double_t GetTriggeredEventTime(const RAT::DS::EV& tev) {
  Double_t day = tev.GetUniversalTime().GetDays() ;
  Double_t second = tev.GetUniversalTime().GetSeconds() ;
  Double_t nanosecond = tev.GetUniversalTime().GetNanoSeconds() ;
  return day*24*3600*1e9 + second*1e9 + nanosecond;  
}

// Function that gets all of the fit results from the 'alberta' fitter
void GetAlbertaFitParameters(const RAT::DS::EV& fev, const RAT::DU::PMTInfo& pmtInfo, Double_t &eventB14, TVector3 pos_cor) {
    
  // Now determine the beta14 and ITR
  const RAT::DS::CalPMTs& calpmts = fev.GetCalPMTs();
  vector<TVector3> pmtDir;
  for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++) {
    TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
    pmtDir.push_back( (pmtpos -pos_cor).Unit() );
  }//for PMT loop
  struct isoVal isoresult = IsoClassifier(pmtDir);
  eventB14 = isoresult.beta14;
  double thij = isoresult.thetaij;
 
}

// Main function 
void antinuSearch( const string& fileName, const string& outName) {
 
  // ********* Define all of the histgrams and ntuples to use as output *********
  TH1D* nHit_all = new TH1D( "nHit_all", "nHit of all triggered events", 10000, 0.5, 10000.5);  
  TH1D* radius_all = new TH1D( "radius_all", "Reconstruced radius of all triggered events", 100, 0., 15000.);      
  TH1D* nHit_withFit = new TH1D( "nHit_withFit", "nHit of all triggered events that successfully reconstruced", 10000, 0.5, 10000.5); 
  TH1D* itr_all = new TH1D( "itr_all", "ITR of all triggered events (with successful fits)", 100, 0., 1.);
  TH1D* qpdt_all = new TH1D( "qpdt_all", "QPDT of all triggered events (with successful fits)", 100, 0., 1.);
  TH1D* beta14_all = new TH1D( "beta14_all", "Beta14 of all triggered events (with successful fits)", 200, -10., 10.);  
  TH1D* timeDifference_us = new TH1D( "timeDifference_us", "Time difference between consecutive triggered events (all events)", 200, 0., 1000.);
  TH1D* timeDifferenceRecon_us = new TH1D( "timeDifferenceRecon_us", "Time difference between consecutive triggered events (only clean, reconstructed events)", 200, 0., 1000.);  
  TH1D* timeDifference_ms = new TH1D( "timeDifference_ms", "Time difference between consecutive triggered events (all events)", 200, 0., 1000.); 
  TH1D* timeDifferenceRecon_ms = new TH1D( "timeDifferenceRecon_ms", "Time difference between consecutive triggered events (only clean, reconstructed events)", 200, 0., 1000.);   
  TH1D* timeDifferenceRecon_5ms = new TH1D( "timeDifferenceRecon_5ms", "Time difference between consecutive triggered events (only clean, reconstructed events)", 200, 0., 5000.);    
  TH1D* timeDifference_inFV = new TH1D( "timeDifference_inFV", "Time difference between event pairs with both in the specified FVs", 200, 0., 1000.);
  TH3D* xyzDistribution = new TH3D("xyzDistribution", "Full position distribution of all clean events that reconstruct", 200, 0., 10000., 200, 0., 10000., 200, 0., 10000.); // 5 cm bins
  TH1D* nHit_cleanWithFit = new TH1D( "nHit_cleanWithFit", "nHit of all clean events that successfully reconstruced", 500, 0.5, 500.5);   
  TH1D* radius_clean = new TH1D( "radius_clean", "Reconstruced radius of all clean triggered events", 100, 0., 10000.);      
  TH1D* positionCut = new TH1D( "positionCut", "Distance of all considered events from the possible neutron event", 100, 0., 10000.);
  TH2D* timeDiff_positionDiff = new TH2D("timeDiff_positionDiff", "Time difference vs position difference of all considered event pairs (after FV cuts)", 200, 0., 1000., 100, 0., 10000.);
  TH2D* itr_pairs = new TH2D( "itr_pairs", "ITR of all considered event pairs", 100, 0., 1., 100., 0., 1.);
  TH2D* qpdt_pairs = new TH2D( "qpdt_pairs", "QPDT of all considered event pairs", 100, 0., 1., 100, 0., 1.);
  TH2D* beta14_pairs = new TH2D( "beta14_pairs", "Beta14 of all considered event pairs", 200, -10., 10., 200, -10., 10.);

  TH2D* recon_time_radius = new TH2D( "recon_time_radius", "Time vs radius of all clean triggered events", 100, 0., 1000., 100, 0., 15000.);
  TH2D* recon_time_nhit = new TH2D( "recon_time_nhit", "Time vs nHit of all clean events that reconstructed", 100, 0., 1000., 500, 0.5, 500.5);
  TH2D* recon_nhit_radius = new TH2D( "recon_nhit_radius", "nHit vs radius of all clean triggered events", 500, 0.5, 500.5, 100, 0., 15000.);
  
  TH1D* timeDifferenceReconNoRetrig_ms = new TH1D( "timeDifferenceReconNoRetrig_ms", "Time difference between consecutive events (only clean, reconstructed events, retriggers removed)", 200, 0., 1000.);
  TH1D* timeDifferenceNoRetrigRecon_us = new TH1D( "timeDifferenceReconNoRetrig_us", "Time difference between consecutive events (only clean, reconstructed events, retriggers removed)", 200, 0., 1000.);  
  TH1D* nHit_cleanWithFitNoRetrig = new TH1D( "nHit_cleanWithFitNoRetrig", "nHit of all clean events that successfully reconstruced, retriggers removed", 500, 0.5, 500.5);   
  TH1D* radius_cleanNoRetrig = new TH1D( "radius_cleanNoRetrig", "Reconstruced radius of all clean triggered events, retriggers removed", 100, 0., 10000.);      
  
  // ntuple to store all considered event pair info in 
  TNtuple* consideredEventPairInfo = new TNtuple("consideredEventPairInfo", "Information for every considered event pair", "runNumber:GTID_eplus:GTID_neutron:time_difference_us:position_difference_mm:nHit_eplus:nHit_neutron:radius_mm_eplus:radius_mm_neutron:x_pos_eplus:y_pos_eplus:z_pos_eplus:x_pos_neutron:y_pos_neutron:z_pos_neutron");
  TNtuple* consideredEventPairExtraInfo = new TNtuple("consideredEventPairExtraInfo", "Fit Quality information for every considered event pair", "runNumber:GTID_eplus:GTID_neutron:trigWord_eplus:trigWord_neutron:isClean_eplus:isClean_neutron::timeSinceLast_eplus::timeSinceLast_neutron::itr_eplus:qpdt_eplus:beta14_eplus:itr_neutron:qpdt_neutron:beta14_neutron");  
  TNtuple* consideredEventPairExtraExtraInfo = new TNtuple("consideredEventPairExtraExtraInfo", "More fit information for every considered event pair", "runNumber:GTID_eplus:GTID_neutron:energy_eplus:energy_neutron:x_dir_eplus:y_dir_eplus:z_dir_eplus:x_dir_neutron:y_dir_neutron:z_dir_neutron:intime_hits20_eplus:intime_hits100_eplus:intime_hits20_neutron:intime_hits100_neutron");  
    
  
  // ********* End of definition of histograms *********
    
  RAT::DU::DSReader dsReader( fileName );
  // To later get the run info
  const RAT::DS::Run& run = dsReader.GetRun();
  
  size_t numEvents = dsReader.GetEntryCount();
 
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
  Double_t previousEventTime = -1.; // Time since the last real event
  Double_t previousReconEventTime = -1.; // Time since the last (clean) event that had a successfull reconstruction

  // Counters for everything
  int numTriggeredEvents = 0;
  int numCandidatePairs = 0;
  
  ULong64_t dataCleaningWord = RAT::GetDataCleaningWord( "analysis_mask" );
  
  // Temporary counters for prescale data
  ////// int numPrescaleTrig = 0;
  ///// int numPrescaleTrigWithFit = 0; 
  
  int failedFits = 0;   
 
  int passWithQualityParam = -1; // Needed to determine which pass has the fit quality parameters
  
  // ********* Loop over all events forwards *********
  for( size_t iEntry = 0; iEntry < numEvents; iEntry++ ) {  
  
    // Get the DS for each entry (contains 0 or more triggered events [1 in the case of data])
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event  
      
    size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this DS entry
    // = 1 for data
    // 0 to N in the case of MC
    
    // Counter to display progress
    if (iEntry%10000 == 0) {
      cout<<">>> "<<iEntry<<endl;    
    }
    
    
    // ********* Looping over triggered events forwards *********
    if (iEVCount) {
      for(size_t iEV=0; iEV < iEVCount; iEV++) //iEVCount
      {
        
        // Get the event only once
        // In MC this sometimes has bugs where it tries to get an ev that doesn't exist
        try { // Attempt to get the event
        const RAT::DS::EV& ev = rDS.GetEV(iEV);
        } catch (...) {
          cout<<"Warning!!     DS entry: "<<iEntry<<" has "<<iEVCount<<" but can't access event "<<iEV<<"!"<<endl;
          goto skipToNextDS; // Bug encountered, skip this event entirely as it doesn't exist
        }
        
        // ********* Code executed for all triggered events *********
        
        Int_t currentEventTriggerWord = rDS.GetEV(iEV).GetTrigType(); // Get the trigger types
        // Determine the triggered time in nanseconds
        Double_t triggeredTime = GetTriggeredEventTime(rDS.GetEV(iEV));
        
        numTriggeredEvents++; // Incerement the total triggered event counter    
            
        // Get the fit information about this event
        Double_t currentEventnHit;
        currentEventnHit = rDS.GetEV(iEV).GetNhitsCleaned();
        // Now get the fit information
        Double_t currentEventRadius; 
        currentEventRadius = -1; // Initialize it to determine if the fit failed
        // Define the variable for energy
        Double_t currentEventEnergy = -1.; // Initialize in the event of failed fit        
        TVector3 currentEventPosition;
        TVector3 currentEventDirection;
        int currentEventGTID;
        currentEventGTID = rDS.GetEV(iEV).GetGTID();
        // First get the event's fit result
        if( rDS.GetEV(iEV).FitResultExists(fitName) ) {  // It needs to exist
          RAT::DS::FitResult fResult = rDS.GetEV(iEV).GetFitResult(fitName);// Get fit result
          if (requireGlobalValidFit) { // Check the global validity of the fit if wanted
            if (!fResult.GetValid() || currentEventnHit < minNHitForFit) {
              goto endOfFitRetrieval; // Failed, skip the it result retrieval
            }
          } // Else, continue with the fit
          RAT::DS::FitVertex fVertex = fResult.GetVertex(0); // Get the vertex 
          if( fVertex.ContainsPosition() ) { // Needs to contain position
             // If wanted check for a valid fit
             if (fitName == "alberta") {
               if (fVertex.ValidPosition()) { // Specialized 
                 TVector3 tempCurrentPosition = fVertex.GetPosition();
                 try {
                   RAT::DS::FitVertex fVertex1 = rDS.GetEV(iEV).GetFitResult("albertadirection").GetVertex(0);  
                   currentEventDirection = fVertex1.GetDirection();
                   // Do the drive correction
                   currentEventPosition = driveP0*tempCurrentPosition + driveP1*currentEventDirection;
                 } catch (...) {
                   cout<<"Event missing direction, no drive correction applied!"<<endl;
                   currentEventPosition = tempCurrentPosition;
                 }
                 currentEventRadius = currentEventPosition.Mag();
               }
             }
             else { // Do the normal waterFitter routine
               if (fVertex.ValidPosition()) {   
                 currentEventPosition = fVertex.GetPosition();
                 currentEventRadius = currentEventPosition.Mag();
               }
             }
          } // End if constains position
          // Now get the fit energy
          if( fVertex.ContainsEnergy() ) {
            if (fVertex.ValidEnergy()) {
              currentEventEnergy = fVertex.GetEnergy();  // If the fit exists and has a valid energy store it, regardless of what the position fit does
            }  
          } // End if contains energy
          
        } // End if fit exists
        
        endOfFitRetrieval:
        
        // Temporary code to see how many prescale events actually have fits
      /////  bitset<27> trigBits = bitset<27>(currentEventTriggerWord);
      ////  if ( trigBits[11] ) {
       //   numPrescaleTrig++;
       //   if (currentEventRadius >= 0.) {
       //     numPrescaleTrigWithFit++;   
       //    // cout<<currentEventTriggerWord<<endl;
      //    }
      //  }       
        
        // Define the variables for the fit quality parameters
        // Default value set to indicate failure (ie. valid check failed)
        Double_t currentEventITR;
        currentEventITR = -1;        
        Double_t currentEventQPDT;
        currentEventQPDT = -1;
        Double_t currentEventBeta14;
        currentEventBeta14 = -20;  
                      
        if (currentEventRadius < 0.) {
          // Indicates a failed fit, increment the counter
          failedFits++;
          nHit_all->Fill(currentEventnHit); // Store the nHit
        }
        else {
          // Successfully fit, store the nHit in both histograms
          nHit_all->Fill(currentEventnHit);
          nHit_withFit->Fill(currentEventnHit);
          radius_all->Fill(currentEventRadius);
          
          // Successful fit, get the fit quality parameters too
          // First determine which pass they are stored in if not already determined, 10 should be sufficient
          bool foundPass = false;
          if (passWithQualityParam < 0) {
            for (int pass = 10; pass >=0; pass--) {
              bool exists;
              try {exists = rDS.GetEV(iEV).GetClassifierResults(pass).ResultExists("ITR:"+fitName);}
              catch (...) {exists = false;}
              if (exists) {
                passWithQualityParam = pass;
                foundPass = true;
                break;
              }
            }
          } // Now the pass number has been determined (observed to be 1 for data, 0 for MC)          
          // Get the ITR, QPDT, and isotropy
          if (!foundPass) {
            passWithQualityParam=-1;   
          }
          if (passWithQualityParam >= 0) {
            const RAT::DS::ClassifierResult &itr_cResult = rDS.GetEV(iEV).GetClassifierResults(passWithQualityParam).GetResult("ITR:"+fitName);
            if (itr_cResult.GetValid()) { 
              currentEventITR = itr_cResult.GetClassification("ITR");
            }       
             const RAT::DS::ClassifierResult &qpdt_cResult = rDS.GetEV(iEV).GetClassifierResults(passWithQualityParam).GetResult("QPDT:"+fitName);
            if (qpdt_cResult.GetValid()) { 
              currentEventQPDT = qpdt_cResult.GetClassification("QPDT");
            } 
            const RAT::DS::ClassifierResult &isotropy_cResult =  rDS.GetEV(iEV).GetClassifierResults(passWithQualityParam).GetResult("isotropy:"+fitName);
            if (isotropy_cResult.GetValid()) { 
              currentEventBeta14 = isotropy_cResult.GetClassification("snobeta14");
            }                     
          }
          else if (fitName == "alberta") {
            GetAlbertaFitParameters(rDS.GetEV(iEV), pmtInfo, currentEventBeta14, currentEventPosition);   
          }
          
          // Fill the histogram for the quality parameters of all events (with successful fits)
          itr_all->Fill(currentEventITR);
          qpdt_all->Fill(currentEventQPDT);
          beta14_all->Fill(currentEventBeta14);
          
        }
        
        int isCurrentEventClean = 1; // All events presumably clean for now
  //      if( RAT::EventIsClean( ev, dataCleaningWord ) ) {
  //        isCurrentEventClean = 1;
  //      }
        
        // Looking at raw time differences between subsequent events
        Double_t absTimeSinceLastEvent = -1;
        if (previousEventTime >= 0.) { // The first event has no previous event, should not store an entry
          absTimeSinceLastEvent = -(previousEventTime-triggeredTime)/1e3;  
          timeDifference_us->Fill(absTimeSinceLastEvent); // Converting ns to us
          timeDifference_ms->Fill(-(previousEventTime-triggeredTime)/1e6); // Converting ns to ms
         // cout<<(previousEventTime-triggeredTime)/1e6<<endl;
          previousEventTime = triggeredTime;
        } else { // Initialize it for the first event
          previousEventTime = triggeredTime;  
        }
        
        // Now do the same, but only for events that reconstructed (these histograms will be used in toy MC)
        if (currentEventRadius >= 0. && isCurrentEventClean) { // Check that it's clean and has reconstructed
          if (previousReconEventTime >= 0.) {
             timeDifferenceRecon_us->Fill(-(previousReconEventTime-triggeredTime)/1e3); // Converting ns to us
             timeDifferenceRecon_ms->Fill(-(previousReconEventTime-triggeredTime)/1e6); // Converting ns to ms
             timeDifferenceRecon_5ms->Fill(-(previousReconEventTime-triggeredTime)/1e3); // Converting ns to us
             // Now the 2D histograms to check for correlations
             recon_time_radius->Fill(-(previousReconEventTime-triggeredTime)/1e3, currentEventRadius);
             recon_time_nhit->Fill(-(previousReconEventTime-triggeredTime)/1e3, currentEventnHit);
             
             previousReconEventTime = triggeredTime;
          } else { // Initialize it 
             previousReconEventTime = triggeredTime; 
          }
          // Also store it's xyx position
          xyzDistribution->Fill(currentEventPosition.X(), currentEventPosition.Y(), currentEventPosition.Z());
          // And nHit
          nHit_cleanWithFit->Fill(currentEventnHit);
          // And just radius
          radius_clean->Fill(currentEventRadius);
          // And one more 2D histogram
          recon_nhit_radius->Fill(currentEventnHit, currentEventRadius);
        }
        
        // ********* End code for all triggered events *********
          
        // ********* Begin the search for antinu candidates *********
        
        // First, only perform this search if the current event can be considered a possible positron
        // ie., has a reconstructed position, is within the positron FV
        if (true) { // Placeholder
          if (((currentEventRadius >= 0) && (currentEventRadius <= FVcutValuePositron)) || (keepAllToPSUP && (currentEventRadius >= 0) && (currentEventRadius <= 8000.) && (currentEventnHit >= minNHitExtWater))) {
              
            ostringstream outStream; // Needs to always be declared, regardless of verbosity   
            if (verbose) {
              // Store the output in a stringstream, only printing if coincidences are considered
              outStream<<"For event with GTID "<<currentEventGTID<<", nHit = "<<currentEventnHit<<", and r = "<<currentEventRadius<<" mm:"<<endl;
              outStream<<"... and quality parameters: "<<currentEventITR<<", "<<currentEventQPDT<<", "<<currentEventBeta14<<endl;  
            }
              
            // ********* Begin loop over future events, looking for pairs *********  
              
            // Start a loop looking at the next events in the data set
            // Break out of the loop as soons as the pairTimeDifference exceeds the cut value (no more possible pairs can occur)
            
            Double_t previousInnerEventTime = -1; // Time since the last event in the inner loop
            
            // Inner loop over all DS entries from this one to the end of the root file...
            for( size_t innerEntry = iEntry; innerEntry < numEvents ; innerEntry++ ) {  
                
              // Get the DS for each entry (contains 0 or more triggered events [1 in the case of data])
              const RAT::DS::Entry& rDSinner = dsReader.GetEntry( innerEntry ); // ds for each event  
      
              size_t innerEVCount = rDSinner.GetEVCount(); // Number of triggered events in this DS entry
              // = 1 for data
              // 0 to N in the case of MC 
    
              // Now loop over all of the triggered events in each DS entry again
              if (innerEVCount) {
                for(size_t innerEV=0; innerEV < innerEVCount; innerEV++) //innerEVCount
                {
                  // Get the event
                  const RAT::DS::EV& inner_ev = rDSinner.GetEV(innerEV);                      
                  // Prevent looking at events with GTIDs that are the same or higher than the current positron candidate event
                  // This will always be the case for at lease the first event here
                  int consideredEventGTID = inner_ev.GetGTID();
                  if (consideredEventGTID > currentEventGTID) {
                    
                    // Get the information about this considered event
                    Int_t consideredEventTriggerWord = inner_ev.GetTrigType(); // Get the trigger types  
                    Double_t inner_triggeredTime = GetTriggeredEventTime(inner_ev);
                    Double_t absTimeSinceLastInnerEvent = -(previousInnerEventTime-inner_triggeredTime)/1e3; // Time since actual last event 
                    previousInnerEventTime = inner_triggeredTime; // Change the time
                    Double_t pairTimeDifference = -(triggeredTime-inner_triggeredTime)/1e3; // in us
                    // If the time difference is greater than specified, immediately stop these inner loops
                    if (pairTimeDifference > timeCutValue) {
                      goto endOfInnerLoop;   
                    }
                    
                    //cout<<absTimeSinceLastEvent<<"\t"<<absTimeSinceLastInnerEvent<<"\t"<<pairTimeDifference<<endl;
                    
                    // Now get the fit information about this event
                    Double_t consideredEventnHit = inner_ev.GetNhitsCleaned();
                    // Now get the fit information
                    Double_t consideredEventRadius = -1; // Initialize it to determine if the fit failed
                    Double_t consideredEventEnergy = -1.; // Initialized
                    TVector3 consideredEventPosition;
                    TVector3 consideredEventDirection;
                    // First get the event's fit result
                    if( inner_ev.FitResultExists(fitName) ) {  // It needs to exist
                    RAT::DS::FitResult inner_fResult = inner_ev.GetFitResult(fitName); // Get fit result
                    if (requireGlobalValidFit) { // Check the global validity of the fit if wanted
                      if (!inner_fResult.GetValid() || consideredEventnHit < minNHitForFit) {
                        goto endOfSecondFitRetrieval; // Failed, skip the it result retrieval
                      }
                    }
                    RAT::DS::FitVertex inner_fVertex = inner_fResult.GetVertex(0); // Get the vertex 
                      if( inner_fVertex.ContainsPosition() ) { // Needs to contain position
                        // If wanted check for a valid fit
                        if (fitName == "alberta") { // Alberta only code
                          if (inner_fVertex.ValidPosition()) { // Check for this flag 
                            TVector3 tempConsideredPosition = inner_fVertex.GetPosition();
                            try {
                              RAT::DS::FitVertex inner_fVertex1 = inner_ev.GetFitResult("albertadirection").GetVertex(0);
                              consideredEventDirection = inner_fVertex1.GetDirection();
                              // Do the drive correction
                              consideredEventPosition = driveP0*tempConsideredPosition + driveP1*consideredEventDirection;
                            } catch (...) {
                              cout<<"Event missing direction, no drive correction applied!"<<endl;
                              consideredEventPosition = tempConsideredPosition;  
                            }
                            consideredEventRadius = consideredEventPosition.Mag();
                          }
                        }
                        else { // Regarless, run the normal waterFitte routine
                          if (inner_fVertex.ValidPosition()) { 
                            consideredEventPosition = inner_fVertex.GetPosition();
                            consideredEventRadius = consideredEventPosition.Mag();
                          }  
                        }
                      } // End if constains position
                      // Now get the fit energy
                      if( inner_fVertex.ContainsEnergy() ) {
                        if (inner_fVertex.ValidEnergy()) {
                          consideredEventEnergy = inner_fVertex.GetEnergy();  // If the fit exists and has a valid energy store it, regardless of what the position fit does
                        }  
                      } // End if contains energy
          
                    } // End if fit exists
                    
                    endOfSecondFitRetrieval:
                    
                    // ********* Begin application of some antineutrino criteria *********  
                        
                    if (verbose) {
                      outStream<<"  Considering GTID = "<<consideredEventGTID<<" with delta_t = "<<pairTimeDifference<<" us, nHit = "<<consideredEventnHit<<", and r = "<<consideredEventRadius<<" mm"<<endl;   
                    }
                      
                    // First require that the candidate neutron event have a successful fit and is within the specified FV
                    if ((consideredEventRadius >= 0) && (consideredEventRadius <= FVcutValueNeutron) || (keepAllToPSUP && (consideredEventRadius >= 0) && (consideredEventRadius <= 8000.))) {
                      
                      // Store time difference information
                      timeDifference_inFV->Fill(pairTimeDifference);
                        
                      // Calculate the distance between the reconstruced positions of the two events
                      Double_t pairPositionDifference = (consideredEventPosition-currentEventPosition).Mag();
                      // And store it
                      positionCut->Fill(pairPositionDifference);
                      // And now store both
                      timeDiff_positionDiff->Fill(pairTimeDifference, pairPositionDifference);
                      
                      // Now take a look at some fit quality parameters
                      // Get the ITR, QPDT, and isotropy
                      Double_t consideredEventITR = -1;
                      Double_t consideredEventBeta14 = -20;
                      Double_t consideredEventQPDT = -1;
                      
                      if (passWithQualityParam >= 0) {
                        const RAT::DS::ClassifierResult &inner_itr_cResult = inner_ev.GetClassifierResults(passWithQualityParam).GetResult("ITR:"+fitName);
                        if (inner_itr_cResult.GetValid()) { 
                          consideredEventITR = inner_itr_cResult.GetClassification("ITR");
                        }       
                      
                        const RAT::DS::ClassifierResult &inner_qpdt_cResult = inner_ev.GetClassifierResults(passWithQualityParam).GetResult("QPDT:"+fitName);
                        if (inner_qpdt_cResult.GetValid()) { 
                          consideredEventQPDT = inner_qpdt_cResult.GetClassification("QPDT");
                        } 
                      
                        const RAT::DS::ClassifierResult &inner_isotropy_cResult = inner_ev.GetClassifierResults(passWithQualityParam).GetResult("isotropy:"+fitName);
                        if (inner_isotropy_cResult.GetValid()) { 
                          consideredEventBeta14 = inner_isotropy_cResult.GetClassification("snobeta14");
                        }               
                      }
                      else if (fitName == "alberta") {
                        GetAlbertaFitParameters(inner_ev, pmtInfo, consideredEventBeta14, consideredEventPosition);   
                      }
                      
                      // Now fill the 2D histograms
                      itr_pairs->Fill(consideredEventITR, currentEventITR);
                      qpdt_pairs->Fill(consideredEventQPDT, currentEventQPDT);
                      beta14_pairs->Fill(consideredEventBeta14, currentEventBeta14);
                      
                      if (verbose) {
                        outStream<<"    Quality parameters: "<<consideredEventITR<<", "<<consideredEventQPDT<<", "<<consideredEventBeta14<<endl;  
                        outStream<<"      Position difference = "<<pairPositionDifference<<" mm"<<endl;   
                      }
                      
                      int isConsideredEventClean = 1; // All considered clean for now
                      // Determine whether these events would pass the analysis_mask data cleaning flags

                      //if( RAT::EventIsClean( inner_ev, dataCleaningWord ) ) {
                      //  isConsideredEventClean = 1;
                      //}
                        
                      // Get some more information about the events
                      double currentEventInTimeHits100 = rDS.GetEV(iEV).GetInTimeHits100();
                      double currentEventInTimeHits20 = rDS.GetEV(iEV).GetInTimeHits20();
                      double consideredEventInTimeHits100 = inner_ev.GetInTimeHits100();
                      double consideredEventInTimeHits20 = inner_ev.GetInTimeHits20();                    
                        
                      // Now store all of the information in the ntuples that can be used in a later analysis
                      consideredEventPairInfo->Fill(run.GetRunID(), currentEventGTID, consideredEventGTID, pairTimeDifference, pairPositionDifference, currentEventnHit, consideredEventnHit, currentEventRadius, consideredEventRadius, currentEventPosition.X(), currentEventPosition.Y(), currentEventPosition.Z(), consideredEventPosition.X(), consideredEventPosition.Y(), consideredEventPosition.Z());
                      consideredEventPairExtraInfo->Fill(run.GetRunID(), currentEventGTID, consideredEventGTID, consideredEventTriggerWord, currentEventTriggerWord, isConsideredEventClean, isCurrentEventClean, absTimeSinceLastEvent, absTimeSinceLastInnerEvent, currentEventITR, currentEventQPDT, currentEventBeta14, consideredEventITR, consideredEventQPDT, consideredEventBeta14);
                      consideredEventPairExtraExtraInfo->Fill(run.GetRunID(), currentEventGTID, consideredEventGTID, currentEventEnergy, consideredEventEnergy, currentEventDirection.X(), currentEventDirection.Y(), currentEventDirection.Z(), consideredEventDirection.X(), consideredEventDirection.Y(), consideredEventDirection.Z(), currentEventInTimeHits20, currentEventInTimeHits100, consideredEventInTimeHits20, consideredEventInTimeHits100);
                        
                    } // End if successfully fit and is within the positron FV
                      
                    // Print the screen output if wanted for this event's considerations
                    if (verbose) {
                      cout<<outStream.str();   
                    }
                      
                    // ********* End application of antineutrino criteria  *********
                      
                  }// End if the GTID is less than the current event
                  else if (consideredEventGTID == currentEventGTID) {
                    previousInnerEventTime = GetTriggeredEventTime(inner_ev); // Initialize the time for the 1st inner loop, the next event should have code executed
                  }
                   
                } // End of inner loop over triggered events
                                
              } // End of if there are triggered events in the inner loops  
                
            } // End of inner loop over all DS entries
        
    
            // Label to end the inner loops, jumps to here
            endOfInnerLoop:
            
            // Do nothing, just need a pace for this label to go to
            if (verbose) {}
              
          } // End of radius exists and less than neutron FV cut
        } // Placeholder
        
        // ********* End of search for antinu candidates *********  
          
      } // End of main loop over trig events    
    } // End of if non-zero triggered events
    // ********* End forward loop over triggered events *********
    
    skipToNextDS:
    // Do nothing, just need a pace for this label to go to
    if (verbose) {}
      
  } // ********* End forward loop over all events ********* 

  // Display the end of run statistics
  cout<<"Total triggered events: "<<numTriggeredEvents<<endl;
  cout<<"Failed fits: "<<failedFits<<endl;
  cout<<"Failure rate: "<<Double_t(failedFits)/Double_t(numTriggeredEvents)*100<<" %"<<endl;
  
  // Temporary code for prescale counters
////  cout<<"Total events with prescale trigger: "<<numPrescaleTrig<<endl;
///  cout<<"  Number with successful fits: "<<numPrescaleTrigWithFit<<endl;  
  
  // ********* Label and write out the histograms and ntuples ********* 
  
  nHit_all->SetXTitle( "nHit" );
  nHit_all->SetYTitle( "Counts" );
  radius_all->SetXTitle( "Reconstructed radius [mm]" );
  radius_all->SetYTitle( "Counts" );  
  radius_clean->SetXTitle( "Reconstructed radius [mm]" );
  radius_clean->SetYTitle( "Counts" );   
  nHit_withFit->SetXTitle( "nHit" );
  nHit_withFit->SetYTitle( "Counts" ); 
  nHit_cleanWithFit->SetXTitle( "nHit" );
  nHit_cleanWithFit->SetYTitle( "Counts" );  
  itr_all->SetXTitle( "ITR" );
  itr_all->SetYTitle( "Counts" );
  qpdt_all->SetXTitle( "QPDT" );
  qpdt_all->SetYTitle( "Counts" );
  beta14_all->SetXTitle( "Beta14" );
  beta14_all->SetYTitle( "Counts" );  
  timeDifference_us->SetXTitle( "Time difference [us]" );
  timeDifference_us->SetYTitle( "Counts" );  
  timeDifferenceRecon_us->SetXTitle( "Time difference [us]" );
  timeDifferenceRecon_us->SetYTitle( "Counts" );    
  timeDifference_ms->SetXTitle( "Time difference [ms]" );
  timeDifference_ms->SetYTitle( "Counts" );  
  timeDifferenceRecon_ms->SetXTitle( "Time difference [ms]" );
  timeDifferenceRecon_ms->SetYTitle( "Counts" );  
  timeDifferenceRecon_5ms->SetXTitle( "Time difference [us]" );
  timeDifferenceRecon_5ms->SetYTitle( "Counts" );    
  timeDifference_inFV->SetXTitle( "Time difference [us]" );
  timeDifference_inFV->SetYTitle( "Counts" );
  xyzDistribution->SetXTitle( "X [mm]" );
  xyzDistribution->SetYTitle( "Y [mm]" );
  xyzDistribution->SetZTitle( "Z [mm]" );
  positionCut->SetXTitle( "Difference in reconstructed position [mm]" );
  positionCut->SetYTitle( "Counts" );
  timeDiff_positionDiff->SetXTitle( "Time difference [us]" );
  timeDiff_positionDiff->SetYTitle( "Position difference [mm]" );
  timeDiff_positionDiff->SetZTitle( "Counts" );    
  itr_pairs->SetXTitle( "e+ ITR" );
  itr_pairs->SetYTitle( "n ITR" );
  itr_pairs->SetZTitle( "Counts" );    
  qpdt_pairs->SetXTitle( "e+ QPDT" );
  qpdt_pairs->SetYTitle( "n QPDT" );
  qpdt_pairs->SetZTitle( "Counts" ); 
  beta14_pairs->SetXTitle( "e+ Beta14" );
  beta14_pairs->SetYTitle( "n Beta14" );
  beta14_pairs->SetZTitle( "Counts" ); 
  recon_time_radius->SetXTitle( "Time difference [us]" );
  recon_time_radius->SetYTitle( "Reconstructed radius [mm]" );
  recon_time_radius->SetZTitle( "Counts" ); 
  recon_time_nhit->SetXTitle( "Time difference [us]" );
  recon_time_nhit->SetYTitle( "nHit" );
  recon_time_nhit->SetZTitle( "Counts" ); 
  recon_nhit_radius->SetXTitle( "nHit" );
  recon_nhit_radius->SetYTitle( "Reconstructed radius [mm]" );
  recon_nhit_radius->SetZTitle( "Counts" );      
   
  TFile *file=new TFile(outName.c_str(),"RECREATE");
  //TFile *file=new TFile(outName,"RECREATE");
  file->cd();
  nHit_all->Write();
  radius_all->Write();
  radius_clean->Write();
  nHit_withFit->Write();
  nHit_cleanWithFit->Write();  
  itr_all->Write();
  qpdt_all->Write();
  beta14_all->Write();  
  timeDifference_us->Write();
  timeDifferenceRecon_us->Write();
  timeDifference_ms->Write();
  timeDifferenceRecon_ms->Write();
  timeDifferenceRecon_5ms->Write();
  timeDifference_inFV->Write();
  xyzDistribution->Write();
  positionCut->Write();  
  timeDiff_positionDiff->Write();
  itr_pairs->Write();
  qpdt_pairs->Write();
  beta14_pairs->Write();
  recon_time_radius->Write();
  recon_time_nhit->Write();
  recon_nhit_radius->Write();
  consideredEventPairInfo->Write();
  consideredEventPairExtraInfo->Write();
  consideredEventPairExtraExtraInfo->Write();    
  
  file->Close();
  
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
