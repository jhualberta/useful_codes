// $Id: deapEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file deapEventAction.cc
/// \brief Implementation of the deapEventAction class

#include "deapEventAction.hh"
#include "deapRunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "deapTrackerHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapEventAction::deapEventAction(deapRunAction* runAction)
: G4UserEventAction()//, fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapEventAction::~deapEventAction()
{
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapEventAction::BeginOfEventAction(const G4Event*)
{
  //fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapEventAction::EndOfEventAction(const G4Event* event)
{
  // accumulate statistics in run action
  //fRunAction->AddEdep(fEdep);

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  /// periodic printing

  G4int eventID = event->GetEventID();
  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
    }
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4cout << "    "
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }

/// added to save root file
  auto evthc = static_cast<deapTrackerHitsCollection*>(event->GetHCofThisEvent()->GetHC(0));
  auto analysisManager = G4AnalysisManager::Instance();
  G4int         ftrackID;
  G4double      fstepLength;
  G4double      fedep;
  G4ThreeVector fpos;

  //G4double EdepHit0 = 0.0;
  if (evthc->GetSize() != 0){
     ftrackId = (*evthc)[0]->GetTrackID();
     fstepLength = (*evthc)[0]->GetStepLength();
     fedep = (*evthc)[0]->GetEdep();
     fpos  = (*evthc)[0]->GetPos();
     G4cout<< EdepHit0 <<G4endl;
  }

//fTrackID, Edep, StepLength, Pos
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, ftrackID);
  analysisManager->FillNtupleDColumn(1, fstepLength);
  analysisManager->FillNtupleDColumn(2, fedep);
  analysisManager->FillNtupleDColumn(3, fpos);
  analysisManager->AddNtupleRow();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
