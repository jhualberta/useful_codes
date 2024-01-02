// $Id: deapEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file deapEventAction.cc
/// \brief Implementation of the deapEventAction class

#include "deapEventAction.hh"
#include "deapRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapEventAction::deapEventAction(deapRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
  //Instantiate analysis manager
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); //using ROOT
  //auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel( 1 );
  analysisManager->SetNtupleMerging( 1 );
  analysisManager->CreateNtuple("tree", "tree");
  analysisManager->CreateNtupleDColumn("EdepHit0");
  analysisManager->FinishNtuple();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapEventAction::~deapEventAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapEventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapEventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
