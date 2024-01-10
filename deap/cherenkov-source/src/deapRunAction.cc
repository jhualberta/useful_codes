// $Id: deapRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file deapRunAction.cc
/// \brief Implementation of the deapRunAction class

#include "G4AnalysisManagerState.hh"
#include "deapRunAction.hh"
// #include "AnalysisManager.hh"
#include "deapPrimaryGeneratorAction.hh"
#include "deapDetectorConstruction.hh"
// #include "deapRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "g4root.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapRunAction::deapRunAction()
: G4UserRunAction()//, fEdep(0.)
{
  //add new units for dose
  // 
  //  const G4double milligray = 1.e-3*gray;    
  //  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  // Register accumulable to the accumulable manager
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->RegisterAccumulable(fEdep);
  //accumulableManager->RegisterAccumulable(fEdep2);
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
  //Instantiate analysis manager
  //
  auto analysisManager = G4AnalysisManager::Instance(); //using ROOT
  analysisManager->SetVerboseLevel( 1 );
  analysisManager->SetNtupleMerging( 1 );
  analysisManager->CreateNtuple("T", "T");
  analysisManager->CreateNtupleDColumn("trackID");
  analysisManager->CreateNtupleDColumn("stepLength");
  analysisManager->CreateNtupleDColumn("edep");
  analysisManager->CreateNtupleDColumn("position");

  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapRunAction::~deapRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapRunAction::BeginOfRunAction(const G4Run*)
{
  // reset accumulables to their initial values
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Reset();

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //Inform G4AnalysisManager of output file
  //
  auto analysisManager = G4AnalysisManager::Instance();
  G4String fileName = "outCheren";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
/*
  // Merge accumulables 
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  //G4double edep  = fEdep.GetValue();
  //G4double edep2 = fEdep2.GetValue();
  
  //G4double rms = edep2 - edep*edep/nofEvents;
  //if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const deapDetectorConstruction* detectorConstruction
   = static_cast<const deapDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  //G4double dose = edep/mass;
  //G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const deapPrimaryGeneratorAction* generatorAction
   = static_cast<const deapPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
  //Inform G4AnalysisManager its time to end
  //
*/
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Additional actions
//void deapRunAction::AddEdep(G4double edep)
//{
//
//  //fEdep  += edep;
//  //fEdep2 += edep*edep;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
