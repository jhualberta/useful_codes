// $Id: deapActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file deapActionInitialization.cc
/// \brief Implementation of the deapActionInitialization class

#include "deapActionInitialization.hh"
#include "deapPrimaryGeneratorAction.hh"
#include "deapRunAction.hh"
#include "deapEventAction.hh"
#include "deapSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapActionInitialization::deapActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapActionInitialization::~deapActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapActionInitialization::BuildForMaster() const
{
  deapRunAction* runAction = new deapRunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapActionInitialization::Build() const
{
  SetUserAction(new deapPrimaryGeneratorAction);

  deapRunAction* runAction = new deapRunAction;
  SetUserAction(runAction);

  deapEventAction* eventAction = new deapEventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new deapSteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
