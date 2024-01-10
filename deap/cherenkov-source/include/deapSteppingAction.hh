// $Id: deapSteppingAction.hh 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file deapSteppingAction.hh
/// \brief Definition of the deapSteppingAction class

#ifndef deapSteppingAction_h
#define deapSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class deapEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class deapSteppingAction : public G4UserSteppingAction
{
  public:
    deapSteppingAction(deapEventAction* eventAction);
    virtual ~deapSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    deapEventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
