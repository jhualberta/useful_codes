// $Id: deapEventAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file deapEventAction.hh
/// \brief Definition of the deapEventAction class

#ifndef deapEventAction_h
#define deapEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class deapRunAction;

/// Event action class
///

class deapEventAction : public G4UserEventAction
{
  public:
    deapEventAction(deapRunAction* runAction);
    virtual ~deapEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }

  private:
    deapRunAction* fRunAction;
    G4double     fEdep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
