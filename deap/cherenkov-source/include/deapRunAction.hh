// $Id: deapRunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file deapRunAction.hh
/// \brief Definition of the deapRunAction class

#ifndef deapRunAction_h
#define deapRunAction_h 1

#include "G4UserRunAction.hh"
//#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed MeV is then printed on the screen.

class deapRunAction : public G4UserRunAction
{
  public:
    deapRunAction();
    virtual ~deapRunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    //void AddEdep (G4double edep); 

  //private:
  //  G4Accumulable<G4double> fEdep;
  //  G4Accumulable<G4double> fEdep2;
};

#endif

