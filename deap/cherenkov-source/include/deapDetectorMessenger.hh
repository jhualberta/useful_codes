// $Id: deapDetectorMessenger.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file deapDetectorMessenger.hh
/// \brief Definition of the deapDetectorMessenger class

#ifndef deapDetectorMessenger_h
#define deapDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class deapDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for deapDetectorConstruction.
///
/// It implements commands:
/// - /B2/det/setTargetMaterial name
/// - /B2/det/setChamberMaterial name
/// - /B2/det/stepMax value unit

class deapDetectorMessenger: public G4UImessenger
{
  public:
    deapDetectorMessenger(deapDetectorConstruction* );
    virtual ~deapDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    deapDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fDeapDirectory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fTargMatCmd;
    //G4UIcmdWithAString*      fChamMatCmd;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
