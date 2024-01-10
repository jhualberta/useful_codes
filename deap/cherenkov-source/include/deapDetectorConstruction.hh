// $Id: deapDetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file deapDetectorConstruction.hh
/// \brief Definition of the deapDetectorConstruction class

#ifndef deapDetectorConstruction_h
#define deapDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class deapDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    deapDetectorConstruction();
    virtual ~deapDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume*  fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
