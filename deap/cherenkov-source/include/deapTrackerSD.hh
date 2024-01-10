// $Id: deapTrackerSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file deapTrackerSD.hh
/// \brief Definition of the deapTrackerSD class

#ifndef deapTrackerSD_h
#define deapTrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "deapTrackerHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// deapTracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class deapTrackerSD : public G4VSensitiveDetector
{
  public:
    deapTrackerSD(const G4String& name, 
                const G4String& hitsCollectionName);
    virtual ~deapTrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    deapTrackerHitsCollection* fHitsCollection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
