// $Id: deapTrackerHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file deapTrackerHit.hh
/// \brief Definition of the deapTrackerHit class

#ifndef deapTrackerHit_h
#define deapTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fStepLength, fEdep, fPos

class deapTrackerHit : public G4VHit
{
  public:
    deapTrackerHit();
    deapTrackerHit(const deapTrackerHit&);
    virtual ~deapTrackerHit();

    // operators
    const deapTrackerHit& operator=(const deapTrackerHit&);
    G4int operator==(const deapTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; };
    //void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetStepLength (G4double dL)      { fStepLength = dL; };
    void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    //G4int GetChamberNb() const   { return fChamberNb; };
    G4double GetStepLength() const { return fStepLength; };
    G4double GetEdep() const     { return fEdep; };
    G4ThreeVector GetPos() const { return fPos; };

  private:

      G4int         fTrackID;
      //G4int         fChamberNb;
      G4double      fStepLength;
      G4double      fEdep;
      G4ThreeVector fPos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<deapTrackerHit> deapTrackerHitsCollection;

extern G4ThreadLocal G4Allocator<deapTrackerHit>* deapTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* deapTrackerHit::operator new(size_t)
{
  if(!deapTrackerHitAllocator)
      deapTrackerHitAllocator = new G4Allocator<deapTrackerHit>;
  return (void *) deapTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void deapTrackerHit::operator delete(void *hit)
{
  deapTrackerHitAllocator->FreeSingle((deapTrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
