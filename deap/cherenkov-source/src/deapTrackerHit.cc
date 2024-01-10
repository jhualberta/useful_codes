// $Id: deapTrackerHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file deapTrackerHit.cc
/// \brief Implementation of the deapTrackerHit class

#include "deapTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<deapTrackerHit>* deapTrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapTrackerHit::deapTrackerHit()
 : G4VHit(),
   fTrackID(-1),
   fEdep(0.),
   fPos(G4ThreeVector()),
   fTrackLength(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapTrackerHit::~deapTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapTrackerHit::deapTrackerHit(const deapTrackerHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fTrackLength = right.fTrackLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const deapTrackerHit& deapTrackerHit::operator=(const deapTrackerHit& right)
{
  fTrackID   = right.fTrackID;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fTrackLength = right.fTrackLength;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int deapTrackerHit::operator==(const deapTrackerHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void deapTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapTrackerHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID
     << "Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
