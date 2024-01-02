// $Id: deapPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file deapPrimaryGeneratorAction.cc
/// \brief Implementation of the deapPrimaryGeneratorAction class

#include "deapPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapPrimaryGeneratorAction::deapPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);

  //dir9 = [-0.0153, -0.3919, 0]
  //dir24 = [-0.1876, -0.5774, 0]
  //dir30 = [-0.3972, -0.5889, 0]
  //dir45 = [-0.5846, -0.5645, 0]
  //dir55 = [-0.428, -0.7414, 0]
  //dir160 = [-0.6393,-0.7100, 0]
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-0.0153, -0.3919, 0.0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-0.1876, -0.5774, 0.0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-0.397,-0.589,0.));//NOTE! (-x,-y, 0)
  fParticleGun->SetParticleMomentumDirection( G4ThreeVector(-0.6393,-0.7100, 0.) );
  fParticleGun->SetParticleEnergy(2.6*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapPrimaryGeneratorAction::~deapPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void deapPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("deapPrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }
  // pmtid9 = [13.67, 349.92, 821.34]
  // pmtid24 = [168.84, 519.66, 715.14]
  // pmtid30 = [357.48, 530.01, 633.51]
  // pmtid45 = [526.14, 508.05, 524.43]
  // pmtid55 = [385.2, 667.26, 465.12]
  // pmtid160 = [575.37 , 639.0 , -265.68]
  //G4double size = 0.8;
  //id9
  //G4double x0 = 13.67/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = 349.92/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = 821.34/10*cm;//-0.5 * envSizeZ;

  //id24 
  //G4double x0 = 168.84/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = 519.66/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = 715.14/10*cm;//-0.5 * envSizeZ;

  ////id30
  G4double x0 = 357.48/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = 530.01/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = 633.51/10*cm;//-0.5 * envSizeZ;

  //id45
  //G4double x0 = 13.67/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = 349.92/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = 821.34/10*cm;//-0.5 * envSizeZ;

  //id55
  //G4double x0 = 13.67/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = 349.92/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = 821.34/10*cm;//-0.5 * envSizeZ;
 
  //id160
  //G4double x0 = 575.37/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = 639.0/10*cm;//size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = -265.68/10*cm;//-0.5 * envSizeZ;
  
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
