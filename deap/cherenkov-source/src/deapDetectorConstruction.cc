// $Id: deapDetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file deapDetectorConstruction.cc
/// \brief Implementation of the deapDetectorConstruction class

#include "deapDetectorConstruction.hh"
#include "deapDetectorMessenger.hh"
#include "deapTrackerSD.hh"
#include "G4SDManager.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
//#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
//#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapDetectorConstruction::deapDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

deapDetectorConstruction::~deapDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* deapDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 780*cm, env_sizeZ = 780*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.03*env_sizeXY;
  G4double world_sizeZ  = 1.03*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.9*world_sizeXY, 0.9*world_sizeXY, 0.9*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.9*env_sizeXY, 0.9*env_sizeXY, 0.9*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 

  // Shape 1, acrylic shell
  G4Material* acrylic = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  //G4Element *elC = nist->FindOrBuildElement("C");
  //G4Element *elH = nist->FindOrBuildElement("H");
  //G4Element *elO = nist->FindOrBuildElement("O");
  //acrylic->AddElement(elC, 5);
  //acrylic->AddElement(elH, 8);
  //acrylic->AddElement(elO, 2);
///  G4Material* vacuum = nist->FindOrBuildMaterial("G4_VACUUM");

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("gasArgon", z=18., a=39.95*g/mole, density = 0.001784*g/cm3);
  // gas Ar 0.001784*g/cm3 at room temperature; 0.004915 g/cm3 at 100K 
  G4Material* medium = G4Material::GetMaterial("liquidArgon");
  G4Material* mediumTop = G4Material::GetMaterial("gasArgon");
  G4ThreeVector pos1 = G4ThreeVector(0*cm, 0*cm, 0*cm);

  // Dimensions of the AV shell
  G4double shellInnerRadius = 85.0*cm;
  G4double shellOuterRadius = shellInnerRadius + 50*cm;
  G4double shellThetaMin = 0.*deg;
  G4double shellThetaMax = 360.*deg;
  G4double shellPhiMin = 0.*deg;
  G4double shellPhiMax = 360.*deg;

//// Dimensions of the inner AV, medium bottom, liquid Ar
//// Sphere geometry
//  G4double innerRadius = 0*cm;
//  G4double outerRadius = 85.0*cm;
//  G4double shellThetaMin1 = 0.*deg;
//  G4double shellThetaMax1 = 58.88*deg;
//
//  G4double shellThetaMin2 = 58.88*deg;
//  G4double shellThetaMax2 = 180.*deg;
  
  // G4double shellHalfHeight = 0.5 * 57.6*cm;

  // Create the shell solid, thickness = 50 cm
  G4Sphere* shellSolid = new G4Sphere("AVshell", shellInnerRadius, shellOuterRadius, shellPhiMin, shellPhiMax, shellThetaMin, shellThetaMax);
  // Create logical volume for the shell
  G4LogicalVolume* shellLogical = new G4LogicalVolume(shellSolid, acrylic, "AVshell");
  // Create visualization attributes for the shell
  //G4VisAttributes* shellVisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue color
  //shellLogical->SetVisAttributes(shellVisAttributes);
  new G4PVPlacement(0, pos1, shellLogical, "AVshell", logicEnv, false, 0, checkOverlaps);

  // Create the shell solid
  //
  G4double r1 = 85.0*cm;
  // cap over a sphere. liquid level = 55.1 cm, then H = 85 - 55.1 = 29.9 cm
  G4Ellipsoid* innerAVtop = new G4Ellipsoid("innerAVtop", r1, r1, r1, r1-29.9*cm, r1);
  G4Ellipsoid* innerAVbot = new G4Ellipsoid("innerAVbot", r1, r1, r1, -r1, r1-29.9*cm);  

//  G4Sphere* innerAVtop = new G4Sphere("innerAVtop", innerRadius, outerRadius, shellPhiMin, shellPhiMax, shellThetaMin1, shellThetaMax1);
// Create logical volume for the shell
  G4LogicalVolume* shellLogical2 = new G4LogicalVolume(innerAVtop, mediumTop, "innerAVtop");
  // Create visualization attributes for the shell
  //G4VisAttributes* shellVisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue color
  //shellLogical->SetVisAttributes(shellVisAttributes);
  new G4PVPlacement(0, pos1, shellLogical2, "innerAVtop", logicEnv, false, 0, checkOverlaps);

//  G4Sphere* innerAVbot = new G4Sphere("innerAVbot", innerRadius, outerRadius, shellPhiMin, shellPhiMax, shellThetaMin2, shellThetaMax2);
 
  G4LogicalVolume* shellLogical3 = new G4LogicalVolume(innerAVbot, medium, "innerAVbot");
  new G4PVPlacement(0, pos1, shellLogical3, "innerAVbot", logicEnv, false, 0, checkOverlaps);

  // Create placement for the shell in the world volume
  //G4VPhysicalVolume* shellPhysical = 
  //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), shellLogical, "AVshell", logicEnv, false, 0, checkOverlaps);

  fScoringVolume = shellLogical;

  G4String trackerAVshell = "AVshell";
  deapTrackerSD* aTrackerSD = new deapTrackerSD(trackerAVshell,"TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);

  return physWorld;
  
  //  return shellPhysical;

//  new G4PVPlacement(0,                       //no rotation
//                    pos1,                    //at position
//                    logicShape1,             //its logical volume
//                    "AVshell",                //its name
//                    logicEnv,                //its mother  volume
//                    false,                   //no boolean operation
//                    0,                       //copy number
//                    checkOverlaps);          //overlaps checking

/*
  //     
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape       
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;      
  G4Trd* solidShape2 =    
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb, 
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape2;
*/
  //
  //always return the physical World
  //
  // return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
