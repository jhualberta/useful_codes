#include <CLHEP/Units/PhysicalConstants.h>
using namespace CLHEP;

#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4ThreeVector.hh>
#include <G4Color.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>

#include <RAT/GeoN16SNOMANSourceFactory.hh>
#include <RAT/Materials.hh>
#include <RAT/Detector.hh>
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <RAT/EnvelopeConstructor.hh>
#include <RAT/PMTConstructorParams.hh>
#include <RAT/CalibPMTSD.hh>

#include <vector>
#include <string>

namespace RAT
{

void
GeoN16SNOMANSourceFactory::Construct( DBLinkPtr table,
                                      const bool checkOverlaps)
{
  std::string prefix = table->GetS( "index" );
  std::string motherName = table->GetS( "mother" );
  int lcn = table->GetI( "lcn" );

  G4Color stainlessSteelCol(0.0,1.0,0.0);
  G4Color acrylicSnoCol(1.0,0.0,0.0);
  G4Color airCol(1.0,1.0,1.0);
  G4Color plasticScintCol(0.0,1.0,1.0);
  G4Color polyethyleneCol(0.0,0.0,1.0);
  G4Color scintCol(1.0,1.0,0.0);

  // First make the physical + logical volumes, from the outside in
  // GEDS 160
  G4VSolid* n16Source = new G4Tubs( prefix + "_solid", 0.0, 57.15, 254.0, 0.0, twopi );
  G4LogicalVolume* n16SourceLog = new G4LogicalVolume( n16Source, G4Material::GetMaterial( "stainless_steel" ), prefix + "_logic" );
  n16SourceLog->SetVisAttributes( stainlessSteelCol );
  // GEDS 161, Redundant if n16Source is same material, kept to match SNOMAN
  G4VSolid* can = new G4Tubs( prefix + "_can_solid", 0.0, 52.39, 219.3, 0.0, twopi );
  G4LogicalVolume* canLog = new G4LogicalVolume( can, G4Material::GetMaterial( "stainless_steel" ), prefix + "_can_logic" );
  canLog->SetVisAttributes( stainlessSteelCol );
  // GEDS 162
  G4VSolid* window = new G4Tubs( prefix + "_window_solid", 0.0, 30.48, 9.5, 0.0, twopi );
  G4LogicalVolume* windowLog = new G4LogicalVolume( window, G4Material::GetMaterial( "acrylic_sno" ), prefix + "_window_logic" );
  windowLog->SetVisAttributes( acrylicSnoCol );
  // GEDS 163
  G4VSolid* gap = new G4Tubs( prefix + "_gap_solid", 0.0, 52.39, 80.35, 0.0, twopi );
  G4LogicalVolume* gapLog = new G4LogicalVolume( gap, G4Material::GetMaterial( "air" ), prefix + "_air_logic" );
  gapLog->SetVisAttributes( airCol );
  // GEDS 164
  G4VSolid* scint = new G4Tubs( prefix + "_scint_solid", 0.0, 50.8, 79.4, 0.0, twopi );
  G4LogicalVolume* scintLog = new G4LogicalVolume( scint, G4Material::GetMaterial( "G4_PLASTIC_SC_VINYLTOLUENE" ), prefix + "_scint_logic" ); // TODO PHIL Fix Scint material - (MATT: possibly correct now)
  // Now make the scintillator sensitive (this is the PMT aspect)
  G4SDManager* sDManager = G4SDManager::GetSDMpointer();
  CalibPMTSD* pmtSD = new CalibPMTSD( table->GetS( "sensitive_detector" ), lcn, 0.1 * MeV, 1.0 ); ///100% efficient?
  sDManager->AddNewDetector( pmtSD );
  scintLog->SetSensitiveDetector( pmtSD );
  scintLog->SetVisAttributes( plasticScintCol );
  // GEDS 165
  G4VSolid* gas = new G4Tubs( prefix + "_gas_solid", 0.0, 47.8, 76.2, 0.0, twopi );
  G4LogicalVolume* gasLog = new G4LogicalVolume( gas, G4Material::GetMaterial( "air" ), prefix + "_gas_logic" );
  gasLog->SetVisAttributes( airCol );
  // GEDS 166
  G4VSolid* instr = new G4Tubs( prefix + "_instr_solid", 0.0, 52.39, 129.35, 0.0, twopi );
  G4LogicalVolume* instrLog = new G4LogicalVolume( instr, G4Material::GetMaterial( "air" ), prefix + "_instr_logic" );
  instrLog->SetVisAttributes( airCol );

  // GEDS 167
  G4VSolid* pmtCase = new G4Tubs( prefix + "_pmt_case_solid", 0.0 * mm, 27.71 * mm, 76.2 * mm, 0.0, twopi );
  G4LogicalVolume* pmtCaseLog = new G4LogicalVolume( pmtCase, G4Material::GetMaterial( "aluminum" ), prefix + "_pmt_case_logic" );
  // GEDS 168
  G4VSolid* pmtVac = new G4Tubs( prefix + "_pmt_vac_solid", 0.0 * mm, 25.71 * mm, 74.2 * mm, 0.0, twopi );
  G4LogicalVolume* pmtVacLog = new G4LogicalVolume( pmtVac, G4Material::GetMaterial( "pmt_vacuum" ), prefix + "_pmt_vac_logic" );

  // Set the visibility attributes for outer/mother volume(s)
  G4VisAttributes* visAttributes = GeoFactory::LoadVisualisation( table );
  n16SourceLog->SetVisAttributes( visAttributes );

  // Now position the entire source in the detector
  // Position by the centre of the gas volume
  // The z numbers come from the relative positioning of the gas volume (i.e. the volume filled by the n16 generator)
  // relative to the n16Source volume (i.e. the volume positioned in the macro).
  std::pair< G4ThreeVector, G4RotationMatrix* > translation = LoadTranslation( table );
  translation.first.setZ( translation.first.z() + ( 138.95 - 0.95 + 11.40 ) * mm );

  new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ), pmtVacLog, prefix + "_pmt_vac", pmtCaseLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, -53.15*mm ), pmtCaseLog, prefix + "_pmt_case", instrLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 89.75 ),   instrLog,  prefix + "_instr_phys", canLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, -49.1 ),   windowLog, prefix + "_window_phys", canLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, -138.95 ), gapLog,    prefix + "_gap_phys", canLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.0 ),     gasLog,    prefix + "_gas_phys", scintLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.95 ),    scintLog,  prefix + "_scint_phys", gapLog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, -11.4 ),   canLog,    prefix + "_can_phys", n16SourceLog, false, 0 );

  // The main (outer/mother) volume
  new G4PVPlacement( 0, translation.first, prefix + "_source_phys", n16SourceLog,  Detector::FindPhysicalVolume( motherName) , false, 0 );

}

} // namespace RAT
