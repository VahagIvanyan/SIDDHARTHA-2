#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Trap.hh"

#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
//#include "G4Transform3D.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fDefaultMaterial(0),fRingMaterial(0),fTitaniumMaterial(0),fGasTarget(0),fVacuumMat(0),fKaptonMaterial(0),fSiliconDetector(0),fPlasticDetector(0),fAluminumMaterial(0),fCarbonFiber(0),fMylarSctoch(0),fLeadShield(0),fCopperCon(0),fPhysiWorld(0),
 fDetectorMessenger(0)
{
  // default parameter values of the absorbers
  fNbOfAbsor = 2;
  fAbsorThickness[0] = 0*mm;        //dummy, for initialization   
  DefineMaterials();
for(G4int amount = 1; amount <= fNbOfAbsor; amount++)  { 
  fAbsorThickness[amount] = 0.75*cm; 
  fAbsorSizeYZ = 10*cm;
  //ComputeParameters();
  // materials
  SetAbsorMaterial(amount,"NE102");
}
  ComputeParameters();

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  man->FindOrBuildMaterial("G4_AIR");

	G4Element* H = man->FindOrBuildElement("H");
	G4Element* C = man->FindOrBuildElement("C");
	G4Element* Cu = man->FindOrBuildElement("Cu");
	G4Element* N = man->FindOrBuildElement("N");
	G4Element* O = man->FindOrBuildElement("O");
	G4Element* Fe = man->FindOrBuildElement("Fe");
	G4Element* Al = man->FindOrBuildElement("Al");
	G4Element* Ti = man->FindOrBuildElement("Ti");
	G4Element* Pb = man->FindOrBuildElement("Pb");
	G4Element* Si = man->FindOrBuildElement("Si");
  

	///////////////////////////////////////////////////////////////////////////////////
	///				  AIR					        ///
	///////////////////////////////////////////////////////////////////////////////////
	G4int ncomponents;
	G4double density,temperature,pressure,fractionmass; 
	
				density = 0.001225*g/cm3;
       				pressure = 101422.34057*pascal;
				temperature = 293*kelvin; 
				G4Material* Air = new G4Material("Air", density, ncomponents=2, kStateGas, temperature, pressure);
				Air->AddElement(N, fractionmass=79.*perCent);
				Air->AddElement(O, fractionmass=21.*perCent);
				G4NistManager *nist_man = G4NistManager::Instance();
				G4Material *AIR_mat = nist_man->FindOrBuildMaterial("Air");
	
//////	Graphite
  G4Material* Graphite = new G4Material("Graphite", 2.266 *g/cm3, 1);
  Graphite->AddElement(C,1);
  G4Material* Carbon = new G4Material("Carbon", 2.26*g/cm3, 1);
  Carbon->AddElement(C,1);	 			  
  G4Material* LiquidCarbon = new G4Material("LiquidCarbon", 1.2 *g/cm3, 1);
  LiquidCarbon->AddElement(C,1);
//////	Hydrogen
  G4Material* Hydrogen = new G4Material("Hydrogen", 0.071 *g/cm3, 1);
  Hydrogen->AddElement(H,1);
//KaptonPolymide
  G4Material* C22H10N2O5  =  new G4Material("KaptonPolymide", 1.43*g/cm3, 4);
  C22H10N2O5->AddElement(C,22);
  C22H10N2O5->AddElement(H,10);
  C22H10N2O5->AddElement(N,2);
  C22H10N2O5->AddElement(O,5);
  C22H10N2O5->GetIonisation()->SetMeanExcitationEnergy(79.6*eV);
//Copper
  G4Material* Copper  =  new G4Material("Cu", 8.96*g/cm3, 1);
  Copper->AddElement(Cu,1);
// Cryogenic Hydrogen
  density = 13.e-4*g/cm3;
  pressure    = 300000*pascal; 
  temperature = 30*kelvin; 
  G4Material* CryogHydro =  new G4Material("CryogHydro", 1., 1.00784*g/mole, density,
		                     kStateGas,temperature,pressure);	
//Water 
G4Material* H2O = 
  new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);  
// Example of Galactic Vacuum				
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* Galactic =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);
  G4Material* Galactic2 =   
  new G4Material("Galactic2", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);
//Vacuum inside of the Chamber
  density     = 1.e-16*g/cm3;   
  pressure    = 3.e-18*pascal;
  temperature = 20*kelvin;
  G4Material* Vacuum =   
  new G4Material("Vacuum", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);
//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//NE102 Plastic Scintillators
	G4Material* NE102 = new G4Material("NE102",1.032*g/cm3,2);
	NE102-> AddMaterial( Hydrogen,8.5*perCent);
 	NE102-> AddMaterial( Carbon,91.5*perCent);
//Silicon
	G4Material* Silicon = new G4Material("Silicon", 2.3290*g/cm3, 1);
	Silicon->AddElement(Si,1);
//Fe 
	G4Material* Iron = new G4Material("Iron", 7.850*g/cm3, 1);
	Iron->AddElement(Fe,1);
//Lead
	G4Material* Lead = new G4Material("Lead", 11.34*g/cm3, 1);
	Lead->AddElement(Pb,1);
//Aluminum
	G4Material* Aluminum = new G4Material("Al", 2.7*g/cm3, 1);
	Aluminum->AddElement(Al,1);
//Titanium
	G4Material* Titanium = new G4Material("Ti", 4.506*g/cm3, 1);
	Titanium->AddElement(Ti,1);
//CarbonFiber
	G4Material* CarbonFiber = new G4Material("CarbonFiber", 1.55*g/cm3, 1);
	CarbonFiber->AddElement(C,1);
//MYLAR
//C10H8O4
	G4Material* MYLAR = new G4Material("C10H8O4", 1.39*g/cm3, 3);
	MYLAR->AddElement(C,10);
	MYLAR->AddElement(H,8);
	MYLAR->AddElement(O,4);
//Low_Carbon_Stainless_Steel
	/*G4Material* Low_Carbon_Steel = new G4Material("Low_Carbon_Steel", 7.8499*g/cm3, 2);
  	Low_Carbon_Steel->AddMaterial(Iron,99.75*perCent);
  	Low_Carbon_Steel->AddMaterial(Carbon,0.25*perCent);*/
//Volume Materials
  fDefaultMaterial = AIR_mat;
  fRingMaterial = Galactic;
//Siddharta
  fLeadShield = Lead;
  fCopperCon = Copper;
  fGasTarget = CryogHydro;
fVacuumMat = Vacuum;
  fKaptonMaterial = C22H10N2O5; 
  fTitaniumMaterial = Titanium;
  fSiliconDetector = Silicon;
  fPlasticDetector = NE102;
  fAluminumMaterial = Aluminum;
  fCarbonFiber = CarbonFiber;
  fMylarSctoch = MYLAR;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeParameters()
{
  // Compute total thickness of absorbers
  fAbsorSizeX = 0.;
  for (G4int iAbs=1; iAbs<=fNbOfAbsor; iAbs++) {
    fAbsorSizeX += fAbsorThickness[iAbs];
  }

  fWorldSizeX = 2.0*m;
  fWorldSizeYZ = 2.0*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // complete the Calor parameters definition
  ComputeParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // World
  //
  G4Box* solidWorld =
    new G4Box("World",                                             //name
               1*m,1*m,1*m);       //size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,              //solid
                        fDefaultMaterial,        //material
                        "World");                //name

  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                      G4ThreeVector(),          //at (0,0,0)
                      logicWorld,               //logical volume
                      "World",                  //name
                       0,                       //mother volume
                       false,                   //no boolean operation
                       0);                      //copy number

/*===================================================================================================================*/
/*----------------------------------------------    SIDDHARTA Geometry    -------------------------------------------*/  
/*===================================================================================================================*/
///From the Bottom Side of SIDDHARTA 

 // rotation
G4RotationMatrix* rotation = new G4RotationMatrix();
rotation->rotateX(-90*deg);
G4RotationMatrix* rotation2 = new G4RotationMatrix();
     rotation2->rotateZ(0*deg);
	 rotation2->rotateX(-90*deg);
	 rotation2->rotateY(180*deg);
G4RotationMatrix* rotation3 = new G4RotationMatrix();
   	rotation3->rotateY(90*deg);
G4ThreeVector ConPlacing(0.,0.,0.);


/*============================================================================================*/
/*----------------------------------   Above the Pipe  ---------------------------------------*/
/*============================================================================================*/

///ConicalPlate
  G4Cons* solidCopperCon = new G4Cons("solidCopperCon", 4.95*cm, 5.2*cm, 3.85*cm, 4.1*cm, 3*cm, 0., twopi);
  G4LogicalVolume* LogicCopperCon = new G4LogicalVolume(solidCopperCon,fCopperCon,"CopperConLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,-32.705*cm,0.),LogicCopperCon,"CopperConPV",logicWorld,false,0);
  G4Cons* solidConPlate = new G4Cons("solidConPlate", 0., 5.4*cm, 0., 4.3*cm, 3.1*cm, 0., twopi);
  G4Cons* solidConRemovalPart = new G4Cons("solidConRemovalPart", 0., 5.2*cm, 0., 4.1*cm, 3.12*cm, 0., twopi);
  G4Box* solidBoxShape = new G4Box("BoxShape",50.*cm,3*cm,29*cm);  
  G4UnionSolid* UnionSolid1 = new G4UnionSolid("BoxShape+solidConPlate",solidBoxShape,solidConPlate,rotation,ConPlacing); 
  G4SubtractionSolid* SubtractionSolid1 = new G4SubtractionSolid("BoxShape+solidConPlate-solidConRemovalPart",
 		      		UnionSolid1,solidConRemovalPart,rotation,ConPlacing); 
  G4LogicalVolume* LogicSubtractionSolid1 = new G4LogicalVolume(SubtractionSolid1,fLeadShield,"SubtractionSolidLV");
  new G4PVPlacement(rotation3,G4ThreeVector(0.,-32.705*cm,0.),LogicSubtractionSolid1,"SubtractionSolidPV",logicWorld,false,0);
/*=============================================================================================*/
/*------------------------------------------	Pipe	---------------------------------------*/
/*=============================================================================================*/
  G4Tubs* solidPipe = new G4Tubs("solidPipe",58.3*mm,58.3*mm+0.650*mm,175*mm,0.,twopi);
  G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe,fAluminumMaterial,"PipeLV");
  new G4PVPlacement(0,G4ThreeVector(0.,-50.*cm,0.), logicPipe,"PipePV",logicWorld,false,0);
//CarbonFiber
  G4Tubs* solidCarbonFiber = new G4Tubs("solidCarbonFiber",59.3*mm,62.1*mm,175*mm,0.,twopi);
  G4LogicalVolume* logicCarbonFiber = new G4LogicalVolume(solidCarbonFiber,fCarbonFiber,"CarbonFiberLV");
  new G4PVPlacement(0,G4ThreeVector(0.,-50.*cm,0.), logicCarbonFiber,"solidCarbonFiberPV",logicWorld,false,0);

/*============================================================================================*/
/*-------------------------------------  Vacuum Chamber  -------------------------------------*/
/*============================================================================================*/
//Aluminum Chamber
  G4Tubs* VacuumChamberDown = new G4Tubs("VacuumChamberDown",20.2*cm,21.2*cm,17.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChamberDown = new G4LogicalVolume(VacuumChamberDown,fAluminumMaterial,"VacuumChamberDownLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,-10.205*cm,0.),logicVacuumChamberDown,"VacuumChamberDownPV",logicWorld,false,0);
  G4Tubs* VacuumChamberUp = new G4Tubs("VacuumChamberUp",35.5*cm,36.3*cm,7.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChamberUp = new G4LogicalVolume(VacuumChamberUp,fAluminumMaterial,"VacuumChamberUpLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,14.795*cm,0.),logicVacuumChamberUp,"VacuumChamberUpPV",logicWorld,false,0);
  G4Tubs* VacuumInsideChamberUp = new G4Tubs("VacuumInsideChamberUp",0*cm,35.5*cm,7.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumInsideChamberUp = new G4LogicalVolume(VacuumInsideChamberUp,fVacuumMat,"VacuumInsideChamberUpLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,14.795*cm,0.),logicVacuumInsideChamberUp,"VacuumInsideChamberUpPV",logicWorld,false,0);
  G4Tubs* solidVacuumChambUp = new G4Tubs("VacuumChambUp",20.2*cm,35.5*cm,0.4*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChambUp = new G4LogicalVolume(solidVacuumChambUp,fAluminumMaterial,"VacuumChambUpLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,6.895*cm,0.),logicVacuumChambUp,"VacuumChambUpPV",logicWorld,false,0);
  G4Tubs* solidVacuumChambUp2 = new G4Tubs("VacuumChambUp2",15*cm,37.5*cm,2*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChambUp2 = new G4LogicalVolume(solidVacuumChambUp2,fAluminumMaterial,"VacuumChambUp2LV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,24.295*cm,0.),logicVacuumChambUp2,"VacuumChambUp2PV",logicWorld,false,0);

  G4Tubs* solidVacuumChambUp3 = new G4Tubs("VacuumChambUp3",13*cm,15*cm,23.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChambUp3 = new G4LogicalVolume(solidVacuumChambUp3,fAluminumMaterial,"VacuumChambUp3LV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,49.795*cm,0.),logicVacuumChambUp3,"VacuumChambUp3PV",logicWorld,false,0);

  G4Tubs* VacuumInsideChambUp3 = new G4Tubs("VacuumInsideChambUp3",0*cm,13*cm,24.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumInsideChambUp3 = new G4LogicalVolume(VacuumInsideChambUp3,fVacuumMat,"VacuumInsideChambUp3LV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,47.795*cm,0.),logicVacuumInsideChambUp3,"VacuumInsideChambUp3PV",logicWorld,false,0);

  G4Tubs* solidVacuumChambUp4 = new G4Tubs("VacuumChambUp4",0*cm,13*cm,1*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChambUp4 = new G4LogicalVolume(solidVacuumChambUp4,fAluminumMaterial,"VacuumChambUp4LV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,74.295*cm,0.),logicVacuumChambUp4,"VacuumChambUp4PV",logicWorld,false,0);
//Vacuum Volume
  G4Tubs* VacuumVolume = new G4Tubs("VacuumVolume",0*cm,20.2*cm,17.5*cm,0.,twopi);
  G4LogicalVolume* logicVacuumVolume = new G4LogicalVolume(VacuumVolume,fVacuumMat,"VacuumVolumeLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,-10.205*cm,0.),logicVacuumVolume,"VacuumVolumePV",logicWorld,false,0);
//Aluminum Bottom Ring
  G4Tubs* solidVacuumChamberBottomRing = new G4Tubs("ChamberBottomRing",5.8*cm,21.2*cm,1*cm,0.,twopi);
  G4LogicalVolume* logicVacuumChamberBottomRing = new G4LogicalVolume(solidVacuumChamberBottomRing,fAluminumMaterial,"ChamberBottomRingLV");
  new G4PVPlacement(rotation,G4ThreeVector(0.,-28.705*cm,0.),logicVacuumChamberBottomRing,"ChamberBottomRingPV",logicWorld,false,0);
/*======================================================================================*/
/*--------------------------      Gaseous Target Entrance         ----------------------*/	
/*======================================================================================*/
//Al Outer Ring(Under the Target) 
  G4Tubs* AlOuterRingBottom = new G4Tubs("AlOuterRingBottom",7.5*cm,8.75*cm,1*cm,0.,twopi);
  G4LogicalVolume* logicAlOuterRingBottom = new G4LogicalVolume(AlOuterRingBottom,fAluminumMaterial,"AlOuterRingBottomLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,16.5*cm),logicAlOuterRingBottom,"AlOuterRingBottomPV",logicVacuumVolume,false,0);
//KaptonWindow
  G4Tubs* solidKaptonWindow = new G4Tubs("KaptonWindow",0.,5.8*cm,0.0625*mm,0.,twopi);
  G4LogicalVolume* logicKaptonWindow = new G4LogicalVolume(solidKaptonWindow,fKaptonMaterial,"KaptonWindowLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,15.5*cm),logicKaptonWindow,"KaptonWindowPV",logicVacuumVolume,false,0);
// Al Inner Ring(Under the Target)
  G4Tubs* AlInnerRingBottom = new G4Tubs("AlInnerRingBottom",5.8*cm,7.5*cm,0.5*cm,0.,twopi);
  G4LogicalVolume* logicAlInnerRingBottom = new G4LogicalVolume(AlInnerRingBottom,fAluminumMaterial,"AlInnerRingBottomLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,15*cm),logicAlInnerRingBottom,"AlInnerRingBottomPV",logicVacuumVolume,false,0);
/*====================================*/
/*-------    Volume for SDD    -------*/	
/*====================================*/
// SDD Upper Ring
  G4Tubs* solidSDDupRing = new G4Tubs("solidSDDupRing",7.275*cm,7.5*cm,2.25*cm,0.,twopi);
  G4LogicalVolume* logicSDDupRing = new G4LogicalVolume(solidSDDupRing,fRingMaterial,"SDDupRingLV");
for (G4int iSDD = 1; iSDD <= 24 ; iSDD++) {
//SDD
  G4Box* solidSDD = new G4Box("solidSDD",0.25*mm,1.7*cm,0.9*cm); 
  G4LogicalVolume* logicSDD = new G4LogicalVolume(solidSDD,fSiliconDetector,"SDDLV");
    G4double phi = iSDD*(twopi/24);
	G4RotationMatrix SDDRot = G4RotationMatrix();
       	SDDRot.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),std::sin(phi),0.);     
    G4ThreeVector SDDposision = 7.325*cm*uz;
    G4Transform3D transform = G4Transform3D(SDDRot,SDDposision);
  new G4PVPlacement(transform,logicSDD,"SDDPV",logicSDDupRing,false,iSDD);
  G4VisAttributes* VisSDD = 0;
  VisSDD = new G4VisAttributes(G4Colour(1.,0.,0.));
  VisSDD->SetVisibility(true);
  logicSDD->SetVisAttributes(VisSDD);
}
  new G4PVPlacement(0,G4ThreeVector(0.,0.,10.95*cm),logicSDDupRing,"SDDupRingPV",logicVacuumVolume,false,0);
// SDD down ring
  G4Tubs* solidSDDdownRing = new G4Tubs("solidSDDdownRing",7.275*cm,7.4*cm,2.25*cm,0.,twopi);
  G4LogicalVolume* logicSDDdownRing = new G4LogicalVolume(solidSDDdownRing,fRingMaterial,"SDDdownRingLV");
for (G4int jSDD = 1; jSDD <= 24 ; jSDD++) {
//SDD
  G4Box* soliddownSDD = new G4Box("soliddownSDD",0.25*mm,1.7*cm,0.9*cm); 
  G4LogicalVolume* logicdownSDD = new G4LogicalVolume(soliddownSDD,fSiliconDetector,"SDDdownLV");
    G4double phi = jSDD*(twopi/24);
	G4RotationMatrix SDDdownRot = G4RotationMatrix();
        //SDDdownRot.rotateX(-90*deg);
       	SDDdownRot.rotateZ(phi);
    G4ThreeVector uzdown = G4ThreeVector(std::cos(phi),std::sin(phi),0.);     
    G4ThreeVector SDDdownposision = 7.325*cm*uzdown;
    G4Transform3D transform1 = G4Transform3D(SDDdownRot,SDDdownposision);
  new G4PVPlacement(transform1,logicdownSDD,"SDDdownPV",logicSDDdownRing,false,jSDD);
  G4VisAttributes* Vis_SDD = 0;
  Vis_SDD = new G4VisAttributes(G4Colour(1., 0., 0.));
  Vis_SDD->SetVisibility(true);
  logicdownSDD->SetVisAttributes(Vis_SDD);
 }
  new G4PVPlacement(0,G4ThreeVector(0.,0.,7.55*cm),logicSDDdownRing,"SDDdownRingPV",logicVacuumVolume,false,0);
//Target Volume(Kapton)
  G4Tubs* solidTargetVolume = new G4Tubs("solidTargetVolume",7.2*cm,7.215*cm,6.25*cm,0.,twopi);
  G4LogicalVolume* logicTargetVolume = new G4LogicalVolume(solidTargetVolume,fKaptonMaterial,"TargetVolumeLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,9.25*cm),logicTargetVolume,"solidTargetVolumePV",logicVacuumVolume,false,0);
//Target
  G4Tubs* GasTarget = new G4Tubs("GasTarget",0.,7.2*cm,6.25*cm,0.,twopi);
  G4LogicalVolume* logicTarget = new G4LogicalVolume(GasTarget,fGasTarget,"TargetLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,9.25*cm),logicTarget,"GasTargetPV",logicVacuumVolume,false,0);
//Al Inner Rings(Above the Target)
  G4Tubs* AlInnerRingAbove = new G4Tubs("AlInnerRingAbove",5.5*cm,7.25*cm,0.4*cm,0.,twopi);
  G4LogicalVolume* logicAlInnerRingAbove = new G4LogicalVolume(AlInnerRingAbove,fAluminumMaterial,"AlInnerRingAboveLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,3.395*cm),logicAlInnerRingAbove,"AlInnerRingAbovePV",logicVacuumVolume,false,0);
//Titanium Exit Window
  G4Tubs* solidTitaniumWindow = new G4Tubs("TitaniumWindow",0.,5.5*cm,0.05*mm,0.,twopi);
  G4LogicalVolume* logicTitaniumWindow = new G4LogicalVolume(solidTitaniumWindow,fTitaniumMaterial,"TitaniumWindowLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,4*cm),logicTitaniumWindow,"TitaniumWindowPV",logicVacuumVolume,false,0);
//Al Outer Rings(Above the Target)
  G4Tubs* AlOuterRingAbove = new G4Tubs("AlOuterRingAbove",5.5*cm,8.75*cm,0.6*cm,0.,twopi);
  G4LogicalVolume* logicAlOuterRingAbove = new G4LogicalVolume(AlOuterRingAbove,fAluminumMaterial,"AlOuterRingAboveLV");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.395*cm),logicAlOuterRingAbove,"AlOuterRingAbovePV",logicVacuumVolume,false,0);

/*============================================================================================================================*/
/*-----------------------------------------	Detector(for Kaons)	------------------------------------------------------*/
/*============================================================================================================================*/
  fXfront[0] = -0.5*fAbsorSizeX;
  G4double u_foil,u,pos=-38*cm;
  //
for(G4int k_foil = 1; k_foil<=2;k_foil++)
{
//Wrapping Scotch (Mylar) for Kaon Detectors
G4Box* solidMylarVolume = new G4Box("solidMylarVolume",5*cm+0.006*mm,4.5*cm+0.003*mm+0.25*mm,0.75*mm+0.006*mm+0.25*mm);
G4Box* solidBoxRemovalMylarFoil = new G4Box("BoxRemovalFoil",5*cm+0.003*mm,4.51*cm,0.75*mm+0.006*mm);
if (k_foil == 1)
 		u_foil = pos;
	 else
 		u_foil = pos - 24*cm;
G4ThreeVector pos_foil = G4ThreeVector(0.,u_foil,0.);
G4SubtractionSolid* SubtractionSolid3 = new G4SubtractionSolid("solidMylarVolume-BoxRemovalFoil",
 		      		solidMylarVolume,solidBoxRemovalMylarFoil,rotation,pos_foil); 
G4LogicalVolume* LogicSubtractionSolid3 = new G4LogicalVolume(SubtractionSolid3,fMylarSctoch,"SubtractionSolid3");
//AlFoil for Kaon Detectors
G4Box* solidAlFoil = new G4Box("AlFoil",5*cm+0.003*mm,4.5*cm,0.75*mm+0.003*mm);
G4Box* solidBoxRemovalAlFoil = new G4Box("solidBoxRemovalAlFoil",5*cm,4.51*cm,0.75*mm);
G4SubtractionSolid* SubtractionSolid4 = new G4SubtractionSolid("solidAlFoil-solidBoxRemovalAlFoil",
 		      		solidAlFoil,solidBoxRemovalAlFoil,rotation,pos_foil); 
G4LogicalVolume* LogicSubtractionSolid4 = new G4LogicalVolume(SubtractionSolid4,fAluminumMaterial,"solidAlFoil-solidBoxRemovalAlFoil");
new G4PVPlacement(rotation,pos_foil,LogicSubtractionSolid4, "SubtractionSolid4", logicWorld, false,k_foil);
new G4PVPlacement(rotation,pos_foil,LogicSubtractionSolid3, "SubtractionSolid3", logicWorld, false,k_foil);

G4VisAttributes* VisKaonDetect = 0;
  VisKaonDetect = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  VisKaonDetect->SetVisibility(true);
  LogicSubtractionSolid3->SetVisAttributes(VisKaonDetect);
  VisKaonDetect = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  VisKaonDetect->SetVisibility(true);
  LogicSubtractionSolid3->SetVisAttributes(VisKaonDetect);
  VisKaonDetect = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisKaonDetect->SetVisibility(true);
  LogicSubtractionSolid4->SetVisAttributes(VisKaonDetect);
}
  for (G4int k=1; k<=fNbOfAbsor; k++) {

    G4Material* material = fAbsorMaterial[k];
    G4String matname = material->GetName();     
	G4Box* solidAbsor = new G4Box(matname,5*cm,5*cm,0.75*mm);
	G4LogicalVolume* logicAbsor = new G4LogicalVolume(solidAbsor,fPlasticDetector,matname);
	 if (k == 1)
 		u = pos;
	 else
 		u = pos - 24*cm;
	 G4ThreeVector position = G4ThreeVector(0.,u,0.);
      	new G4PVPlacement(rotation,        //no rotation
                        position,        //position
                        logicAbsor,      //logical volume        
                        matname,         //name
                        logicWorld,  	 //mother
                        false,                 //no boulean operat
                        k);                    //copy number
G4VisAttributes* VisKaonDetector= 0;
  VisKaonDetector = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  VisKaonDetector->SetVisibility(true);
  logicAbsor->SetVisAttributes(VisKaonDetector);
  }
//Visualisation
  G4VisAttributes* VisAtt = 0;
  VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  VisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
  VisAtt->SetVisibility(true);
  logicPipe->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  VisAtt->SetVisibility(true);
  logicCarbonFiber->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.396,0.263,0.129));
  VisAtt->SetVisibility(true);
  LogicCopperCon->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  VisAtt->SetVisibility(true);
  LogicSubtractionSolid1->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
  VisAtt->SetVisibility(true);
  logicTarget->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.8515625, 0.64453125, 0.125));
  VisAtt->SetVisibility(true);
  logicTargetVolume->SetVisAttributes(VisAtt);
  VisAtt = new G4VisAttributes(G4Colour(0.8515625, 0.64453125, 0.125));
  VisAtt->SetVisibility(true);
  logicKaptonWindow->SetVisAttributes(VisAtt);
  logicVacuumVolume ->SetVisAttributes(G4VisAttributes::Invisible);
  logicSDDupRing ->SetVisAttributes(G4VisAttributes::Invisible);
  logicSDDdownRing ->SetVisAttributes(G4VisAttributes::Invisible);
  logicVacuumInsideChamberUp ->SetVisAttributes(G4VisAttributes::Invisible);
  logicVacuumInsideChambUp3 ->SetVisAttributes(G4VisAttributes::Invisible);

  VisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
  VisAtt->SetVisibility(true);
  logicVacuumChamberDown->SetVisAttributes(VisAtt);
  logicVacuumChamberUp->SetVisAttributes(VisAtt);
  logicVacuumChambUp->SetVisAttributes(VisAtt);
  logicVacuumChambUp2->SetVisAttributes(VisAtt);
  logicVacuumChambUp3->SetVisAttributes(VisAtt);
  logicVacuumChambUp4->SetVisAttributes(VisAtt);
  logicVacuumChamberBottomRing->SetVisAttributes(VisAtt);
  PrintParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Absorber is " << fNbOfAbsor << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << kMaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int iabs,const G4String& material)
{
  // search the material by its name
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorMaterial: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
      fAbsorMaterial[iabs] = pttoMaterial;
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int iabs,G4double val)
{
  // change Absorber thickness
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorThickness: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[iabs] = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorSizeYZ = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
