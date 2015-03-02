//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: pixelTBDetectorConstruction.cc,v 1.22 2010-01-22 11:57:03 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "pixelTBDetectorConstruction.hh"
#include "pixelTBDetectorMessenger.hh"
#include "pixelTBChamberParameterisation.hh"
#include "pixelTBMagneticField.hh"
#include "pixelTBTrackerSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
pixelTBDetectorConstruction::pixelTBDetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0),
 solidTracker(0),
 solidChamber(0), 
 TargetMater(0), ChamberMater(0),chamberParam(0),
 stepLimit(0), fpMagField(0),
 fWorldLength(0.),  fTargetLength(0.), fTrackerLength(0.),
 NbOfChambers(0) ,  ChamberWidth(0.),  ChamberSpacing(0.)
{
  fpMagField = new pixelTBMagneticField();
  detectorMessenger = new pixelTBDetectorMessenger(this);
   for(size_t ii=0; ii<16; ii++){
     logicTracker[ii]=0;
     physiTracker[ii]=0;
     logicChamber[ii]=0;
     physiChamber[ii]=0;
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
pixelTBDetectorConstruction::~pixelTBDetectorConstruction()
{
  delete fpMagField;
  delete stepLimit;
  delete chamberParam;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* pixelTBDetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);

  G4double a, z;
  G4double density;
  G4int nel;
  G4int natoms;

  //Elements
  G4Element* H  = manager->FindOrBuildElement(1);
  G4Element* C  = manager->FindOrBuildElement(6);
  G4Element* N  = manager->FindOrBuildElement(7);
  G4Element* O  = manager->FindOrBuildElement(8);
  G4Element* Si = manager->FindOrBuildElement(14);
  G4Element* Cl  = manager->FindOrBuildElement(17);

  //Air
  G4Material* Air = new G4Material("Air", density= 1.29*CLHEP::mg/CLHEP::cm3, nel=2);
  Air->AddElement(N, 70*CLHEP::perCent);
  Air->AddElement(O, 30*CLHEP::perCent);

  //Air
  G4Material* Vacuum = new G4Material("InterGalactic", 1, 1.008*CLHEP::g/CLHEP::mole, 1.e-25*CLHEP::g/CLHEP::cm3, kStateGas, 2.73*CLHEP::kelvin, 3.e-18*CLHEP::pascal);


  //Silicon
  G4Material* Silicon = 
  new G4Material("Silicon", z=14., a=28.0855*CLHEP::g/CLHEP::mole, density=2.330*CLHEP::g/CLHEP::cm3);

  //G10
  G4Material* G10 = 
    new G4Material("G10", density= 1.700*CLHEP::g/CLHEP::cm3, nel=4);
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);

  //Aluminum
  G4Material* Alu = 
  new G4Material("Alu", z=13., a= 26.98*CLHEP::g/CLHEP::mole, density= 2.700*CLHEP::g/CLHEP::cm3);

  //Plastic -- polivenylclorid
  G4Material* PVC = 
    new G4Material("PVC", density = 1.40*CLHEP::g/CLHEP::cm3, nel=3);
  PVC->AddElement(Cl, natoms=1);
  PVC->AddElement(C , natoms=2);
  PVC->AddElement(H , natoms=3);

  //Insulation-- polystyrol
  G4Material* insulation = 
    new G4Material("insulation", density = 0.03*CLHEP::g/CLHEP::cm3, nel=2);
  insulation->AddElement(C , natoms=8);
  insulation->AddElement(H , natoms=8);

  G4Material* mylar = 
    new G4Material("mylar", density = 1.38*CLHEP::g/CLHEP::cm3, nel=3);
  mylar->AddElement(O, natoms=4); 
  mylar->AddElement(C, natoms=10); 
  mylar->AddElement(H, natoms=8); 



  // Print all the materials defined.
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;



//--------- Sizes of the principal geometrical components (solids)  ---------
  
//  NbOfChambers = 52*80;
  NbOfChambers = 52*80*9;
  //ChamberSpacing = 12.91*CLHEP::mm;
  ChamberSpacing = 16.0*CLHEP::mm;
  ChamberWidth = 0.285*CLHEP::mm;
  ChamberSize = 1.0*CLHEP::cm;

  ROCSize     = 1.*CLHEP::cm;
  ROCWidth     = 0.200*CLHEP::mm;
  BoardWidth   = 1.2*CLHEP::mm;

  G4double ChamberXShalfwidth = 0.05*CLHEP::mm;


  fTrackerLength = (8+8+3)*ChamberSpacing; // Full length of Tracker
  fTargetLength  = 50*CLHEP::cm;               // Full length of Target
  
  TargetMater  = Vacuum;
  ChamberMater = Silicon;
  BoardMater   = G10;
  
  fWorldLength= 2.0*(fTargetLength+fTrackerLength);
   
  G4double targetSize  = 0.5*fTargetLength;    // Half length of the Target  
  G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker

  G4double BoardSize = 4.2*CLHEP::cm;

  G4double GeneralTranslationZ = 5.0*CLHEP::cm;

  G4double BoxX = 50.0*CLHEP::cm;
  G4double BoxY = 33.0*CLHEP::cm;
  G4double BoxZ = 51.0*CLHEP::cm;

  G4double WallThickness = 1.5*CLHEP::cm; 
  G4double InsulationThickness = 3.0*CLHEP::cm; 


  G4double ScintInnerRad = 0.9*CLHEP::cm;
  G4double ScintOuterRad = 1.0*CLHEP::cm;
  G4double ScintLength   = 10.0*CLHEP::cm;
  G4double ScintStartAngle  = 0;
  G4double ScintSegmentAngle  = CLHEP::pi2;

  G4double ScintPosition = GeneralTranslationZ + WallThickness + InsulationThickness + 4.0*CLHEP::cm + 0.5*ScintOuterRad;


  G4double FloorX = 24.0*CLHEP::cm;
  G4double FloorY = 2.0*CLHEP::cm;
  G4double FloorZ = 30.0*CLHEP::cm;

  G4double FloorPositionY = -0.5*BoardSize - 0.5*FloorY;
  G4double FloorPositionZ = ScintPosition + ScintOuterRad + 0.5*FloorZ;

  G4double tiltedTelescopePosition = GeneralTranslationZ + WallThickness + InsulationThickness + 6.0*CLHEP::cm;

  G4double WallOutX = BoxX;
  G4double WallInX = WallOutX - 2*WallThickness;
  G4double WallOutY = BoxY;
  G4double WallInY = WallOutY - 2*WallThickness;
  G4double WallOutZ = BoxZ;
  G4double WallInZ = WallOutZ - 2*WallThickness;
  
  G4double InsulationOutX = BoxX - 2*WallThickness;
  G4double InsulationInX = InsulationOutX - 2*InsulationThickness;
  G4double InsulationOutY = BoxY - 2*WallThickness;
  G4double InsulationInY = InsulationOutY - 2*InsulationThickness;
  G4double InsulationOutZ = BoxZ - 2*WallThickness;
  G4double InsulationInZ = InsulationOutZ - 2*InsulationThickness;
  

  G4double PeltierX = BoxX - 2* WallThickness -2*InsulationThickness;
  G4double PeltierY = 2.0*CLHEP::cm;
  G4double PeltierZ = BoxZ - 2* WallThickness -2*InsulationThickness;

  G4double BoxCenterXPosition = 0.0*CLHEP::cm;
  G4double BoxCenterYPosition = +0.5*BoxY - WallThickness - InsulationThickness - PeltierY - FloorY - 0.5*BoardSize;
  G4double BoxCenterZPosition = 0.5*BoxZ + GeneralTranslationZ;

  G4double capOutSize = 2.0*CLHEP::cm; 
  G4double capOutHight = 0.5*CLHEP::cm; 

  G4double capInSize = 1.95*CLHEP::cm; 
  G4double capInHight = 0.45*CLHEP::cm; 





//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  


  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
 
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/CLHEP::mm
         << " CLHEP::mm" << G4endl;

  solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother volume
                                 false,           // no boolean operations
                                 0);              // copy number




		 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetSize+trackerSize));
   
  solidTarget = new G4Box("target",targetSize,targetSize,targetSize);
  logicTarget = new G4LogicalVolume(solidTarget,TargetMater,"Target",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume				  
				  "Target",        // its name
				  logicWorld,      // its mother volume
				  false,           // no boolean operations
				  0);              // copy number 

  G4cout << "Target is " << fTargetLength/CLHEP::cm << " cm of " 
         << TargetMater->GetName() << G4endl;

  


  //------------------------------ 
  // Tracker (Pixels mother volumes)
  //------------------------------
  
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
  G4RotationMatrix rotationTracker = G4RotationMatrix(0.,0.,0.);
  G4double angle = 0.0;
  angle = -20.0;
  G4RotationMatrix *mirror_rot = new G4RotationMatrix(G4ThreeVector(-1.5, 1.0, 0.0),angle*CLHEP::deg);
  
  solidTracker = new G4Box("tracker",0.5*ChamberSize, 0.5*ChamberSize, 0.5*ChamberWidth);

  G4String trackerName[16];
  trackerName[0] = "Tracker-0";
  trackerName[1] = "Tracker-1";
  trackerName[2] = "Tracker-2";
  trackerName[3] = "Tracker-3";
  trackerName[4] = "Tracker-4";
  trackerName[5] = "Tracker-5";
  trackerName[6] = "Tracker-6";
  trackerName[7] = "Tracker-7";
  trackerName[8] = "Tracker-8";
  trackerName[9] = "Tracker-9";
  trackerName[10] = "Tracker-10";
  trackerName[11] = "Tracker-11";
  trackerName[12] = "Tracker-12";
  trackerName[13] = "Tracker-13";
  trackerName[14] = "Tracker-14";
  trackerName[15] = "Tracker-15";

  G4ThreeVector posT[16];

  for(size_t ii=0; ii<16; ii++) posT[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, tiltedTelescopePosition + ChamberSpacing*G4double(ii+1));
  //for(size_t ii=4; ii<8; ii++) posT[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, ChamberSpacing*G4double(ii+2));

  for(size_t ii=0; ii<8; ii++){

    logicTracker[ii] = new G4LogicalVolume(solidTracker , Air, "Tracker",0,0,0);
    physiTracker[ii] = new G4PVPlacement(mirror_rot,//0,              // rotated
					 posT[ii],           // at (x,y,z)
					 logicTracker[ii],   // its logical volume				  
					 trackerName[ii],    // its name
					 logicWorld,         // its mother  volume
					 false,              // no boolean operations
					 ii);                // copy number 
  }

  for(size_t ii=8; ii<16; ii++){

    logicTracker[ii] = new G4LogicalVolume(solidTracker , Air, "Tracker",0,0,0);
    physiTracker[ii] = new G4PVPlacement(0,              // no rotation
					 posT[ii],           // at (x,y,z)
					 logicTracker[ii],   // its logical volume				  
					 trackerName[ii],    // its name
					 logicWorld,         // its mother  volume
					 false,              // no boolean operations
					 ii);                // copy number 
  }




  //------------------------------ 
  // Tracker segments
  //------------------------------

  
  solidChamber = new G4Box("chamber", 0.5*ChamberSize, 0.5*ChamberSize, 0.5*ChamberWidth);

  for(size_t ii=0; ii<16; ii++) logicChamber[ii] = new G4LogicalVolume(solidChamber,ChamberMater,"Chamber",0,0,0);
  
  G4double firstPosition = -trackerSize + 0.5*ChamberWidth;
  G4double firstLength = ChamberXShalfwidth;
  //G4double lastLength  = ChamberXShalfwidth;
  fCSZonePercentageX = 0.2;
  fCSZonePercentageY = 0.2;
			   
  // dummy value : kZAxis -- modified by parameterised volume
  //
  chamberParam = new pixelTBChamberParameterisation(  
		         NbOfChambers,          // NoChambers 
			 firstPosition,         // Z of center of first 
			 ChamberSpacing,        // Z spacing of centers
			 ChamberWidth,          // Width Chamber 
			 firstLength,           // lengthInitial 
			 0,                     // lengthFinal
			 fCSZonePercentageX,
			 fCSZonePercentageY);          


  for(size_t ii=0; ii<16; ii++){
    physiChamber[ii] = new G4PVParameterised(
			 "Chamber",       // their name
			 logicChamber[ii],    // their logical volume
			 logicTracker[ii],    // Mother logical volume
			 kZAxis,          // Are placed along this axis 
			 NbOfChambers,    // Number of chambers
			 chamberParam);   // The parametrisation
  }





  G4cout << "There are " << NbOfChambers << " chambers in the tracker region. "
         << "The chambers are " << ChamberWidth/CLHEP::mm << " mm of " 
         << ChamberMater->GetName() << "\n The distance between chamber is "
	 << ChamberSpacing/CLHEP::cm << " cm" << G4endl;



  //------------------------------ 
  // Mother volumes for ROC and board
  //------------------------------



  solidRandB = new G4Box("RandB",0.5*BoardSize, 0.5*BoardSize, 0.5*(BoardWidth+ROCWidth));

  G4String RandBName[16];
  RandBName[0] = "RandB-0";
  RandBName[1] = "RandB-1";
  RandBName[2] = "RandB-2";
  RandBName[3] = "RandB-3";
  RandBName[4] = "RandB-4";
  RandBName[5] = "RandB-5";
  RandBName[6] = "RandB-6";
  RandBName[7] = "RandB-7";
  RandBName[8] = "RandB-8";
  RandBName[9] = "RandB-9";
  RandBName[10] = "RandB-10";
  RandBName[11] = "RandB-11";
  RandBName[12] = "RandB-12";
  RandBName[13] = "RandB-13";
  RandBName[14] = "RandB-14";
  RandBName[15] = "RandB-15";

  G4ThreeVector posT2[16];

  //original  for(size_t ii=0; ii<16; ii++) posT2[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, tiltedTelescopePosition + ChamberSpacing*G4double(ii+1)+(0.5*(BoardWidth+ROCWidth+ChamberWidth)));
  for(size_t ii=0; ii<16; ii++) posT2[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, tiltedTelescopePosition + ChamberSpacing*G4double(ii+1)-(0.5*(BoardWidth+ROCWidth+ChamberWidth)));
  //for(size_t ii=4; ii<8; ii++) posT2[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, ChamberSpacing*G4double(ii+2)+(0.5*(BoardWidth+ROCWidth+ChamberWidth)));

  for(size_t ii=0; ii<8; ii++){

    logicRandB[ii] = new G4LogicalVolume(solidRandB , Air, "RandB",0,0,0);
    physiRandB[ii] = new G4PVPlacement(mirror_rot,              // rotated
					 posT2[ii],           // at (x,y,z)
					 logicRandB[ii],   // its logical volume				  
					 RandBName[ii],    // its name
					 logicWorld,         // its mother  volume
					 false,              // no boolean operations
					 ii);                // copy number 
  }

  for(size_t ii=8; ii<16; ii++){

    logicRandB[ii] = new G4LogicalVolume(solidRandB , Air, "RandB",0,0,0);
    physiRandB[ii] = new G4PVPlacement(0,              // no rotation
				       posT2[ii],           // at (x,y,z)
				       logicRandB[ii],   // its logical volume				  
				       RandBName[ii],    // its name
				       logicWorld,         // its mother  volume
				       false,              // no boolean operations
				       ii);                // copy number 
  }



  //------------------------------ 
  // ROC
  //------------------------------



  solidROC = new G4Box("roc", 0.5*ROCSize, 0.5*ROCSize, 0.5*ROCWidth);

  for(size_t ii=0; ii<16; ii++) logicROC[ii] = new G4LogicalVolume(solidROC,ChamberMater,"ROC",0,0,0);

  G4String ROCName[16];
  ROCName[0] = "ROC-0";
  ROCName[1] = "ROC-1";
  ROCName[2] = "ROC-2";
  ROCName[3] = "ROC-3";
  ROCName[4] = "ROC-4";
  ROCName[5] = "ROC-5";
  ROCName[6] = "ROC-6";
  ROCName[7] = "ROC-7";
  ROCName[8] = "ROC-8";
  ROCName[9] = "ROC-9";
  ROCName[10] = "ROC-10";
  ROCName[11] = "ROC-11";
  ROCName[12] = "ROC-12";
  ROCName[13] = "ROC-13";
  ROCName[14] = "ROC-14";
  ROCName[15] = "ROC-15";


  for(size_t ii=0; ii<16; ii++){
    physiROC[ii] = new G4PVPlacement(0,//mirror_rot,//0,              // no rotation
				     //original					 G4ThreeVector(0,0,-0.5*BoardWidth), // at (x,y,z)
					 G4ThreeVector(0,0,0.5*BoardWidth), // at (x,y,z)
					 logicROC[ii],    // its logical volume				  
					 ROCName[ii],       // its name
					 logicRandB[ii],      // its mother  volume
					 false,           // no boolean operations
					 ii);              // copy number 
  }



  //------------------------------ 
  // G10 Board
  //------------------------------


  solidBoard = new G4Box("board", 0.5*BoardSize, 0.5*BoardSize, 0.5*BoardWidth);

  for(size_t ii=0; ii<16; ii++) logicBoard[ii] = new G4LogicalVolume(solidBoard,BoardMater,"Board",0,0,0);

  G4String boardName[16];
  boardName[0] = "Board-0";
  boardName[1] = "Board-1";
  boardName[2] = "Board-2";
  boardName[3] = "Board-3";
  boardName[4] = "Board-4";
  boardName[5] = "Board-5";
  boardName[6] = "Board-6";
  boardName[7] = "Board-7";
  boardName[8] = "Board-8";
  boardName[9] = "Board-9";
  boardName[10] = "Board-10";
  boardName[11] = "Board-11";
  boardName[12] = "Board-12";
  boardName[13] = "Board-13";
  boardName[14] = "Board-14";
  boardName[15] = "Board-15";

  for(size_t ii=0; ii<16; ii++){
    physiBoard[ii] = new G4PVPlacement(0,//mirror_rot,//0,              // no rotation
				       //original				       G4ThreeVector(0,0,0.5*ROCWidth), // at (x,y,z)
				       G4ThreeVector(0,0,-0.5*ROCWidth), // at (x,y,z)
				       logicBoard[ii],    // its logical volume				  
				       boardName[ii],       // its name
				       logicRandB[ii],      // its mother  volume
				       false,           // no boolean operations
				       ii);              // copy number 
  }


  //------------------------------ 
  // protection caps
  //------------------------------


  solidOutCap = new G4Box("OutCap",0.5*capOutSize,0.5*capOutSize,0.5*capOutHight);
  solidInCap = new G4Box("InCap",0.5*capInSize,0.5*capInSize,0.5*capInHight);

  
  //  G4RotationMatrix *cap_rot = new G4RotationMatrix(G4ThreeVector(0., 0., 0.0),angle*CLHEP::deg);
  G4RotationMatrix cap_rot(G4ThreeVector(0., 0., 0.0),angle*CLHEP::deg);
  G4ThreeVector translation(0, 0, - 0.5*(capOutHight - capInHight));

  G4Transform3D transform(cap_rot,translation);
  cap = new G4SubtractionSolid("cap",solidOutCap,solidInCap,transform);

  G4String CAPName[16];
  CAPName[0] = "CAP-0";
  CAPName[1] = "CAP-1";
  CAPName[2] = "CAP-2";
  CAPName[3] = "CAP-3";
  CAPName[4] = "CAP-4";
  CAPName[5] = "CAP-5";
  CAPName[6] = "CAP-6";
  CAPName[7] = "CAP-7";
  CAPName[8] = "CAP-8";
  CAPName[9] = "CAP-9";
  CAPName[10] = "CAP-10";
  CAPName[11] = "CAP-11";
  CAPName[12] = "CAP-12";
  CAPName[13] = "CAP-13";
  CAPName[14] = "CAP-14";
  CAPName[15] = "CAP-15";

  G4ThreeVector posT3[16];


  for(size_t ii=0; ii<16; ii++) posT3[ii] = G4ThreeVector(0.0*CLHEP::mm, 0.0*CLHEP::mm, tiltedTelescopePosition + ChamberSpacing*G4double(ii+1)+(0.5*capOutHight));


  for(size_t ii=0; ii<8; ii++){

    logicCap[ii] = new G4LogicalVolume(cap , Air, "Cap",0,0,0);
    physiCap[ii] = new G4PVPlacement(mirror_rot,              // rotated
				     posT3[ii],           // at (x,y,z)
				     logicCap[ii],   // its logical volume				  
				     CAPName[ii],    // its name
				     logicWorld,         // its mother  volume
				     false,              // no boolean operations
				     ii);                // copy number 
  }

  for(size_t ii=8; ii<16; ii++){

    logicCap[ii] = new G4LogicalVolume(cap , Air, "Cap",0,0,0);
    physiCap[ii] = new G4PVPlacement(0,              // no rotation
				     posT3[ii],           // at (x,y,z)
				     logicCap[ii],   // its logical volume				  
				     CAPName[ii],    // its name
				     logicWorld,         // its mother  volume
				     false,              // no boolean operations
				     ii);                // copy number 
  }




  //------------------------------ 
  // Wall
  //------------------------------
  solidOutWall = new G4Box("outwall",0.5*WallOutX,0.5*WallOutY,0.5*WallOutZ);
  solidInWall = new G4Box("inwall",0.5*WallInX,0.5*WallInY,0.5*WallInZ);


  solidWall = new G4SubtractionSolid("wall",solidOutWall,solidInWall);
  logicWall = new G4LogicalVolume(solidWall,PVC,"Wall",0,0,0);
  physiWall = new G4PVPlacement(0,
				G4ThreeVector(BoxCenterXPosition,BoxCenterYPosition,BoxCenterZPosition),
				logicWall,
				"WALL",
				logicWorld,
				false,0);


  //------------------------------ 
  // Insulation
  //------------------------------
  solidOutInsulation = new G4Box("outinsulation",0.5*InsulationOutX,0.5*InsulationOutY,0.5*InsulationOutZ);
  solidInInsulation = new G4Box("ininsulation",0.5*InsulationInX,0.5*InsulationInY,0.5*InsulationInZ);
  
  solidInsulation = new G4SubtractionSolid("insulation",solidOutInsulation,solidInInsulation);
  logicInsulation = new G4LogicalVolume(solidInsulation,insulation,"Insulation",0,0,0);
  physiInsulation = new G4PVPlacement(0,
				      G4ThreeVector(0,0,0),
				      logicInsulation,
				      "INSULATION",
				      logicWall,
				      false,0);


  //------------------------------ 
  // Peltiers
  //------------------------------


  solidPeltier = new G4Box("peltier",0.5*PeltierX,0.5*PeltierY,0.5*PeltierZ);

  G4String PELTIERName[2];
  PELTIERName[0] = "PELTIER-Up";
  PELTIERName[1] = "PELTIER-Down";
 
  G4ThreeVector posT4[2];

  posT4[0] = G4ThreeVector(0,0.5*(InsulationInY - PeltierY) ,0);
  posT4[1] = G4ThreeVector(0,-0.5*(InsulationInY - PeltierY) ,0);



  for(size_t ii=0; ii<2; ii++){

    logicPeltier[ii] = new G4LogicalVolume(solidPeltier , Alu, "Peltier");
    physiPeltier[ii] = new G4PVPlacement(0,              // no rotation
				     posT4[ii],           // at (x,y,z)
				     logicPeltier[ii],   // its logical volume	    		  
				     PELTIERName[ii],    // its name
				     logicWall,         // its mother  volume
				     false,              // no boolean operations
				     ii);                // copy number 
  }




  //------------------------------
  // floor
  //------------------------------

  solidFloor = new G4Box("floor",0.5*FloorX,0.5*FloorY,0.5*FloorZ);
  logicFloor = new G4LogicalVolume(solidFloor,Alu,"Floor",0,0,0);
  physiFloor = new G4PVPlacement(0,
				 G4ThreeVector(0, FloorPositionY,FloorPositionZ),
				 logicFloor,
				 "FLOOR",
				 logicWorld,
				 false,0);


  //------------------------------
  // scintillator
  //------------------------------
  G4double rotangle = 90.;
  G4RotationMatrix *ScintRot = new G4RotationMatrix(G4ThreeVector(0, 1.0, 0.0),rotangle*CLHEP::deg);


  solidScint = new G4Tubs("scintilator",ScintInnerRad, ScintOuterRad, ScintLength, ScintStartAngle, ScintSegmentAngle);
  logicScint = new G4LogicalVolume(solidScint,mylar,"Scintilator");
  physiScint = new G4PVPlacement(ScintRot,
				 G4ThreeVector(0,0,ScintPosition),
				 logicScint,
				 "SCINTILATOR",
				 logicWorld,
				 false,0);

 
	 
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "pixelTB/TrackerChamberSD";
  pixelTBTrackerSD* aTrackerSD = new pixelTBTrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  for(size_t ii=0; ii<16; ii++) logicChamber[ii]->SetSensitiveDetector( aTrackerSD );




//--------- Visualization attributes -------------------------------


  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour::White());
  BoxVisAtt->SetForceWireframe(true);
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
  //for(size_t ii=0; ii<8; ii++)  logicTracker[ii]->SetVisAttributes(BoxVisAtt);
  for(size_t ii=0; ii<8; ii++)  logicRandB[ii]->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* WallVisAtt = new G4VisAttributes(G4Colour::Magenta());
  WallVisAtt->SetForceWireframe(true);
  logicWall->SetVisAttributes(WallVisAtt);

  G4VisAttributes* InsulationVisAtt = new G4VisAttributes(G4Colour::Yellow());
  InsulationVisAtt->SetForceWireframe(true);
  logicInsulation->SetVisAttributes(InsulationVisAtt);

  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour::Yellow());
  ChamberVisAtt->SetForceWireframe(true);
  for(size_t ii=0; ii<16; ii++)    logicChamber[ii]->SetVisAttributes(ChamberVisAtt);

  G4VisAttributes* ROCVisAtt = new G4VisAttributes(G4Colour::Yellow());
  for(size_t ii=0; ii<16; ii++)    logicROC[ii]->SetVisAttributes(ROCVisAtt);

  G4VisAttributes* BoardVisAtt = new G4VisAttributes(G4Colour::Green());
  for(size_t ii=0; ii<16; ii++)    logicBoard[ii]->SetVisAttributes(BoardVisAtt);

  G4VisAttributes* CapVisAtt = new G4VisAttributes(G4Colour::Magenta());
  CapVisAtt->SetForceWireframe(true);
  for(size_t ii=0; ii<16; ii++)    logicCap[ii]->SetVisAttributes(CapVisAtt);
  
  G4VisAttributes* AluVisAtt = new G4VisAttributes(G4Colour::Cyan());
  for(size_t ii=0; ii<2; ii++)    logicPeltier[ii]->SetVisAttributes(AluVisAtt);
  
  logicFloor->SetVisAttributes(AluVisAtt);

//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //

  G4double maxStep = 0.5*ChamberWidth;
  stepLimit = new G4UserLimits(maxStep);
  for(size_t ii=0; ii<16; ii++)    logicTracker[ii]->SetUserLimits(stepLimit);
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void pixelTBDetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target is " << fTargetLength/CLHEP::cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBDetectorConstruction::setChamberMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {ChamberMater = pttoMaterial;
       for(size_t ii=0; ii<16; ii++)  logicChamber[ii]->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The chambers are " << ChamberWidth/CLHEP::cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void pixelTBDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetMagFieldValue(fieldValue);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
