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
// $Id: pixelTBDetectorConstruction.hh,v 1.10 2008-09-22 16:41:20 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef pixelTBDetectorConstruction_h
#define pixelTBDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "pixelTBMagneticField.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class pixelTBDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class pixelTBDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     pixelTBDetectorConstruction();
    ~pixelTBDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const 
     //G4VPhysicalVolume* GetTracker() {return physiTracker;};
     G4double GetTrackerFullLength() {return fTrackerLength;};
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     
     void setTargetMaterial (G4String);
     void setChamberMaterial(G4String);
     void SetMagField(G4double);
     void SetMaxStep (G4double);     
  
     G4double fCSZonePercentageX;
     G4double fCSZonePercentageY;

     G4ThreeVector origROC[16];      //this will store all the origins of the ROCs
     G4ThreeVector upleftROC[16];      //this will store all the origins of the ROCs
     G4ThreeVector lowrightROC[16];      //this will store all the origins of the ROCs
     G4ThreeVector GetNormalVector(G4int ROC);


  private:

     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
     
     G4Box*             solidTarget;   // pointer to the solid Target
     G4LogicalVolume*   logicTarget;   // pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
               
     G4Box*             solidTracker;     // pointer to the solid Tracker
     G4LogicalVolume*   logicTracker[16];  // pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker[16];  // pointer to the physical Tracker
     
     G4Box*             solidChamber;     // pointer to the solid Chamber
     G4LogicalVolume*   logicChamber[16];  // pointer to the logical Chamber
     G4VPhysicalVolume* physiChamber[16];  // pointer to the physical Chamber

     G4Box*             solidRandB;     // pointer to the solid ROC and Board volume
     G4LogicalVolume*   logicRandB[16];  // pointer to the logical ROC and Board volume
     G4VPhysicalVolume* physiRandB[16];  // pointer to the physical ROC and Board volume

     G4Box* solidOutCap;
     G4Box* solidInCap;
     G4SubtractionSolid* cap;
     G4LogicalVolume* logicCap[16];
     G4VPhysicalVolume* physiCap[16];

     G4Box*             solidROC;     // pointer to the solid ROC
     G4LogicalVolume*   logicROC[16];  // pointer to the logical ROC
     G4VPhysicalVolume* physiROC[16];  // pointer to the physical ROC

     G4Box*             solidBoard;     // pointer to the solid Board
     G4LogicalVolume*   logicBoard[16];  // pointer to the logical Board
     G4VPhysicalVolume* physiBoard[16];  // pointer to the physical Board

     G4Box* solidOutWall;
     G4Box* solidInWall;
     G4SubtractionSolid* solidWall;    // pointer to the solid Wall
     G4LogicalVolume*   logicWall;     // pointer to the logical Wall
     G4VPhysicalVolume* physiWall;     // pointer to the physical Wall

     G4Box* solidOutInsulation;
     G4Box* solidInInsulation;
     G4SubtractionSolid* solidInsulation;    // pointer to the solid Wall
     G4LogicalVolume*   logicInsulation;     // pointer to the logical Insulation
     G4VPhysicalVolume* physiInsulation;     // pointer to the physical Insulation

     G4Box* solidPeltier;
     G4LogicalVolume* logicPeltier[2];
     G4VPhysicalVolume* physiPeltier[2];

     G4Box*             solidFloor;     // pointer to the solid Floor
     G4LogicalVolume*   logicFloor;     // pointer to the logical Floor
     G4VPhysicalVolume* physiFloor;     // pointer to the physical Floor

     G4Tubs*            solidScint;     // pointer to the solid Scintilator
     G4LogicalVolume*   logicScint;     // pointer to the logical Scintilator
     G4VPhysicalVolume* physiScint;     // pointer to the physical Scintilator


     G4Material*         TargetMater;  // pointer to the target  material
     G4Material*         ChamberMater; // pointer to the chamber material
     G4Material*         BoardMater;   // pointer to the board material

     G4VPVParameterisation* chamberParam; // pointer to chamber parameterisation
     G4UserLimits* stepLimit;             // pointer to user step limits

     pixelTBMagneticField* fpMagField;   // pointer to the magnetic field 
     
     pixelTBDetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume
     G4double fTargetLength;           // Full length of Target
     G4double fTrackerLength;          // Full length of Tracker
     G4int    NbOfChambers;            // Nb of chambers in the tracker region
     G4double ChamberWidth;            // width of the chambers
     G4double ChamberSize;            // size of the chambers
     G4double ROCSize;                // size of the ROC
     G4double ROCWidth;                // width of the ROC
     G4double BoardWidth;                // width of the boards
     G4double ChamberSpacing;	       // distance between chambers
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
