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
// $Id: pixelTBChamberParameterisation.cc,v 1.9 2006-06-29 17:47:58 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "pixelTBChamberParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pixelTBChamberParameterisation::pixelTBChamberParameterisation(  
        G4int    NoChambers, 
        G4double startZ,          //  Z of center of first 
        G4double spacingZ,        //  Z spacing of cesnters
        G4double widthChamber, 
        G4double lengthInitial, 
        G4double lengthFinal )
{
   fNoChambers =  NoChambers; 
   fStartZ     =  startZ; 
   fHalfWidth  =  widthChamber*0.5;
   fSpacing    =  spacingZ;
//    fHalfLengthFirst = 0.5 * lengthInitial; 
   fHalfLengthFirst = lengthInitial; 
//    fHalfLengthLast = lengthFinal;
   fHalfLengthIncr = 0.;
//    if( NoChambers > 0 ){
//       fHalfLengthIncr =  0.9 * (lengthFinal-lengthInitial)/NoChambers;
//       if (spacingZ < widthChamber) {
//        G4Exception("pixelTBChamberParameterisation construction: Width>Spacing");
//       }
//    }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pixelTBChamberParameterisation::~pixelTBChamberParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBChamberParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  const G4int   nCol = 52;
  const G4int   nRow = 80;
  const G4int   npixel = nCol*nRow;
  G4double      Zposition = 0;
  G4double      Xposition = -0.15*CLHEP::mm*nCol/2 -0.075*CLHEP::mm + ((copyNo%npixel)%nCol+1)*0.15*CLHEP::mm; //offset 5 pixel of 150um
  G4double      Yposition = -0.1*CLHEP::mm*nRow/2 -0.05*CLHEP::mm + ((copyNo%npixel)/nCol+1)*0.1*CLHEP::mm;    //offset 5 pixel of 100um

  G4ThreeVector origin(Xposition, Yposition, Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBChamberParameterisation::ComputeDimensions
(G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
{
  G4double  halfLength= fHalfLengthFirst;// + copyNo * fHalfLengthIncr;
  trackerChamber.SetXHalfLength(halfLength*1.5);
  trackerChamber.SetYHalfLength(halfLength);
  trackerChamber.SetZHalfLength(fHalfWidth);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
