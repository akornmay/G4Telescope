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
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pixelTBChamberParameterisation::pixelTBChamberParameterisation(  
        G4int    NoChambers, 
        G4double startZ,          //  Z of center of first 
        G4double spacingZ,        //  Z spacing of cesnters
        G4double widthChamber, 
        G4double lengthInitial, 
        G4double lengthFinal,
	G4double CSZoneX,
	G4double CSZoneY )
{
   fNoChambers =  NoChambers; 
   fStartZ     =  startZ; 
   fHalfWidth  =  widthChamber*0.5;
   fSpacing    =  spacingZ;
//    fHalfLengthFirst = 0.5 * lengthInitial; 
   fHalfLengthFirst = lengthInitial; 
//    fHalfLengthLast = lengthFinal;
   fHalfLengthIncr = 0.;
   fCSZonePercentageX = CSZoneX;
   fCSZonePercentageY = CSZoneY;
   G4cout << sizeof(myPosVectors) << G4endl;
   for(G4int ii = 0; ii < (4160*9) ; ++ii)
     {
       G4cout << ii << G4endl;
       myPosVectors[ii] = NULL;
       myChamberSizes[ii].set = false;

     }




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

  if (myPosVectors[copyNo]) {
    physVol->SetTranslation(*(myPosVectors[copyNo]));
    physVol->SetRotation(0);
  } else {

  //old
  //const G4int   nCol = 52;
  //const G4int   nRow = 80;
  const G4int   nCol = 52 * 3;
  const G4int   nRow = 80 * 3;
  const G4int   npixel = nCol*nRow;
  const G4double      Zposition = 0;

  // Assuming copyNo < npixel
  //  G4int& pixelIndex = copyNo; // otherwise pixelIndex = copyNo%npixel
  //  G4int rowIndex = copyNo/nCol;
  //  G4int colIndex = copyNo%nCol;
  div_t rowCol = div(copyNo, nCol);

  G4int rowIndex = rowCol.quot; // only if it will be optimized away
  G4int colIndex = rowCol.rem;  // only if it will be optimized away
 
  G4double Xposition = -0.15*CLHEP::mm*(nCol/3)/2 -0.075*CLHEP::mm + ((colIndex)/3+1)*0.15*CLHEP::mm;
      
  G4double Yposition = -0.10*CLHEP::mm*(nRow/3)/2 -0.050*CLHEP::mm + ((rowIndex)/3+1)*0.1*CLHEP::mm;

  // old version
  /*
  G4double      Xposition = -0.15*CLHEP::mm*nCol/2 -0.075*CLHEP::mm + (pixelIndex%nCol+1)*0.15*CLHEP::mm; //offset 5 pixel of 150um
  G4double      Yposition = -0.1*CLHEP::mm*nRow/2 -0.05*CLHEP::mm + (pixelIndex/nCol+1)*0.1*CLHEP::mm;    //offset 5 pixel of 100um
  */

  if(colIndex%3 == 0) //check the column
    {
      // position in X stays the same 
    }
  else if(colIndex%3 == 1)
    {
      Xposition = Xposition + 0.5*(1-fCSZonePercentageX)*0.15*CLHEP::mm;     
    }
  else if(colIndex%3 == 2)
    {
      Xposition = Xposition + (1-fCSZonePercentageX)*0.15*CLHEP::mm;;
    }

  if(rowIndex%3 == 0)  //check the row
    {
      // position in Y also stays the same
    }
  else if(rowIndex%3 == 1)
    {
      Yposition = Yposition + 0.5*(1-fCSZonePercentageY)*0.1*CLHEP::mm;
    }
  else if(rowIndex%3 == 2)
    {
      Yposition = Yposition + (1-fCSZonePercentageY)*0.1*CLHEP::mm; 
    }


  myPosVectors[copyNo] = new G4ThreeVector(Xposition, Yposition, Zposition);
  physVol->SetTranslation(*(myPosVectors[copyNo]));
  physVol->SetRotation(0);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBChamberParameterisation::ComputeDimensions
(G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
{

  if (myChamberSizes[copyNo].set) {
    trackerChamber.SetXHalfLength(myChamberSizes[copyNo].halfX);
    trackerChamber.SetYHalfLength(myChamberSizes[copyNo].halfY);
    trackerChamber.SetZHalfLength(myChamberSizes[copyNo].halfZ);
  } else {

    G4double  halfLengthY= fHalfLengthFirst;
    G4double  halfLengthX= fHalfLengthFirst * 1.5;

    const G4int   nCol = 52 * 3;
    const G4int   nRow = 80 * 3;
    
    myChamberSizes[copyNo].halfZ =  fHalfWidth;

    
    if((copyNo%nCol)%3 == 0) //check the column
      {
	myChamberSizes[copyNo].halfX = halfLengthX*fCSZonePercentageX;
      }
    else if((copyNo%nCol)%3 == 1)
      {
	myChamberSizes[copyNo].halfX = halfLengthX*(1-2*fCSZonePercentageX);
      }
    else if((copyNo%nCol)%3 == 2)
      {
	myChamberSizes[copyNo].halfX = halfLengthX*fCSZonePercentageX;
      }
    
    if((copyNo/nCol)%3 == 0)  //check the row
      {
	myChamberSizes[copyNo].halfY = halfLengthY*fCSZonePercentageY;
      }
    else if((copyNo/nCol)%3 == 1)
      {
	myChamberSizes[copyNo].halfY = halfLengthY*(1-2*fCSZonePercentageY);
      }
    else if((copyNo/nCol)%3 == 2)
      {
	myChamberSizes[copyNo].halfY = halfLengthY*fCSZonePercentageY;
      }

    trackerChamber.SetXHalfLength(myChamberSizes[copyNo].halfX);
    trackerChamber.SetYHalfLength(myChamberSizes[copyNo].halfY);
    trackerChamber.SetZHalfLength(myChamberSizes[copyNo].halfZ);
    myChamberSizes[copyNo].set = true;
      

 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*

 -----------------------------------
|      |                     |      |
|  02  |                     |      |
|      |                     |      |
 -----------------------------------
|      |                     |      |
|      |                     |      |
|      |                     |      |
|      |                     |      |
|  01  |                     |      |
|      |                     |      |
|      |                     |      |
|      |                     |      |
 -----------------------------------
|      |                     |      |
|  00  |         10          |  20  |
|      |                     |      |
 -----------------------------------
 

 */
