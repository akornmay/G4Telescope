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
/// \file analysis/AnaEx02/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 74272 2013-10-02 14:48:50Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <TH2D.h>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TFile;
class TTree;
class TH1D;
class TChain;
class pixelTBTrackerHit;

const G4int MaxHisto = 5;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
public:
  
  HistoManager();
  ~HistoManager();
  
  void book();
  void open();
  void save();

  void AddHit(pixelTBTrackerHit* Hit, G4int EventNumber, G4int &returnROC);
  void AddEmptyEvent(G4int EventNumber, G4int PixelTestBoard);

  void CollectHits(pixelTBTrackerHit* Hit, G4int EventNumber, std::map<int,pixelTBTrackerHit>& mapofHits);
  void AddHits(G4double pixelArray[16][52][80], G4int EventNumber);


  void GetNumbersOfParticlesInEvent();
 
  G4int NumberOfP;
  G4long TotalNumberOfParticles; 

  void SetFirstTurn(G4int first);
  G4int GetFirstTurn();

  void SetRunNumber(G4int runNumber);
  G4int GetRunNumber();

  void SetQieDir(G4String qieDir);

 void SetROCVectors(G4ThreeVector orROC[16], G4ThreeVector upleftROC[16], G4ThreeVector lowrightROC[16]);
  G4double CheckHitDistance(pixelTBTrackerHit* Hit);
  //std::pair<G4double,G4double> CheckPixel(pixelTBTrackerHit* Hit);
  void CheckPixel(pixelTBTrackerHit* Hit);
  void getHitmap(pixelTBTrackerHit* Hit);  

  G4ThreeVector projectOnSensorPlane(pixelTBTrackerHit* Hit);

  double BEAMINTENSITY_PARAM1;		   
  double BEAMINTENSITY_PARAM2;		   
  double BEAMINTENSITY_SCALING;		   

  void SetTransparentMode(G4bool transparentMode){fTransparentMode = transparentMode;};
  void SetTriggerBucket(G4int TRIGGER_BUCKET);
  G4int GetTriggerBucket();

  void AddHits3by3(pixelTBTrackerHit* Hit, G4int EventNumber, G4int &returnROC);

private:

  TFile* fOutFileTilted;
  TFile* fOutFileStraight;

  TTree* fTreeTilted;
  TTree* fTreeStraight;
  
  G4int roc, col, row, vcal;// adc;
  G4double pulseHeight, energy;
  G4float flux;
  G4int eventNr;
  G4double x, y, z;

  TChain* fChain;                      //this is the root file we want to read
  G4float BeamIntensity[588];         //this is where the beam intensity of 
  G4int fTurn;
  G4int fBucket;

  G4int firstTurn;
  G4int fRunNumber;
  G4String fQieDir;

  G4bool fTransparentMode;
  G4int fTriggerBucket;

  G4ThreeVector origROC[16];
  G4ThreeVector UnormROC[16]; //normed vector pointing from row 0 to row 80
  G4ThreeVector VnormROC[16]; //normed vector pointing from col 0 to col 52
  G4ThreeVector NnormROC[16]; //normed vector pointing out of the plane of the roc

  TH1D* eta_straight;
  TH1D* eta_tilted;

  TH2D* reducedHitmap;



};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

