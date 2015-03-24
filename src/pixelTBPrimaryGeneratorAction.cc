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
// $Id: pixelTBPrimaryGeneratorAction.cc,v 1.7 2006-06-29 17:48:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "pixelTBPrimaryGeneratorAction.hh"
#include "pixelTBDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "HistoManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pixelTBPrimaryGeneratorAction::pixelTBPrimaryGeneratorAction(pixelTBDetectorConstruction* myDC, HistoManager * histo )
  :myDetector(myDC),n_particle(0),fHistoManager(histo)
{

  G4int n_particle_per_shot = 1;
  
  //G4cout << "n_particle is " << n_particle << G4endl;

  particleGun = new G4ParticleGun(n_particle_per_shot);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  
  particleGun->SetParticleDefinition(particle);
//   particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(120.0*CLHEP::GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

pixelTBPrimaryGeneratorAction::~pixelTBPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void pixelTBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //  G4int n_particle;
  //  n_particle = CLHEP::RandLandau::shoot();
  //n_particle = 2;
  n_particle = 1;//fHistoManager->NumberOfP;
  //G4cout << "xxn_particle is " << n_particle << G4endl;

  for(int i =0; i<n_particle; i++){
  G4double position = -0.5*(myDetector->GetWorldFullLength());
  particleGun->SetParticlePosition(G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,position));
  

  G4double x0 = -1.*CLHEP::mm; //-0.5*(Detector->GetWorldSizeX());
  G4double y0 = -1.*CLHEP::mm;
  G4double z0 = 0.*CLHEP::mm;

  x0 = CLHEP::cm*(G4UniformRand()-0.5);
  y0 = CLHEP::cm*(G4UniformRand()-0.5);
  while(x0 > 5. || x0 < -5.){
    x0 = 2*CLHEP::mm*(G4RandGauss::shoot(-2.21,6.713));
  }
  while(y0 > 5. || y0 < -5.){
    y0 = 2*CLHEP::mm*(G4RandGauss::shoot(0.19,5.851));
  }
    particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  //  particleGun->SetNumberOfParticles(4);

  G4double xdir = 0.;
  G4double ydir = 0.;
  G4double telescopelength = 7 * 16*CLHEP::mm;
  xdir = 1*CLHEP::mm * G4RandGauss::shoot(0,0.1);
  ydir = 1*CLHEP::mm * G4RandGauss::shoot(0,0.055);
  //  G4cout << " particle dir " << xdir << " / " << ydir << G4endl;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(xdir, ydir, telescopelength));
 
  particleGun->GeneratePrimaryVertex(anEvent);
  }
  //G4cout << "Done" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

