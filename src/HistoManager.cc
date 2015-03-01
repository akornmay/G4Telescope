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
/// \file analysis/AnaEx02/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 74272 2013-10-02 14:48:50Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "pixelTBTrackerHit.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fChain(0),fTurn(0),fBucket(0)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
  G4String fileName = "simdataTree_PixelTestBoard1_RUN";
  fileName.append(std::to_string(fRunNumber));
  fileName.append("_TURNSTART");
  fileName.append(std::to_string(firstTurn));
  fileName.append(".root");

  fOutFileStraight = new TFile(fileName,"RECREATE");
  if(!fOutFileStraight) {
  G4cout << " HistoManager::book :" 
	 << " problem creating the ROOT TFile "
	 << G4endl;
  return;
  }

  G4String title = "Telescope";
  fTreeStraight = new TTree(title,title);

  fTreeStraight->Branch("Event", &eventNr, "Event/i");
  fTreeStraight->Branch("roc",&roc,"roc/I");
  fTreeStraight->Branch("pulseHeight", &pulseHeight, "pulseHeight/D");
  fTreeStraight->Branch("energy", &energy,"energy/D");
  fTreeStraight->Branch("vcal", &vcal, "vcal/I");
  fTreeStraight->Branch("col", &col, "col/I");
  fTreeStraight->Branch("row", &row, "row/I");
  fTreeStraight->Branch("flux", &flux, "flux/F");
  fTreeStraight->Branch("x", &x, "x/F");
  fTreeStraight->Branch("y", &y, "y/F");

  G4cout << "Booked output file for straight telescope" << G4endl;

  fileName = "simdataTree_PixelTestBoard2_RUN";
  fileName.append(std::to_string(fRunNumber));
  fileName.append("_TURNSTART");
  fileName.append(std::to_string(firstTurn));
  fileName.append(".root");



  fOutFileTilted = new TFile(fileName,"RECREATE");
  if(!fOutFileTilted) {
  G4cout << " HistoManager::book :" 
	 << " problem creating the ROOT TFile "
	 << G4endl;
  return;
  }

  
  fTreeTilted = new TTree(title,title);

  fTreeTilted->Branch("Event", &eventNr, "Event/i");
  fTreeTilted->Branch("roc",&roc,"roc/I");
  fTreeTilted->Branch("pulseHeight", &pulseHeight, "pulseHeight/D");
  fTreeTilted->Branch("energy", &energy,"energy/D");
  fTreeTilted->Branch("vcal", &vcal, "vcal/I");
  fTreeTilted->Branch("col", &col, "col/I");
  fTreeTilted->Branch("row", &row, "row/I");
  fTreeTilted->Branch("flux", &flux, "flux/F");
  fTreeTilted->Branch("x", &x, "x/F");
  fTreeTilted->Branch("y", &y, "y/F");

  G4cout << "Booked output file for tilted telescope" << G4endl;
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::open()
{
  G4String QIEFileName = fQieDir + "RawData_spill";
  QIEFileName.append(std::to_string(fRunNumber));
  QIEFileName.append(".bin.root");
  //check if file exists
  //std::ifstream is(QIEFileName,std::ios::in);
  std::istringstream is(QIEFileName.c_str());
  if(!is)
    {
      G4cout << "Error: File \"" << QIEFileName << "\" doesn't exist."  << G4endl;
    }


  G4cout << "Reading from ROOT tree file with the name: " << QIEFileName << G4endl; 
  //open the right root tree file
  G4cout << "Creating new TChain " << G4endl; 
  fChain = new TChain("tree_QIE");
  G4cout << "Adding file to chain" << G4endl; 
  fChain->Add(QIEFileName.c_str());
  //disable all branches
  G4cout << "Doing brach business " << G4endl; 
  fChain->SetBranchStatus("*",false);
  //only load the branch we really need
  fChain->SetBranchStatus("BeamIntensity",true);
  fChain->SetBranchAddress("BeamIntensity",BeamIntensity);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::GetNumbersOfParticlesInEvent(){

  fChain->GetEntry(fTurn);
  G4double temp;
  if(fTransparentMode)
    {
      // when simulating the transparent mode we only need to simulate one bucket of every turn
      
      temp = BeamIntensity[fTriggerBucket];
      ++fTurn;
    }
  else
    {
      //when running in full simulation mode we need to simulate every bucket of a turn
  
      temp = BeamIntensity[fBucket];
      //G4cout << "Evaluating Turn: " << fTurn << " Bucket: " << fBucket << " BeamIntensity: " << temp <<  G4endl;
            ++fBucket;
      if(fBucket >588)
	{
	  fBucket = 0; ++fTurn;
	}
    }
      // NumberOfP = (int)((temp - BeamIntensity[32] +260.)/320.);
  NumberOfP = (int)((temp - BeamIntensity[32] +260.)/320. + 0.5) * 10;

  TotalNumberOfParticles += NumberOfP;
  
  // G4cout << "number of particles " <<  NumberOfP << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
 if (fOutFileTilted) {
    fOutFileTilted->Write();       // Writing the histograms to the file
    fOutFileTilted->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Tilted tree saved \n" << G4endl;
  }

 if (fOutFileStraight) {
    fOutFileStraight->Write();       // Writing the histograms to the file
    fOutFileStraight->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Straight tree saved \n" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddHit(pixelTBTrackerHit* Hit,G4int EventNumber, G4int &returnROC)
{
  //energy
  energy = Hit->GetEdep();

  //  G4double xray_offset = -600. ;
  //  G4double xray_slope = 50. ;
  //  G4double ionisation_energy = 3.62; //[eV/e-]
  vcal = (int)(energy / 3.62 * 1000.0 - (-600.0)) / 50.0;
  
  //  double adc = (p3 + p2 * TMath::TanH(p0*vcal - p1));
  //p0=+3.122759e-03, p1=+1.085111e+00, p2=+1.036756e+02, p3=+1.520615e+02);
  pulseHeight = (+1.520615e+02) + (+1.036756e+02) * TMath::TanH((+3.122759e-03)*vcal - (+1.085111e+00)) ;

  //  position
  G4ThreeVector pos = Hit->GetPos();
  x = pos.getX(); //get's returned in mm
  y = pos.getY(); //get's returned in mm
  z = pos.getZ(); //get's returned in mm
  //
  //  col = 52 / 2 + x / 0.15;
  //  row = 80 / 2 + y / 0.1 ;
  col = Hit->GetChamberNb()%52;
  row = Hit->GetChamberNb()/52;
  returnROC = (int) (z/16.0 + 0.5) - 11;
  eventNr = EventNumber;

  if(returnROC < 8)
    {
      roc = returnROC;
      //      G4cout << "Filling tilted" << G4endl;
      fTreeTilted->Fill();
    }
  else
    {
      roc = returnROC%8;  
      //      G4cout << "Filling straight" << G4endl;
      fTreeStraight->Fill();
    }

}


void HistoManager::AddEmptyEvent(G4int EventNumber, G4int PixelTestBoard)
{

 //energy
  pulseHeight = energy = -1;
  vcal = -1;
  //  position
  x = y = z = -1;
  //
  col = -1;
  row = -1;
  roc = -1;

  //flux
  flux = -1;

  eventNr = EventNumber;

 if(PixelTestBoard == 2)
    {
      G4cout << "Filling tilted" << G4endl;
      fTreeTilted->Fill();
    }
 if(PixelTestBoard == 1)
    {
      G4cout << "Filling straight" << G4endl;
      fTreeStraight->Fill();
    }

}


G4int HistoManager::GetFirstTurn()
{
  return firstTurn;
}

void HistoManager::SetFirstTurn(G4int first)
{

  fTurn = firstTurn = first;
  eventNr = firstTurn * 588;

}

G4int HistoManager::GetRunNumber()
{
  return fRunNumber;
}

void HistoManager::SetRunNumber(G4int runNumber)
{
  fRunNumber = runNumber;
}

void HistoManager::SetQieDir(G4String qieDir)
{
  fQieDir = qieDir;
}

void HistoManager::SetTriggerBucket(G4int TRIGGER_BUCKET)
{
  fTriggerBucket = TRIGGER_BUCKET;

  if(fTransparentMode &&  (fTriggerBucket < 0 || fTriggerBucket > 588))
    {
      G4cout << "The configuration points to a non-existend bucket number" << G4endl;
      G4cout << "Ending simulation here " << G4endl;
      exit(0);
    }
}

G4int HistoManager::GetTriggerBucket()
{
  return fTriggerBucket;
}
