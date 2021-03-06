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
  fTreeStraight->Branch("flux", &flux, "flux/D");
  fTreeStraight->Branch("x", &x, "x/D");
  fTreeStraight->Branch("y", &y, "y/D");

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
  fTreeTilted->Branch("flux", &flux, "flux/D");
  fTreeTilted->Branch("x", &x, "x/D");
  fTreeTilted->Branch("y", &y, "y/D");

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
      //  G4cout << "in Turn " << fTurn << " BeamIntensity[" << fTriggerBucket << "]=" << temp << G4endl;
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
  NumberOfP = (int)(((temp - BeamIntensity[32] + BEAMINTENSITY_PARAM1)/BEAMINTENSITY_PARAM2) * BEAMINTENSITY_SCALING);

  TotalNumberOfParticles += NumberOfP;
  
  //  G4cout << "number of particles " <<  NumberOfP << G4endl;
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
  energy = Hit->GetEdep()/CLHEP::keV;

  //  G4double xray_offset = -600. ;
  //  G4double xray_slope = 50. ;
  //  G4double ionisation_energy = 3.62; //[eV/e-]
  vcal = (int)(energy / 3.62 * 1000.0 - (-600.0)) / 50.0;
  
  //  double adc = (p3 + p2 * TMath::TanH(p0*vcal - p1));
  //p0=+3.122759e-03, p1=+1.085111e+00, p2=+1.036756e+02, p3=+1.520615e+02);
  pulseHeight = (+1.520615e+02) + (+1.036756e+02) * TMath::TanH((+3.122759e-03)*vcal - (+1.085111e+00)) ;

  //  position
  G4ThreeVector pos = Hit->GetPos();
  z = pos.getZ(); //get's returned in mm
  returnROC = (int) (z/16.0 + 0.5) - 11;
  
  

  G4ThreeVector onPlaneHit = projectOnSensorPlane(Hit);
  x = onPlaneHit.getX(); //get's returned in mm
  y = onPlaneHit.getY(); //get's returned in mm
  z = onPlaneHit.getZ(); //get's returned in mm

  eventNr = EventNumber;

  col = Hit->GetChamberNb()%(52*3)/3;
  row = Hit->GetChamberNb()/(52*3)/3;


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

void HistoManager::AddHits3by3(pixelTBTrackerHit* Hit, G4int EventNumber, G4int &returnROC)
{
  //first we check where the hit was
  G4ThreeVector pos = Hit->GetPos();
  z = pos.getZ(); //get's returned in mm
  returnROC = (int) (z/16.0 + 0.5) - 11;

  G4ThreeVector onPlaneHit = projectOnSensorPlane(Hit);


  x = onPlaneHit.getX(); //get's returned in mm
  y = onPlaneHit.getY(); //get's returned in mm
  z = onPlaneHit.getZ(); //get's returned in mm

  //now let's get the pixel number






}


void HistoManager::AddHits(G4double pixelArray[16][52][80], G4int EventNumber)
{


  G4int minRoc = 16;
  G4int maxRoc = -1;
  //loop over all pixels
  for(roc = 0; roc < 16; ++roc)
    {

      for(col = 0; col < 52; ++col)
	{

	  for(row = 0; row < 80; ++row)
	    {

	      if(pixelArray[roc][col][row] != 0)
		{
		  energy = pixelArray[roc][col][row];
		  vcal = (int)(energy / 3.62 * 1000.0 - (-600.0)) / 50.0;
		  pulseHeight = (+1.520615e+02) + (+1.036756e+02) * TMath::TanH((+3.122759e-03)*vcal - (+1.085111e+00)) ;
		  x=0;
		  y=0;
		  eventNr = EventNumber;

		  if(roc < 8)
		    {
		      fTreeTilted->Fill();
		      if(roc < minRoc) minRoc = roc;
		      if(roc > maxRoc) maxRoc = roc;
		    }
		  else
		    {
		      fTreeStraight->Fill();
		      if(roc < minRoc) minRoc = roc;
		      if(roc > maxRoc) maxRoc = roc;
		    }

		}
	    }
	}
    }

  //
  if(minRoc > 7) AddEmptyEvent(EventNumber,2);  //we have an empty event in testboard2/tilted telescope
  if(maxRoc < 8) AddEmptyEvent(EventNumber,1);  //we have an empty event in testboard1/straight telescope 

}



void HistoManager::CollectHits(pixelTBTrackerHit* Hit, G4int EventNumber, std::map<int,pixelTBTrackerHit>& mapofHits)
{
  // let's manipulate all the pixels in the array that are hit
  // I assume that the array has been initialized to be all 0

  // we need the number of the chamber that was hit

  pixelTBTrackerHit originalHit = *Hit;
  G4double totalEnergy = originalHit.GetEdep();

  G4int chamberNo = originalHit.GetChamberNb();
  const G4int nCol = 52 * 3;
  //const G4int nRow = 80 * 3;
  //const G4int npixel = nCol * nRow;

  G4ThreeVector pos = originalHit.GetPos();
  G4int ROC = (int)(pos.getZ()/16.0 +0.5) - 11;
  G4int Col = (chamberNo%nCol)/3;
  G4int Row = (chamberNo/nCol)/3;
  
  G4double CSFraction = 1/3.;
  G4double CSFractionCorner = 1/9.;
 //let's save everything in a map

 if((chamberNo%nCol)%3 == 0) // this gives the COL
   {
     if((chamberNo/nCol)%3 == 1) // this gives the ROW
       {
	 //this is a left-center hit, we share charge to the pixel on the left
	 if(Col != 0)
	   {
	     //pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
	     //pixelArray[ROC][Col-1][Row] += CSFraction*Hit->GetEdep();

	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-CSFraction)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-CSFraction)*totalEnergy);
	       }
	     //same for the pixel we share the charge with

	     if(mapofHits.count(ROC*10000+(Col-1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFraction*totalEnergy);
		 Hit->SetChamberNb(chamberNo -1);
		 mapofHits[ROC*10000+(Col-1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col-1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col-1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
	       }

	   }
	 else
	   {
	     // pixelArray[ROC][Col][Row] += Hit->GetEdep();
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
	       }

	   }

       }
     if((chamberNo/nCol)%3 == 0) // this gives the ROW
       {
	 //this is a left-bottom hit, we have to share some charge to the neigboring pixels
	 if(Col != 0 && Row != 0)
	   {
	     //pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-5*CSFractionCorner)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-5*CSFractionCorner)*totalEnergy);
	       }
	     //chargeSharing Col-1
	     if(mapofHits.count(ROC*10000+(Col-1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo -1);
		 mapofHits[ROC*10000+(Col-1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col-1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col-1)*100+Row].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row-1
	     if(mapofHits.count(ROC*10000+Col*100+(Row-1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo - nCol );
		 mapofHits[ROC*10000+Col*100+(Row-1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row-1)].GetEdep();
		 mapofHits[ROC*10000+Col-1*100+(Row-1)].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row-1 and col-1
	     if(mapofHits.count(ROC*10000+(Col-1)*100+(Row-1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo - nCol -1 );
		 mapofHits[ROC*10000+(Col-1)*100+(Row-1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col-1)*100+(Row-1)].GetEdep();
		 mapofHits[ROC*10000+(Col-1)*100+(Row-1)].SetEdep(tempE + CSFractionCorner*totalEnergy);
	       }
	   }
	 else
	   {
	     if(Col == 0 && Row != 0)
	       {
		 // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-2*CSFraction)*totalEnergy);
		   }
		 //		 pixelArray[ROC][Col][Row-1] += CSFraction*Hit->GetEdep();
		 //chargeSharing row-1
		 if(mapofHits.count(ROC*10000+Col*100+(Row-1)) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo - nCol );
		     mapofHits[ROC*10000+Col*100+(Row-1)] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+(Row-1)].GetEdep();
		     mapofHits[ROC*10000+Col-1*100+(Row-1)].SetEdep(tempE + CSFraction*totalEnergy);
		   }
		 
	       }
	     else if(Col != 0 && Row == 0)
	       {
		 //		 pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-2*CSFraction)*totalEnergy);
		   }

		 //		 pixelArray[ROC][Col-1][Row] += CSFraction*Hit->GetEdep();
		 //chargeSharing Col-1
		 if(mapofHits.count(ROC*10000+(Col-1)*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo -1);
		     mapofHits[ROC*10000+(Col-1)*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+(Col-1)*100+Row].GetEdep();
		     mapofHits[ROC*10000+(Col-1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
		   }
	       }
	     else
	       {
		 //	 pixelArray[ROC][Col][Row] += Hit->GetEdep();
		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
		   }
	       }
	   }
       }

     if((chamberNo/nCol)%3 == 2) // this gives the ROW
       {
	 //this is a left-top hit, we have to share some charge to the neighboring pixel
	 if(Col != 0 && Row != 79)
	   {
	     //pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-5*CSFractionCorner)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-5*CSFractionCorner)*totalEnergy);
	       }
	     //chargeSharing Col-1
	     if(mapofHits.count(ROC*10000+(Col-1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo -1);
		 mapofHits[ROC*10000+(Col-1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col-1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col-1)*100+Row].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row+1
	     if(mapofHits.count(ROC*10000+Col*100+(Row+1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo + nCol );
		 mapofHits[ROC*10000+Col*100+(Row+1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row+1)].GetEdep();
		 mapofHits[ROC*10000+Col-1*100+(Row+1)].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row+1 and col-1
	     if(mapofHits.count(ROC*10000+(Col-1)*100+(Row+1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo + nCol -1 );
		 mapofHits[ROC*10000+(Col-1)*100+(Row+1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col-1)*100+(Row+1)].GetEdep();
		 mapofHits[ROC*10000+(Col-1)*100+(Row+1)].SetEdep(tempE + CSFractionCorner*totalEnergy);
	       }
	   }
	 else
	   {
	     if(Col == 0 && Row != 79)
	       {
		 // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 // pixelArray[ROC][Col][Row+1] += CSFraction*Hit->GetEdep();

		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing row+1
		 if(mapofHits.count(ROC*10000+Col*100+(Row+1)) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo + nCol );
		     mapofHits[ROC*10000+Col*100+(Row+1)] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+(Row+1)].GetEdep();
		     mapofHits[ROC*10000+Col*100+(Row+1)].SetEdep(tempE + CSFraction*totalEnergy);
		   }



	       }
	     else if(Col != 0 && Row == 79)
	       {
		 //pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 //pixelArray[ROC][Col-1][Row] += CSFraction*Hit->GetEdep();

		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing Col-1
		 if(mapofHits.count(ROC*10000+(Col-1)*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo -1);
		     mapofHits[ROC*10000+(Col-1)*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+(Col-1)*100+Row].GetEdep();
		     mapofHits[ROC*10000+(Col-1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
		   }
	  
	       }
	     else
	       {
		 //pixelArray[ROC][Col][Row] += Hit->GetEdep();

		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
		   }
	   
	       }
	   }
       }
   }
 else if((chamberNo%nCol)%3 == 1) // this gives the COL
   {
     if((chamberNo/nCol)%3 == 1) // this gives the ROW
       {
	 //this is a center hit, all the charge stays in one pixel
	 //pixelArray[ROC][Col][Row] += Hit->GetEdep();
	 //pixelHit
	 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	   {
	     //if the pixel is not hit yet, we add the hit to our collection
	     mapofHits[ROC*10000+Col*100+Row] = *Hit;
	   }
	 else
	   {
	     //the pixel was already hit, we just add the energy to the already hit pixel
	     G4double tempE;
	     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
	     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
	   }
	 
       }
     if((chamberNo/nCol)%3 == 0) // this gives the ROW
       {
	 //this is a bottom hit, we have to share some charge to the pixel underneath
	 if(Row != 0)
	   {
	     // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
	     //pixelArray[ROC][Col][Row - 1] += CSFraction*Hit->GetEdep();
	  
//pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
	       }
	     //chargeSharing row-1
	     if(mapofHits.count(ROC*10000+Col*100+(Row-1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFraction*totalEnergy);
		 Hit->SetChamberNb(chamberNo - nCol );
		 mapofHits[ROC*10000+Col*100+(Row-1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row-1)].GetEdep();
		 mapofHits[ROC*10000+Col*100+(Row-1)].SetEdep(tempE + CSFraction*totalEnergy);
	       }

	   }
	 else
	   {
	     //pixelArray[ROC][Col][Row] += Hit->GetEdep();
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
	       }
	   }
       }
     if((chamberNo/nCol)%3 == 2) // this gives the ROW
       {
	 //this is a top hit, we have to share some charge to the pixel above
	 if(Row != 79)
	   {
	     // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
	     //pixelArray[ROC][Col][Row + 1] += CSFraction*Hit->GetEdep();

	     //pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
	       }
	     //chargeSharing row+1
	     if(mapofHits.count(ROC*10000+Col*100+(Row+1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFraction*totalEnergy);
		 Hit->SetChamberNb(chamberNo + nCol );
		 mapofHits[ROC*10000+Col*100+(Row+1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row+1)].GetEdep();
		 mapofHits[ROC*10000+Col*100+(Row+1)].SetEdep(tempE + CSFraction*totalEnergy);
	       }
	   }
	 else
	   {
	     //pixelArray[ROC][Col][Row] += Hit->GetEdep();
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
	       }

	   }
       }
   }
 else if((chamberNo%nCol)%3 == 2) // this gives the COL
   {
     if((chamberNo/nCol)%3 == 1) // this gives the ROW
       {
	 //this is a right-center hit, we share charge to the pixel on the right
	 if(Col != 51)
	   {
	     // pixel hit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-CSFraction)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-CSFraction)*totalEnergy);
	       }
	     //same for the pixel we share the charge with

	     if(mapofHits.count(ROC*10000+(Col+1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFraction*totalEnergy);
		 Hit->SetChamberNb(chamberNo + 1);
		 mapofHits[ROC*10000+(Col+1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col+1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col+1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
	       }

	   }
	 else
	   {
	     //	     pixelArray[ROC][Col][Row] += Hit->GetEdep();

	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
	       }
	   }

       }
     if((chamberNo/nCol)%3 == 2) // this gives the ROW
       {
	 //this is a right-top hit, we have to share some charge to the neigboring pixels
	 if(Col != 51 && Row != 79)
	   {
	     //pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-5*CSFractionCorner)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-5*CSFractionCorner)*totalEnergy);
	       }
	     //chargeSharing Col+1
	     if(mapofHits.count(ROC*10000+(Col+1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo +1);
		 mapofHits[ROC*10000+(Col+1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col+1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col+1)*100+Row].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row+1
	     if(mapofHits.count(ROC*10000+Col*100+(Row+1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo + nCol );
		 mapofHits[ROC*10000+Col*100+(Row+1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row+1)].GetEdep();
		 mapofHits[ROC*10000+Col*100+(Row+1)].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row+1 and col+1
	     if(mapofHits.count(ROC*10000+(Col+1)*100+(Row+1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo + nCol + 1 );
		 mapofHits[ROC*10000+(Col+1)*100+(Row+1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col+1)*100+(Row+1)].GetEdep();
		 mapofHits[ROC*10000+(Col+1)*100+(Row+1)].SetEdep(tempE + CSFractionCorner*totalEnergy);
	       }

	   }
	 else
	   {
	     if(Col == 51 && Row != 79)
	       {
		 // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 // pixelArray[ROC][Col][Row+1] += CSFraction*Hit->GetEdep();
	       
		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing row+1
		 if(mapofHits.count(ROC*10000+Col*100+(Row+1)) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo + nCol );
		     mapofHits[ROC*10000+Col*100+(Row+1)] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+(Row+1)].GetEdep();
		     mapofHits[ROC*10000+Col*100+(Row+1)].SetEdep(tempE + CSFraction*totalEnergy);
		   }


	       }
	     else if(Col != 51 && Row == 79)
	       {
		 // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 // pixelArray[ROC][Col+1][Row] += CSFraction*Hit->GetEdep();
		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing Col+1
		 if(mapofHits.count(ROC*10000+(Col+1)*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo +1);
		     mapofHits[ROC*10000+(Col+1)*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+(Col+1)*100+Row].GetEdep();
		     mapofHits[ROC*10000+(Col+1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
		   }

	       }
	     else
	       {
		 //	 pixelArray[ROC][Col][Row] += Hit->GetEdep();
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
		   }
	       }
	   }
       }

     if((chamberNo/nCol)%3 == 0) // this gives the ROW
       {
	 //this is a right-bottom hit, we have to share some charge to the neighboring pixel
	 if(Col != 51 && Row != 0)
	   {
	     //pixelHit
	     if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep((1-5*CSFractionCorner)*totalEnergy);
		 mapofHits[ROC*10000+Col*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		 mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-5*CSFractionCorner)*totalEnergy);
	       }
	     //chargeSharing Col+1
	     if(mapofHits.count(ROC*10000+(Col+1)*100+Row) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo +1);
		 mapofHits[ROC*10000+(Col+1)*100+Row] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col+1)*100+Row].GetEdep();
		 mapofHits[ROC*10000+(Col+1)*100+Row].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing row-1
	     if(mapofHits.count(ROC*10000+Col*100+(Row-1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(2*CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo - nCol );
		 mapofHits[ROC*10000+Col*100+(Row-1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+Col*100+(Row-1)].GetEdep();
		 mapofHits[ROC*10000+Col*100+(Row-1)].SetEdep(tempE + 2*CSFractionCorner*totalEnergy);
	       }
	     //chargeSharing col+1 and row-1 
	     if(mapofHits.count(ROC*10000+(Col+1)*100+(Row-1)) == 0)
	       {
		 //if the pixel is not hit yet, we add the hit to our collection
		 Hit->SetEdep(CSFractionCorner*totalEnergy);
		 Hit->SetChamberNb(chamberNo - nCol + 1 );
		 mapofHits[ROC*10000+(Col+1)*100+(Row-1)] = *Hit;
	       }
	     else
	       {
		 //the pixel was already hit, we just add the energy to the already hit pixel
		 G4double tempE;
		 tempE = mapofHits[ROC*10000+(Col+1)*100+(Row-1)].GetEdep();
		 mapofHits[ROC*10000+(Col+1)*100+(Row-1)].SetEdep(tempE + CSFractionCorner*totalEnergy);
	       }
	     ///////
	   }
	 else
	   {
	     if(Col == 51 && Row != 0)
	       {
		 //pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 //pixelArray[ROC][Col][Row-1] += CSFraction*Hit->GetEdep();

		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing row-1
		 if(mapofHits.count(ROC*10000+Col*100+(Row-1)) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo - nCol );
		     mapofHits[ROC*10000+Col*100+(Row-1)] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+(Row-1)].GetEdep();
		     mapofHits[ROC*10000+Col*100+(Row-1)].SetEdep(tempE + CSFraction*totalEnergy);
		   }

	       }
	     else if(Col != 51 && Row == 0)
	       {
		 // pixelArray[ROC][Col][Row] += (1-CSFraction)*Hit->GetEdep();
		 //pixelArray[ROC][Col+1][Row] += CSFraction*Hit->GetEdep();

		 //pixelHit
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep((1-1*CSFraction)*totalEnergy);
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + (1-1*CSFraction)*totalEnergy);
		   }
		 //chargeSharing Col+1
		 if(mapofHits.count(ROC*10000+(Col+1)*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     Hit->SetEdep(CSFraction*totalEnergy);
		     Hit->SetChamberNb(chamberNo +1);
		     mapofHits[ROC*10000+(Col+1)*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+(Col+1)*100+Row].GetEdep();
		     mapofHits[ROC*10000+(Col+1)*100+Row].SetEdep(tempE + CSFraction*totalEnergy);
		   }


	       }
	     else
	       {
		 //pixelArray[ROC][Col][Row] += Hit->GetEdep();
		 if(mapofHits.count(ROC*10000+Col*100+Row) == 0)
		   {
		     //if the pixel is not hit yet, we add the hit to our collection
		     mapofHits[ROC*10000+Col*100+Row] = *Hit;
		   }
		 else
		   {
		     //the pixel was already hit, we just add the energy to the already hit pixel
		     G4double tempE;
		     tempE = mapofHits[ROC*10000+Col*100+Row].GetEdep();
		     mapofHits[ROC*10000+Col*100+Row].SetEdep(tempE + totalEnergy);
		   }

	       }
	   }
       }
   }

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


void HistoManager::SetROCVectors(G4ThreeVector orROC[16], G4ThreeVector upleftROC[16], G4ThreeVector lowrightROC[16])
{
  //let's calculate the normed vectors

  for(int ii = 0; ii < 16; ++ii)
    {
      origROC[ii] = orROC[ii];    //the vector of the origin does not need to be modified

      // u is the vector pointing from row 0 to row 80
      G4ThreeVector U = upleftROC[ii] - orROC[ii]; 
      // v is the vector pointing from col 0 to col 52
      G4ThreeVector V = lowrightROC[ii] - orROC[ii]; 

      //now we have to norm these vectors 

      UnormROC[ii] = U.unit();
      //      UnormROC[ii] = UnormROC[ii] * 0.1*CLHEP::mm;
      VnormROC[ii] = V.unit();
      //VnormROC[ii] = VnormROC[ii] * 0.15*CLHEP::mm;

      NnormROC[ii] =  UnormROC[ii].cross(VnormROC[ii]);


      G4cout <<"ROC " << ii << " Unorm" << UnormROC[ii] <<" Mag " << UnormROC[ii].mag() << G4endl;
      G4cout <<"ROC " << ii << " Vnorm" << VnormROC[ii] <<" Mag " << VnormROC[ii].mag() << G4endl;
      G4cout <<"ROC " << ii << " Nnorm" << NnormROC[ii] <<" Mag " << NnormROC[ii].mag() << G4endl;

      G4cout << "u x v is orthogonal? " <<  ((UnormROC[ii].isOrthogonal(VnormROC[ii])) ? "true" : "false") << G4endl;

    }

}

G4double HistoManager::CheckHitDistance(pixelTBTrackerHit* Hit)
{
  G4ThreeVector pos = Hit->GetPos();
  G4int ROC = (int) (z/16.0 + 0.5) - 11;

  G4double distance = (NnormROC[ROC].dot(pos - origROC[ROC])) / NnormROC[ROC].mag(); 
  return distance;

}

void HistoManager::CheckPixel(pixelTBTrackerHit* Hit)
{
  // let's check what pixel the hit is associated to by chamber number and by the position
  G4int chamberNo = Hit->GetChamberNb();
  const G4int nCol = 52;
  //const G4int nRow = 80 * 3;
  //const G4int npixel = nCol * nRow;
  
  G4ThreeVector pos = Hit->GetPos();
  G4int ROC = (int)(pos.getZ()/16.0 +0.5) - 11;
  G4int Col = (chamberNo%nCol);
  G4int Row = (chamberNo/nCol);
  

  //now let's calculate the position of the hit projected on the plane

  G4ThreeVector onPlaneHit = pos - CheckHitDistance(Hit)*NnormROC[ROC];
  //  G4ThreeVector onPlaneVector = origROC[ROC] - onPlaneHit;
  G4ThreeVector onPlaneVector = onPlaneHit - origROC[ROC];

  G4double planeV = onPlaneVector.dot(VnormROC[ROC]) / (VnormROC[ROC].mag2());
  G4double planeU = onPlaneVector.dot(UnormROC[ROC]) / (UnormROC[ROC].mag2());

  G4int onPlaneCol = planeV / (0.15*CLHEP::mm);
  G4int onPlaneRow = planeU / (0.10*CLHEP::mm);

  G4cout << "COL>>> from Chamber(" << Col << ") from coord(" << onPlaneCol << ")" ;
  (Col != onPlaneCol)? G4cout << "!!!!!!!!!!" << G4endl : G4cout << G4endl;
  G4cout << "ROW>>> from Chamber(" << Row << ") from coord(" << onPlaneRow << ")" ;
  (Row != onPlaneRow)? G4cout << "!!!!!!!!!!" << G4endl : G4cout << G4endl;

  ///let's say we only spread charge to the next pixels in col and row direction, not diagonally
  //in Col we spread the charge among pixel

  G4int onPlaneColA = (planeV - 0.5 * 0.15*CLHEP::mm)/(0.15*CLHEP::mm);
  G4int onPlaneColB = (planeV + 0.5 * 0.15*CLHEP::mm)/(0.15*CLHEP::mm);

  G4cout << "Charge spread between Col " << onPlaneColA << " and Col " << onPlaneColB << G4endl;
  G4double centerColA = 0.15*CLHEP::mm * onPlaneColA + 0.075*CLHEP::mm;
  G4double centerColB = 0.15*CLHEP::mm * onPlaneColB + 0.075*CLHEP::mm;


}


G4ThreeVector HistoManager::projectOnSensorPlane(pixelTBTrackerHit* Hit){

  G4ThreeVector pos = Hit->GetPos();
  G4ThreeVector dir = Hit->GetDir();
  G4int ROC = (int)(pos.getZ()/16.0 +0.5) - 11;

  G4double d;

  d = (origROC[ROC] - pos).dot(NnormROC[ROC]) / (dir.dot(NnormROC[ROC]));

  return (d*dir + pos);

  //  return  pos - CheckHitDistance(Hit)*NnormROC[ROC];

}





void HistoManager::getHitmap(pixelTBTrackerHit* Hit)
{

  // we only look at a single vertical plane, let's say ROC 8

  G4int testROC = 8;
  //the ROC the hit is on is
  G4ThreeVector pos = Hit->GetPos();
  G4int ROC = (int)(pos.getZ()/16.0 +0.5) - 11;

  G4int posx = (int)(pos.getX() * 1000);   //in um
  G4int posy = (int)(pos.getY() * 1000);   //in um


  if( ROC == testROC)
    {
      //      G4cout << "Filling hitmap" << G4endl;
      reducedHitmap->Fill(posx%300,posy%200);
    }

}

