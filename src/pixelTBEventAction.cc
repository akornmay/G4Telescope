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
// $Id: pixelTBEventAction.cc,v 1.11 2006-06-29 17:48:05 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "pixelTBEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "pixelTBTrackerSD.hh"
#include "pixelTBTrackerHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
pixelTBEventAction::pixelTBEventAction(HistoManager* histManager)
  :fHistManager(histManager)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
pixelTBEventAction::~pixelTBEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void pixelTBEventAction::BeginOfEventAction(const G4Event*)
{
  fHistManager->GetNumbersOfParticlesInEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void pixelTBEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  //
  if (event_id < 100 || event_id%100 == 0) 
    {
     G4cout << ">>> Event " << evt->GetEventID() << G4endl;
     //     G4cout << "    " << n_trajectories 
     //	   << " trajectories stored in this event." << G4endl;
    }

  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID("trackerCollection"); 
  //  HCE->AddHitsCollection( HCID, trackerCollection ); 

  G4HCofThisEvent* HCOTE = evt->GetHCofThisEvent();

  pixelTBTrackerHitsCollection* trackerCollection = 0;
  
  if(HCOTE)
    {
      trackerCollection = (pixelTBTrackerHitsCollection*)(HCOTE->GetHC(HCID));
    }
  if(trackerCollection)
    {
      G4int NbHits = trackerCollection->entries();
      //G4cout << "XXXX There are " << NbHits << " in the collection" <<  G4endl;

      if(NbHits == 0)
	{
	  fHistManager->AddEmptyEvent(evt->GetEventID(),1);
	  fHistManager->AddEmptyEvent(evt->GetEventID(),2);

	}

      if(NbHits != 0)
	{
	  //initialize a empty array
	  G4double pixelArray[16][52][80] = {{{0}}};

	  for(G4int i = 0; i < NbHits; i++)
	    {
	      fHistManager->CollectHits((*trackerCollection)[i], evt->GetEventID(),pixelArray);
	    }

	  fHistManager->AddHits(pixelArray, evt->GetEventID());
	}

    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
