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
// $Id: exampleN02.cc,v 1.16 2009-10-30 14:59:59 allison Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
// 
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "pixelTBDetectorConstruction.hh"
#include "pixelTBPhysicsList.hh"
#include "pixelTBPrimaryGeneratorAction.hh"
#include "pixelTBRunAction.hh"
#include "pixelTBEventAction.hh"
#include "pixelTBSteppingAction.hh"
#include "pixelTBSteppingVerbose.hh"
#include "HistoManager.hh"


#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"



void ReadSettings(char* confFile, HistoManager* histo)
{
  //  G4cout << "Reading configuration from file" << G4endl;

  if(confFile[0]=='\0')
    {
      G4cout << "Using default parameters" << G4endl;
      
      
      histo->SetFirstTurn(0);

    }
  else
    {
      std::ifstream istr(confFile,std::ios::in);
      if(!istr)
	{
	  G4cout << "Error: File \"" << confFile << "\" doesn't exist."  << G4endl;
      	}
      char buf[255];
      G4String Parameter;
      G4String Value;

      char equal;

      while(!istr.eof() && !istr.fail())
	{
	  istr>>Parameter;
	  if(Parameter[0]=='#')
	    {
	      istr.getline(buf,255);
	      continue;
	    }
	  istr >> equal >> Value;
	  istr.getline(buf,255);

	  if(equal!='=')
	    {
	      G4cout << "Error: Syntax error for parameter " << Parameter << G4endl;
	      continue;
	    }

	  if(Parameter=="FIRST_TURN")
	    {
	      histo->SetFirstTurn(atol(Value.c_str()));
	      continue;
	    }

	  if(Parameter=="LAST_TURN")
	    {
	      continue;
	    }
	}

    }

  G4cout << "RUNNING simulation from turn " << histo->GetFirstTurn() << G4endl;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new pixelTBSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // set an HistoManager
  //
  HistoManager*  histo = new HistoManager();

  // read config file if provided
  //
  char confName[255];
  if(argc<3)
    {
      confName[0]='\0';
    }
  else
    {
      strcpy(confName,argv[2]);
    }
  ReadSettings(confName, histo);

  // Run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  pixelTBDetectorConstruction* detector = new pixelTBDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new pixelTBPhysicsList;
  runManager->SetUserInitialization(physics);
   
  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new pixelTBPrimaryGeneratorAction(detector,histo);
  runManager->SetUserAction(gen_action);
  //
  G4UserRunAction* run_action = new pixelTBRunAction(histo);
  runManager->SetUserAction(run_action);
  //
  G4UserEventAction* event_action = new pixelTBEventAction(histo);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new pixelTBSteppingAction;
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    
     
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

