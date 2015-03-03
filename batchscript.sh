#!/bin/bash

#. $G4WORKDIR/bin/Darwin-clang/pixelTB vis.mac conf/conf0.txt
/afs/cern.ch/user/a/akornmay/geant4_workdir/bin/Linux-g++/pixelTB novis.mac conf/conf0.txt &
/afs/cern.ch/user/a/akornmay/geant4_workdir/bin/Linux-g++/pixelTB novis.mac conf/conf1.txt &
/afs/cern.ch/user/a/akornmay/geant4_workdir/bin/Linux-g++/pixelTB novis.mac conf/conf2.txt &
/afs/cern.ch/user/a/akornmay/geant4_workdir/bin/Linux-g++/pixelTB novis.mac conf/conf3.txt &


