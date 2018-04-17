#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/MCFlat

# integrated luminosity for data
LUMI=2215

root -l -q selectMC.C+\(\"samplesMC.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d
