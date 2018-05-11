#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/work/x/xuyan/scratch/CMSSW_8_0_27/src/Ntuples

# integrated luminosity for data
LUMI=2215

root -l -q selectDATA.C+\(\"samplesDATA.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d *.pcm
