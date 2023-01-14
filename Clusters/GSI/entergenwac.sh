#!/bin/bash
export ALIEN_SITE=GSI
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@Python::v3.9.12-10)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@AliGenerators::v20221125-1)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@ROOT::v6-26-04-patches-alice2-22)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@CMake::v3.23.1-12)
source set-WAC-GSI

