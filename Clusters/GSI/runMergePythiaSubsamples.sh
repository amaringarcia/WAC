#!/bin/bash

if [ $# -gt 1 ]; then
  echo "usage: runMergePythiaSubsamples filename"
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "usage: runMergePythiaSubsamples filename"
  exit 1
fi

FILENAME=$1

# setting the root and pythia scenario
export ALIEN_SITE=GSI
source /cvmfs/alice.cern.ch/etc/login.sh
LATEST="VO_ALICE@ROOT::v6-26-04-patches-alice2-22"
export ALIPHYSICS_VERSION=$LATEST

eval $(alienv printenv $LATEST) 
echo $LATEST

# merge the subsamples in a single file
hadd ${FILENAME}.root BUNCH??/Output/${FILENAME}.root

