#!/bin/bash

# setting the root and pythia scenario
export ALIEN_SITE=GSI
source /cvmfs/alice.cern.ch/etc/login.sh
LATEST="VO_ALICE@AliGenerators::v20221125-1"
export ALIPHYSICS_VERSION=$LATEST

eval $(alienv printenv $LATEST) 
echo $LATEST

export PYTHIA8=/cvmfs/alice.cern.ch/el7-x86_64/Packages/pythia/v8304-44

####################################################################################################
echo "Setting up WAC"
####################################################################################################
export WAC_ROOT=/lustre/alice/users/$USER/CLUSTERMODELWAC
export WAC_SOURCE="$WAC_ROOT"
export WAC_BIN="$WAC_ROOT/bin"
export WAC_LIB="$WAC_ROOT/lib"

export PATH="$WAC_BIN:$PATH"
export DYLD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$DYLD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$LD_LIBRARY_PATH"

TASKIX=$SLURM_ARRAY_TASK_ID
SEED=$(( (SLURM_ARRAY_TASK_ID + SLURM_ARRAY_JOB_ID*1000) % 900000000 ))
echo "The seed is $SEED"

# Execute application code
RunPythiaSimulationTwoParticlesDiff $TASKIX $SEED

