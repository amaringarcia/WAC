#!/bin/bash

if [ $# -gt 2 ]; then
  echo "usage: batchRunStatUncertain basedirectory productiondirectory"
  exit 1
fi

if [ $# -lt 2 ]; then
  echo "usage: batchRunStatUncertain basedirectory productiondirectory"
  exit 1
fi

BASEDIRECTORY=$1
PRODUCTIONDIRECTORY=$2

if [ ! -d $BASEDIRECTORY/$PRODUCTIONDIRECTORY ]
then
  echo Production directory $PRODUCTIONDIRECTORY does NOT exist. WE CANNOT PROCEED!!!!
  exit 1	
fi

# the configuration file
CONFIGURATIONFILE=$BASEDIRECTORY/$PRODUCTIONDIRECTORY/configuration.json

# extract needed information
RAPIDITIES=(`sed -n '/"abs\_y"\s*:\s*\[\(.*\)\],/p' configuration.json | sed 's/\s*"abs\_y"\s*:\s*\[\(.*\)\],/\1/' | sed 's/,/0/g' | sed 's/0\.//g'`)
EVENTFILTERS=(`sed -n '/"teventfilter"\s*:\s*\[\(.*\)\],/p' configuration.json | sed 's/\s*"teventfilter"\s*:\s*\[\(.*\)\],/\1/' | tr ',' ' '`)
NRAPIDITIES=${#RAPIDITIES[@]}
NEVENTFILTERS=${#EVENTFILTERS[@]}

# the results production tag
PRODUCTIONTIME=`date +%Y%m%d_%H%M%S`

# submit the extraction of results with statistical uncertainties
# one job per rapidity (parallel) and per multiplicity class (sequential)
for (( i=0; i<NRAPIDITIES; i++ ))
do
  PRODUCTIONTAG=`printf "%s_Rap%03d" ${PRODUCTIONTIME} ${RAPIDITIES[i]}`
  echo Submitting for production tag ${PRODUCTIONTAG}
  JOBID=
  for (( j=0; j<NEVENTFILTERS; j++ ))
  do
    if [ ${j} -eq 0 ]
    then
      DEPENDENCY=
    else
      DEPENDENCY="-d afterany:${JOBID}"
    fi
    echo Submitting multiplicity class ${EVENTFILTERS[j]}
    JOBNAME=`printf "waitStatsUncertain_%03d_%03d" ${i} ${j}`
    cmd="sbatch -J ${JOBNAME} ${DEPENDENCY} --workdir=${BASEDIRECTORY}/${PRODUCTIONDIRECTORY} --mem-per-cpu=8000 --time=07:00:00 -o ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.out -e ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runScriptInSingularity.sh /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runStatsUncertain.sh ${PRODUCTIONTAG} ${i} ${j}"
    JOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
    echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
    echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
  done
done


