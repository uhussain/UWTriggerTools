#!/bin/bash

# Submit all UCT ntuple jobs
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 DATE"
  exit 1
fi

#OUTPUTDIR=/scratch/$LOGNAME
OUTPUTDIR=../data/LSB50

hadd -f $OUTPUTDIR/uct_mc_efficiency_pu50.root $hdfs/${1}-PU50-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_mc_efficiency_pu35.root $hdfs/${1}-PU35-makeEfficiencyTree_cfg/*root &

wait
