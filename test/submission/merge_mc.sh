#!/bin/bash

# Submit all UCT ntuple jobs
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 DATE"
  exit 1
fi

#OUTPUTDIR=/scratch/$LOGNAME
OUTPUTDIR=../data/LSB50/Apr29

hadd -f $OUTPUTDIR/uct_mc_efficiency_pu50_hcal.root $hdfs/${1}-PU50-HCALONLY-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_mc_efficiency_pu50_norm.root $hdfs/${1}-PU50-makeEfficiencyTree_cfg/*root &

wait
