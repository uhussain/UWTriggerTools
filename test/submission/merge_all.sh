#!/bin/bash

# Submit all UCT ntuple jobs
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 DATE"
  exit 1
fi

#OUTPUTDIR=/scratch/$LOGNAME
OUTPUTDIR=../data/LSB50/nohcal

hadd -f $OUTPUTDIR/uct_rates_stage1b.root $hdfs/2013-04-23-Stage1c-Norm-makeRateTrees_cfg/*root &
hadd -f $OUTPUTDIR/uct_rates_stage1c.root $hdfs/2013-04-23-Stage1c-HCAL-makeRateTrees_cfg/*root &

wait
