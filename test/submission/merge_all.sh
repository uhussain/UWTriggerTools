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

hadd -f $OUTPUTDIR/uct_eg_efficiency_eic3.root $hdfs/${1}-EIC3-EGEfficiency-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_eg_efficiency_eic4.root $hdfs/${1}-EIC4-EGEfficiency-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_tau_efficiency.root $hdfs/${1}-TauEfficiency-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_jet_efficiency.root $hdfs/${1}-JetEfficiency-makeEfficiencyTree_cfg/*root &
hadd -f $OUTPUTDIR/uct_rates_eic3.root $hdfs/${1}-EIC3-Rates-makeRateTrees_cfg/*root &
hadd -f $OUTPUTDIR/uct_rates_eic4.root $hdfs/${1}-EIC4-Rates-makeRateTrees_cfg/*root &

wait
