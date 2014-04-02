#!/bin/bash

# Submit rate ntuple jobs on the ZeroBias3 dataset (2012C)
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

farmoutAnalysisJobs $1-EVAN \
  --infer-cmssw-path \
  --input-file-list=NeutrinoGun13_Pu4025ns.txt \
  --input-files-per-job=5 \
  ../makeRateTrees_cfg.py  isMC=1  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' \
  eicIsolationThreshold=3

farmoutAnalysisJobs $1-DPG \
  --infer-cmssw-path \
  --input-file-list=NeutrinoGun13_Pu4025ns.txt  \
  --input-files-per-job=5 \
  ../runUCTToGTinterface.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 
