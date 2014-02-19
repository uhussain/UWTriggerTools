#!/bin/bash

# Submit rate ntuple jobs on the ZeroBias3 dataset (2012C)
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

#farmoutAnalysisJobs $1-HCAL \
#  --infer-cmssw-path \
#  --input-file-list=zero_bias_files.txt \
#  --input-files-per-job=5 \
#  --job-count=10 \
#  ../makeRateTrees_cfg.py  isMC=0 stage1B=0 stage1=1 eicCardHcalOnly=1 \
#  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' \
#  eicIsolationThreshold=3

farmoutAnalysisJobs $1-Norm \
  --infer-cmssw-path \
  --input-file-list=zero_bias_files.txt \
  --input-files-per-job=1 \
  --job-count=150 \
  makeRateTrees_Jets_cfg.py  isMC=0 stage1B=0 stage1=1 eicCardHcalOnly=0 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' \
  eicIsolationThreshold=3
