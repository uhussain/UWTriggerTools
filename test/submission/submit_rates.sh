#!/bin/bash

# Submit rate ntuple jobs on the ZeroBias3 dataset (2012C)
EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME EICTHRESH"
  exit 1
fi

farmoutAnalysisJobs $1 \
  --infer-cmssw-path \
  --input-file-list=zero_bias_files.txt \
  ../makeRateTrees_cfg.py  stage1B=1 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' \
  eicIsolationThreshold=$2
