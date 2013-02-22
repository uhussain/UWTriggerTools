#!/bin/bash

# Submit efficiency ntuple jobs on the Zee skim.
EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME EIC_THRESHOLD"
  exit 1
fi

farmoutAnalysisJobs $1 \
  --infer-cmssw-path \
  --input-file-list=zee_skim_files.txt \
  ../makeEfficiencyTree_cfg.py  stage1B=1\
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' \
  eicIsolationThreshold=$2
