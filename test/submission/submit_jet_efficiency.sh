#!/bin/bash

# Submit efficiency ntuple jobs on the mu+jet skim
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

farmoutAnalysisJobs $1 \
  --infer-cmssw-path \
  --input-file-list=jet_skim_files.txt \
  ../makeEfficiencyTree_cfg.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
