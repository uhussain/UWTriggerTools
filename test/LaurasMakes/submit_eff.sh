#!/bin/bash

# Submit rate ntuple jobs on the ZeroBias3 dataset (2012C)
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

farmoutAnalysisJobs $1-WithHF \
  --infer-cmssw-path \
  --input-file-list=jet_skim_files.txt \
  --input-files-per-job=1 \
  makeEfficiencyTree_Jets_cfg.py isMC=0 \
    'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

#farmoutAnalysisJobs $1-NoHF \
#  --infer-cmssw-path \
#  --input-file-list=jet_skim_files.txt \
#  --input-files-per-job=1 \
#  makePuMultTrees_cfg.py isMC=0 inclHF=0 \
#    'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 
