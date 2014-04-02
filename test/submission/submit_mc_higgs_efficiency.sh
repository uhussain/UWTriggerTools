#!/bin/bash

# Submit efficiency ntuple jobs on the mu-tau skim.
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

farmoutAnalysisJobs $1-PU50 \
  --infer-cmssw-path \
  --input-file-list=htt_14TeV_mc_pu50.txt \
  ../makeEfficiencyTree_cfg.py  stage1B=1 isMC=1 stage1=0 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

farmoutAnalysisJobs $1-PU35 \
  --infer-cmssw-path \
  --input-file-list=htt_14TeV_mc_pu35.txt \
  ../makeEfficiencyTree_cfg.py  stage1B=1 isMC=1 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

farmoutAnalysisJobs $1-PU50-HCALONLY \
  --infer-cmssw-path \
  --input-file-list=htt_14TeV_mc_pu50.txt \
  ../makeEfficiencyTree_cfg.py  eicCardHcalOnly=1 stage1B=1 stage1=0 isMC=1 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

farmoutAnalysisJobs $1-PU35-HCALONLY \
  --infer-cmssw-path \
  --input-file-list=htt_14TeV_mc_pu35.txt \
  ../makeEfficiencyTree_cfg.py  eicCardHcalOnly=1 stage1B=1 stage1=0 isMC=1 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 
