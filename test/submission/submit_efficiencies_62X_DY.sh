#!/bin/bash

# Submit efficiency ntuple jobs on the mu-tau skim.
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

#farmoutAnalysisJobs $1-20PU25ns \
#  --infer-cmssw-path \
#  --input-file-list=DY20bs25.txt \
#  ../makeEfficiencyTree_cfg.py  stage1B=0 isMC=1 \
#  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

farmoutAnalysisJobs $1-DY-40PU25ns \
 --infer-cmssw-path \
 --input-file-list=DY40bs25.txt \
 ../makeEfficiencyTree_cfg.py  isMC=1 \
 'inputFiles=$inputFileNames' 'outputFile=$outputFileName'  
# --resubmit-failed-jobs




farmoutAnalysisJobs $1-TT-40PU25ns \
  --infer-cmssw-path \
  --input-file-list=TTbar25ns.txt \
  ../makeEfficiencyTree_cfg.py  isMC=1 \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'  
# --resubmit-failed-jobs

farmoutAnalysisJobs $1-DYGEN-40PU25ns \
  --infer-cmssw-path \
  --input-file-list=DY25ns40PUGENRAW.txt \
  ../makeGenEffi_cfg.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
