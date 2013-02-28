#!/bin/bash
voms-proxy-init --voms cms --valid 100:00

# Submit all UCT ntuple jobs
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 DATE"
  exit 1
fi

./submit_rates.sh $1-EIC3-Rates 3
./submit_eg_efficiency.sh $1-EIC3-EGEfficiency 3
./submit_tau_efficiency.sh $1-TauEfficiency
./submit_jet_efficiency.sh $1-JetEfficiency
