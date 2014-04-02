'''

Module for computing RECO level total sumET and sumHT
from PF objects.

Produces two vectors of LeafCandidates (1 object in each vector), with labels
"set" and "sht".

Author: Evan K. Friis, UW Madison

'''

import FWCore.ParameterSet.Config as cms

pfSumET = cms.EDProducer(
    "PFSumET",
    src = cms.InputTag("particleFlow"),
    maxEta = cms.double(3),
    minForHt = cms.double(5),
    excludeMuons = cms.bool(True)
)
