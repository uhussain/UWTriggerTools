#flake8: noqa
'''

Emulate the L1 and UCT upgrade primitives, and put them in the event.

Authors: Isobel Ojalvo, Sridhara Dasu (kludger)

'''

import FWCore.ParameterSet.Config as cms


from Configuration.StandardSequences.RawToDigi_Data_cff import *
from L1Trigger.UCT2015.Lut import *

# Modify the HCAL TPGs according to the proposed HTR modification.  If the HCAL
# is above a given energy threshold, set the MIP bit.
hackHCALMIPs = cms.EDProducer(
    "HcalTpgMipEmbedder",
    src = cms.InputTag("hcalDigis"),
    threshold = cms.double(3), # In GeV
    rawThreshold = cms.uint32(3), # In TPG rank
    cutOnRawBits = cms.bool(False), # What to cut on
)

uctDigis = cms.EDProducer(
    "L1RCTProducer",
    #hcalDigis = cms.VInputTag(cms.InputTag("hcalDigis")),
    hcalDigis = cms.VInputTag(cms.InputTag("hackHCALMIPs")),
    useEcal = cms.bool(True),
    useHcal = cms.bool(True),
    ecalDigis = cms.VInputTag(cms.InputTag("ecalDigis:EcalTriggerPrimitives")),
    BunchCrossings = cms.vint32(0),
    getFedsFromOmds = cms.bool(False),
    queryDelayInLS = cms.uint32(10),
    queryIntervalInLS = cms.uint32(100)#,
)

UCT2015EClusterProducer = cms.EDProducer(
    "UCT2015EClusterProducer",
    debug = cms.bool(False),
    puCorrect = cms.bool(True),
    puETMax = cms.uint32(7),
    eClusterSeed = cms.uint32(5),
    # Transparency correction calibration
    # calib_v4 = 2012 data.  Use calib_v1 (ideal) for MC.
    ecalCalibration = cms.vdouble(eg_calib_v4),
    ecalLSB = cms.double(0.5),
    ecalDigis = cms.VInputTag(cms.InputTag("ecalDigis:EcalTriggerPrimitives"))
)

UCT2015Producer = cms.EDProducer(
    "UCT2015Producer",
    puCorrect = cms.bool(True),
    useUICrho = cms.bool(True),
    useHI = cms.bool(False),
    # All of these uint32 thresholds are in GeV.
    puETMax = cms.uint32(7),
    regionETCutForHT = cms.uint32(5),
    regionETCutForMET = cms.uint32(0),
    minGctEtaForSums = cms.uint32(4),
    maxGctEtaForSums = cms.uint32(17),
    jetSeed = cms.uint32(5),
    egtSeed = cms.uint32(5),
    relativeIsolationCut = cms.double(1.0),
    relativeJetIsolationCut = cms.double(1.0),
    egammaLSB = cms.double(1.0), # This has to correspond with the value from L1CaloEmThresholds
    regionLSB = RCTConfigProducers.jetMETLSB,
)

UCTStage1BProducer = cms.EDProducer(
    "UCTStage1BProducer",
    puCorrect = cms.bool(True),
    egSeed = cms.uint32(5),
    tauSeed = cms.uint32(5),
    regionalHoECut = cms.double(0.05),
    egRelativeRgnIsolationCut = cms.double(0.1),
    egRelativeJetIsolationCut = cms.double(0.1),
    egRelativeEMRgnIsolationCut = cms.double(0.1),
    egRelativeEMJetIsolationCut = cms.double(0.1),
    tauRelativeRgnIsolationCut = cms.double(0.1),
    tauRelativeJetIsolationCut = cms.double(0.1),
    tauRelativeEMRgnIsolationCut = cms.double(0.1),
    tauRelativeEMJetIsolationCut = cms.double(0.1),
    egLSB = cms.double(0.5),
    tauLSB = cms.double(1.0),  # This has to correspond with the value from L1CaloEmThresholds
    regionLSB = RCTConfigProducers.jetMETLSB
)

from L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi import *

uctDigiStep = cms.Sequence(
    # Only do the digitization of objects that we care about
    #RawToDigi
    gctDigis
    * gtDigis
    * ecalDigis
    * hcalDigis
)

uctEmulatorStep = cms.Sequence(
    hackHCALMIPs
    # Now make UCT and L1 objects
    * uctDigis
    * UCT2015EClusterProducer
    * UCT2015Producer
    * UCTStage1BProducer
    * l1extraParticles
)

emulationSequence = cms.Sequence(uctDigiStep * uctEmulatorStep)
