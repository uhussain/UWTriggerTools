import FWCore.ParameterSet.Config as cms

uct2015L1ExtraParticles = cms.EDProducer(
    "UCT2015L1ExtraProducer",

    produceMuonParticles = cms.bool(False),
    muonSource = cms.InputTag("gtDigis"),

    produceCaloParticles = cms.bool(True),
    egLabel   = cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked"),
    tauLabel  = cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked"),
    jetLabel  = cms.InputTag("UCT2015Producer", "JetUnpacked"),

    etMissLabel = cms.InputTag("UCT2015Producer", "METUnpacked"),
    etTotLabel  = cms.InputTag("UCT2015Producer", "SETUnpacked"),
    htMissLabel = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    htTotLabel  = cms.InputTag("UCT2015Producer", "SHTUnpacked"),

    egIso = cms.double(0.2),
    tauIso = cms.double(0.2),

    regionLabel   = cms.InputTag("uctDigis"),
    rgnThreshold  = cms.double(0.),
    regionLSB     = cms.double(0.5),
    muFixedIso    = cms.double(5.),
    muRelIso      = cms.double(-999.),
    
    centralBxOnly = cms.bool(True)
)


