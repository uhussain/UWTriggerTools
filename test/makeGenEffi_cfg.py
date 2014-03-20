'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("ReRunningL1")

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
#"/store/mc/Fall13dr/QCD_Pt-10to15_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v3/00000/163B68E6-959E-E311-912A-003048C692D8.root",
"/store/user/ldodd/TT_Tune4C_13TeV-pythia8-tauola/TT_Tune4C_13TeV-pythia8-tauola-tsg_PU40bx25_POSTLS162_V2-v1/fb508503c16d6e4b02bc25104d11f7c2/skim_109_1_e51.root"

                             )   
                             )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
#process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")

process.l1ExtraTreeProducerUCT = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticlesUCT:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticlesUCT:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticlesUCT"),
   metLabel = cms.untracked.InputTag("l1extraParticlesUCT:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticlesUCT:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticlesUCT"),
   maxL1Extra = cms.uint32(20)
)

process.l1ExtraTreeProducer = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticles:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticles:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticles:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticles:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticles:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticles"),
   metLabel = cms.untracked.InputTag("l1extraParticles:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticles:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticles"),
   maxL1Extra = cms.uint32(20)
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('genEffis.root')
)

process.createGenParticlesEle =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(0),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(100)
)


process.createGenParticlesEleIso =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(0),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(1)
)

common_ntuple_branches = cms.PSet(
    # Run, lumi, event number
    run = cms.string("id.run"),
    lumi = cms.string("id.luminosityBlock"),
    evt = cms.string("id.event"),

    recoPt = cms.string("reco.pt"),
    recoEta = cms.string("reco.eta"),
    recoPhi = cms.string("reco.phi"),

    # Whether there exists a L1/UCT object matched to reco
    l1Match = cms.string("l1Match"),
    l1gMatch = cms.string("l1gMatch"),

    l1Pt = cms.string("? l1Match ? l1.pt : 0"),
    l1Eta = cms.string("? l1Match ? l1.eta : 0"),
    l1Phi = cms.string("? l1Match ? l1.phi : 0"),
    l1Type = cms.string("? l1Match ? l1.type() : -1"),
    # TODO add L1extra eta/phi indices

    l1DPhi = cms.string("? l1Match ? deltaPhi(l1.phi, reco.phi) : -1"),
    l1DR = cms.string("? l1Match ? deltaR(l1.eta, l1.phi, reco.eta, reco.phi) : -1"),

    l1gPt = cms.string("? l1gMatch ? l1g.pt : 0"),
    l1gEta = cms.string("? l1gMatch ? l1g.eta : 0"),
    l1gPhi = cms.string("? l1gMatch ? l1g.phi : 0"),

    l1gDPhi = cms.string("? l1gMatch ? deltaPhi(l1g.phi, reco.phi) : -1"),
    l1gDEta = cms.string("? l1gMatch ? l1g.eta - reco.eta : -10"),
    l1gDR = cms.string("? l1gMatch ? deltaR(l1g.eta, l1g.phi, reco.eta, reco.phi) : -1"),
)

# Specific to EG tau objects
egtau_branches = cms.PSet(
    l1gSecondRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedSecondRegionEt', -4) : -2"),
    l1gThirdRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedThirdRegionEt', -4) : -2"),
    l1gJetPt = cms.string("? l1gMatch ? l1g.getFloat('associatedJetPt', -4) : -2"),
    l1gEllIso = cms.string("? l1gMatch ? l1g.getInt('ellIsolation', -4) : -2"),
    l1gTauVeto = cms.string("? l1gMatch ? l1g.getInt('tauVeto', -4) : -2"),
    l1gMIP = cms.string("? l1gMatch ? l1g.getInt('mipBit', -4) : -2"),
    l1gIsEle = cms.string("? l1gMatch ? l1g.getInt('isEle', -4) : -2"),
)

process.rlxEGEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "NonIsolated"),
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","NonIsolated")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)

process.isoEGEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEleIso"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Isolated")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)



process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra *
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerUCT
    *process.createGenParticlesEle
    *process.rlxEGEfficiency
    *process.createGenParticlesEleIso
    *process.isoEGEfficiency

)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_*_*_ReRunningL1',
          'keep *_l1extraParticles*_*_*') 
)

process.out = cms.EndPath(process.output)


