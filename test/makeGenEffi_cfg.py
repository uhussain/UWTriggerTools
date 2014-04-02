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
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/144FD0F7-BA75-E311-8A18-00266CFFB390.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/2A1BB24D-B075-E311-A5E9-003048F0E18C.root"
                             )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V2::All'

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
    fileName = cms.string('genEffis2.root')
)

process.createGenParticlesEle =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(10),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(100)
)


process.createGenParticlesEleIso =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(10),
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

process.isoEGEfficiencyGENRLX = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
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

process.jetEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("cleanGenJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "Central"), cms.InputTag("l1extraParticlesUCT","Forward")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)

process.cleanGenJets = cms.EDProducer("GenJetClean")

process.l1SumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(False),
    l1MHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticles", "MET"),
    # Evan said change l1METSigSrc to match recoMETSigSrc
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    #l1METSigSrc = cms.InputTag("metsignificance"),
    # fixme
    l1SHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticles", "MET"),
    recoMHTSrc = cms.InputTag("genMetCalo"),
    recoMETSrc = cms.InputTag("genMetCalo"), # calomet
    recoMETSigSrc  = cms.InputTag("genMetCalo"),
    recoSHTSrc = cms.InputTag("genMetCalo"),
    recoSETSrc = cms.InputTag("genMetCalo"),
    recoPFMETSrc = cms.InputTag("genMetTrue"), # pfmet
)
process.uctSumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(False),
    l1MHTSrc = cms.InputTag("l1extraParticlesUCT", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticlesUCT", "MET"),
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("l1extraParticlesUCT", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticlesUCT", "MET"),
    recoMHTSrc = cms.InputTag("genMetCalo"),
    recoMETSrc = cms.InputTag("genMetCalo"), # calomet
    recoMETSigSrc  = cms.InputTag("genMetCalo"),
    recoSHTSrc = cms.InputTag("genMetCalo"),
    recoSETSrc = cms.InputTag("genMetCalo"),
    recoPFMETSrc = cms.InputTag("genMetTrue"), # pfmet
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
    *process.isoEGEfficiencyGENRLX
    * process.cleanGenJets
    *process.jetEfficiency            
    #*process.uctSumsEfficiency
    #*process.l1SumsEfficiency
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

#process.out = cms.EndPath(process.output)


