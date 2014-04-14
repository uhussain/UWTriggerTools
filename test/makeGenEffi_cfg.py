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
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/666531DD-9D75-E311-B95D-0025907DC9DC.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/D83DBA64-DE75-E311-88CD-0025907DC9CC.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/8EDF073A-8275-E311-80F7-00266CFFA25C.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/3415AE4D-B075-E311-865F-003048F0E194.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/B8AF4854-8275-E311-886D-00266CF327C4.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/9E9E1E8A-B075-E311-BFA2-00266CFFA048.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/88B7F378-9E75-E311-9F87-0025904B0FB6.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/606195D5-AF75-E311-8399-0025907DCA7E.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/C611F153-D375-E311-9436-003048F0E5B4.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/D28900FF-9E75-E311-9425-00266CF32F00.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/E29B3934-8275-E311-9202-0025901D42BC.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/C2BA874A-B075-E311-8964-0025901D4C92.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/A808CBDF-8F75-E311-8365-003048F0EBBE.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/7C3B05F1-BA75-E311-9495-003048D4610C.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/9E2AD326-B075-E311-B2AD-002481E0D6A6.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/06C609FB-B075-E311-8F3C-003048F0E59E.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/5A2B6818-DB75-E311-96CA-003048D4365C.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/866703E3-8F75-E311-A836-0025904B12B2.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/B6DD64ED-AF75-E311-8B96-00266CF32A00.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/F61F7047-8275-E311-A2D2-0025904B12F0.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/D285F57B-A475-E311-B032-00266CFFB868.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/F285AC91-DA75-E311-803F-003048C6903C.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/2A1BB24D-B075-E311-A5E9-003048F0E18C.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/144FD0F7-BA75-E311-8A18-00266CFFB390.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20000/CCFD4893-4075-E311-8A60-0025907FD2BA.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20000/BC8576EE-3575-E311-8AA7-00266CF33340.root",
"/store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20000/D8EE3250-9775-E311-9096-002481E0D9BC.root",
                             )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V2::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
#process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")
process.gctUCTDigis.jetSource =cms.InputTag("UCT2015Producer","CorrJetUnpacked")


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
        MinPtThreshold=cms.untracked.double(5),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(100)
)


process.createGenParticlesEleIso =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(5),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(1)
)

process.createGenParticlesTau =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(10),
        GenLevelSelect=cms.untracked.int32(15),
        GenLevelStatus=cms.untracked.int32(3),
        MaxIsolation=cms.untracked.double(100)
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


process.rlxTauEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesTau"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Tau"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Tau")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)


process.rlxTauL1TauPlusJetsEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesTau"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Central"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Tau")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)



process.compareEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticlesUCT","NonIsolated")
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
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
    recoFile  =cms.bool(False)
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
    recoFile  =cms.bool(False)
)

process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra *
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerUCT
    *process.createGenParticlesEle
    *process.rlxEGEfficiency
    *process.compareEfficiency
    *process.createGenParticlesEleIso
    *process.isoEGEfficiency
    *process.isoEGEfficiencyGENRLX
    * process.cleanGenJets
    *process.jetEfficiency            
    *process.uctSumsEfficiency
    *process.l1SumsEfficiency
    *process.createGenParticlesTau
    *process.rlxTauEfficiency
    *process.rlxTauL1TauPlusJetsEfficiency
    *process.rlxTauEfficiency  
)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_*_*_ReRunningL1',
          'keep *_l1extraParticles*_*_*') 
)

#process.out = cms.EndPath(process.output)


