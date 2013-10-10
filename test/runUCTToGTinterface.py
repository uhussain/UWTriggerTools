'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("ReRunningL1")

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring("file:/hdfs/store/mc/Summer12/GluGluToHToTauTau_M-125_14TeV-powheg-pythia6/GEN-SIM-RAW-RECO/PU50_POSTLS161_V12-v1/10000/FE408615-3156-E211-8296-0026189438B9.root")
                             )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
#process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.gctReEmulDigis =cms.EDProducer("UCT2015GctCandsProducer",
            egRelaxed = cms.InputTag("UCT2015Producer","RelaxedEGUnpacked"),
    egIsolated  = cms.InputTag("UCT2015Producer","IsolatedEGUnpacked"),
            tauRelaxed = cms.InputTag("UCT2015Producer","RelaxedTauUnpacked"), # this collection is ignored in the final output, GT constraints 
    tauIsolated  = cms.InputTag("UCT2015Producer","IsolatedTauUnpacked"),
    jetSource  = cms.InputTag("UCT2015Producer","CorrJetUnpacked"), # default are corrected jets 
#    jetSource  = cms.InputTag("UCT2015Producer","JetUnpacked"),
    setSource  = cms.InputTag("UCT2015Producer","SETUnpacked"),
    metSource  = cms.InputTag("UCT2015Producer","METUnpacked"),
    shtSource  = cms.InputTag("UCT2015Producer","SHTUnpacked"),
    mhtSource  = cms.InputTag("UCT2015Producer","MHTUnpacked")
)

process.l1extraParticlesReEmul = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("gtDigis"),
    etTotalSource = cms.InputTag("gctReEmulDigis"),
    nonIsolatedEmSource = cms.InputTag("gctReEmulDigis","rlxEm"),
    etMissSource = cms.InputTag("gctReEmulDigis"),
    htMissSource = cms.InputTag("gctReEmulDigis"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("gctReEmulDigis","forJets"),
    centralJetSource = cms.InputTag("gctReEmulDigis","cenJets"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("gctReEmulDigis","isoTau"),
    isolatedEmSource = cms.InputTag("gctReEmulDigis","isoEm"),
    etHadSource = cms.InputTag("gctReEmulDigis"),
    hfRingEtSumsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    hfRingBitCountsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    centralBxOnly = cms.bool(True),
    ignoreHtMiss = cms.bool(False)
)

process.l1ExtraTreeProducerReEmul = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmul:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmul:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticlesReEmul:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticlesReEmul:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticlesReEmul:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticlesReEmul"),
   metLabel = cms.untracked.InputTag("l1extraParticlesReEmul:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticlesReEmul:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticlesReEmul"),
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
    fileName = cms.string('L1TreeReRun.root')
)

from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis

process.gtReEmulDigis   = gtDigis.clone()

process.gtReEmulDigis.GmtInputTag  = cms.InputTag("gtDigis") # this is original GMT info from DATA (GMT is read out by GT FED)
process.gtReEmulDigis.GctInputTag  = cms.InputTag("gctReEmulDigis")
process.gtReEmulDigis.EmulateBxInEvent = cms.int32(1)


process.p1 = cms.Path(
    process.emulationSequence 
    *process.gtEvmDigis
    *process.dttfDigis
    *process.csctfDigis
    *process.gctReEmulDigis
    *process.gtReEmulDigis
    *process.scalersRawToDigi*
    process.l1extraParticlesReEmul*    
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerReEmul
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


