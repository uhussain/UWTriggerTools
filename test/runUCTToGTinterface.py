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

#"/store/user/tapas/2012-08-01-CRAB_ZEESkim/skim_11_1_2Z3.root",

"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/04F0C0E3-72CA-E111-A802-003048F117EC.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/08361F06-71CA-E111-B50E-003048F1C424.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/0AB01403-71CA-E111-B705-003048F110BE.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/0CC9AA6C-74CA-E111-B56E-BCAEC518FF52.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/0E2082FE-70CA-E111-AA9C-00215AEDFD98.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/102F8715-71CA-E111-A9EF-003048D2BC5C.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/128B8EB8-72CA-E111-A703-BCAEC518FF40.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/12C53B00-71CA-E111-9741-001D09F253C0.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/1467473D-75CA-E111-A572-003048F118C4.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/18F82005-71CA-E111-AC70-001D09F34488.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/1A9A53BE-72CA-E111-B173-003048D2C0F4.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/1CCDBF70-74CA-E111-8334-0025901D62A0.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/22E5C4B8-72CA-E111-894E-003048F11DE2.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/2421D4B6-72CA-E111-B939-003048CFB40C.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/24CA1463-73CA-E111-A700-003048F11C28.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/26A74020-74CA-E111-9778-5404A63886B9.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/2AA64072-74CA-E111-BEE7-BCAEC518FF5F.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/2AB1ACB6-72CA-E111-8361-5404A638869B.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/3848706F-74CA-E111-AEAE-001D09F2527B.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/3C4A28B8-72CA-E111-827A-003048F118AA.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/3CA3E220-74CA-E111-A7A9-5404A6388699.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/44482AFD-70CA-E111-9F5D-003048F117B6.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/446FB0B7-72CA-E111-9432-003048F11C5C.root",
"/store/data/Run2012C/ZeroBias1/RAW/v1/000/198/609/4A9F3A1E-74CA-E111-A126-BCAEC518FF8F.root"
                             )   
                             )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# Load emulation and RECO sequences
#process.load("L1Trigger.UCT2015.emulationMC_cfi") 
process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
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



process.gctReEmulDigisNOTPUCORRECTED=cms.EDProducer("UCT2015GctCandsProducer",
            egRelaxed = cms.InputTag("UCT2015Producer","RelaxedEGUnpacked"),
    egIsolated  = cms.InputTag("UCT2015Producer","IsolatedEGUnpacked"),
            tauRelaxed = cms.InputTag("UCT2015Producer","RelaxedTauUnpacked"), # this collection is ignored in the final output, GT constraints 
    tauIsolated  = cms.InputTag("UCT2015Producer","IsolatedTauUnpacked"),
#    jetSource  = cms.InputTag("UCT2015Producer","CorrJetUnpacked"), # default are corrected jets 
    jetSource  = cms.InputTag("UCT2015Producer","JetUnpacked"),
    setSource  = cms.InputTag("UCT2015Producer","SETUnpacked"),
    metSource  = cms.InputTag("UCT2015Producer","METUnpacked"),
    shtSource  = cms.InputTag("UCT2015Producer","SHTUnpacked"),
    mhtSource  = cms.InputTag("UCT2015Producer","MHTUnpacked")
)



process.l1extraParticlesReEmul = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("gtDigis"),
    etTotalSource = cms.InputTag("gctReEmulDigis"),
    nonIsolatedEmSource = cms.InputTag("gctReEmulDigis","nonIsoEm"),
    etMissSource = cms.InputTag("gctReEmulDigis"),
    htMissSource = cms.InputTag("gctReEmulDigis"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("gctReEmulDigis","forJets"),
    centralJetSource = cms.InputTag("gctReEmulDigis","cenJets"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("gctReEmulDigis","tauJets"),
    isolatedEmSource = cms.InputTag("gctReEmulDigis","isoEm"),
    etHadSource = cms.InputTag("gctReEmulDigis"),
    hfRingEtSumsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    hfRingBitCountsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    centralBxOnly = cms.bool(True),
    ignoreHtMiss = cms.bool(False)
)


process.l1extraParticlesReEmulNOTPUCORRECTED = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("gtDigis"),
    etTotalSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"),
    nonIsolatedEmSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED","nonIsoEm"),
    etMissSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"),
    htMissSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED","forJets"),
    centralJetSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED","cenJets"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED","tauJets"),
    isolatedEmSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED","isoEm"),
    etHadSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"),
    hfRingEtSumsSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"), # these are empty
    hfRingBitCountsSource = cms.InputTag("gctReEmulDigisNOTPUCORRECTED"), # these are empty
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


process.l1ExtraTreeProducerReEmulNOTPUCORRECTED = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED"),
   metLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticlesReEmulNOTPUCORRECTED"),
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
    fileName = cms.string('L1TreeReRun_NOPATCH_WithoutIDAfter63.root')
)

from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis

process.gtReEmulDigis   = gtDigis.clone()

process.gtReEmulDigis.GmtInputTag  = cms.InputTag("gtDigis") # this is original GMT info from DATA (GMT is read out by GT FED)
process.gtReEmulDigis.GctInputTag  = cms.InputTag("gctReEmulDigis")
process.gtReEmulDigis.EmulateBxInEvent = cms.int32(1)

process.gtReEmulDigisNOTPUCORRECTED   = gtDigis.clone()

process.gtReEmulDigisNOTPUCORRECTED.GmtInputTag  = cms.InputTag("gtDigis") # this is original GMT info from DATA (GMT is read out by GT FED)
process.gtReEmulDigisNOTPUCORRECTED.GctInputTag  = cms.InputTag("gctReEmulDigisNOTPUCORRECTED")
process.gtReEmulDigisNOTPUCORRECTED.EmulateBxInEvent = cms.int32(1)


process.testMET = cms.EDAnalyzer("TestEtHt",
        originalL1Extra=cms.InputTag("l1extraParticles","MET"),
        uctCands=cms.InputTag("UCT2015Producer","METUnpacked"),
        uctCands2=cms.InputTag("UCT2015Producer","SETUnpacked"),
        modifiedL1Extra=cms.InputTag("l1extraParticlesReEmul","MET"),
)

process.testMHT = cms.EDAnalyzer("TestEtHt",
        originalL1Extra=cms.InputTag("l1extraParticles","MHT"),
        uctCands=cms.InputTag("UCT2015Producer","MHTUnpacked"),
        uctCands2=cms.InputTag("UCT2015Producer","SHTUnpacked"),
        modifiedL1Extra=cms.InputTag("l1extraParticlesReEmul","MHT"),
)




process.p1 = cms.Path(
    process.emulationSequence 
    *process.gtEvmDigis
    *process.dttfDigis
    *process.csctfDigis
    *process.gctReEmulDigis
    *process.gctReEmulDigisNOTPUCORRECTED
    *process.gtReEmulDigis
    *process.gtReEmulDigisNOTPUCORRECTED
    *process.scalersRawToDigi*
    process.l1extraParticles*
    process.l1extraParticlesReEmul*  
    process.l1extraParticlesReEmulNOTPUCORRECTED*
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerReEmul*
    process.l1ExtraTreeProducerReEmulNOTPUCORRECTED

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


