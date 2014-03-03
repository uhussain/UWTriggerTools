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
#"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/1C84AE42-AD86-E311-8D05-00304867C1BA.root",
#"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/240A6A68-BB86-E311-B7C1-0025905A612E.root",
#"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/40070548-AF86-E311-B743-003048679164.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/BAAE7B3B-A186-E311-A40C-00304867926C.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/CE457DFF-B786-E311-8405-0025905A48D6.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/CE9BB547-AF86-E311-8E14-0025905A4964.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/9EAC35F8-A286-E311-B7FA-002618943948.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/AAD09E92-A186-E311-A4C5-002590596484.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/96C03CD0-B286-E311-9564-0025905A608C.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/2CBB06DF-B286-E311-87FD-0025905A60A6.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/D6034621-BC86-E311-B42A-0025905938D4.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/CC87AA2D-A486-E311-B043-001A92971B7C.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/D6F1788B-B886-E311-9016-003048678B34.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/7298FC54-A186-E311-9C1C-003048679006.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/2EED3630-A886-E311-AB10-0025905A611C.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/9A9B71B0-B486-E311-B89A-0025905964CC.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/D2582D20-A886-E311-B111-00261894388D.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/5A2A338A-B186-E311-9F99-00261894389A.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/6CC1BE27-BC86-E311-AB4B-0025905A48D8.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx50_POSTLS162_V1-v1/00003/64571E2D-B086-E311-9F9E-003048D15E2C.root"
                             )   
                             )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
#process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
process.load("Configuration.Geometry.GeometryIdeal_cff")


process.gctReEmulDigisOldTaus =cms.EDProducer("UCT2015GctCandsProducer",
            egRelaxed = cms.InputTag("UCT2015Producer","RelaxedEGUnpacked"),
    egIsolated  = cms.InputTag("UCT2015Producer","IsolatedEGUnpacked"),
            tauRelaxed = cms.InputTag("UCT2015Producer","RelaxedEcalSeedTauUnpacked"), # this collection is ignored in the final output, GT constraints 
    tauIsolated  = cms.InputTag("UCT2015Producer","IsolatedTauEcalSeedUnpacked"),
    jetSource  = cms.InputTag("UCT2015Producer","CorrJetUnpacked"), # default are corrected jets 
#    jetSource  = cms.InputTag("UCT2015Producer","JetUnpacked"),
    setSource  = cms.InputTag("UCT2015Producer","SETUnpacked"),
    metSource  = cms.InputTag("UCT2015Producer","METUnpacked"),
    shtSource  = cms.InputTag("UCT2015Producer","SHTUnpacked"),
    mhtSource  = cms.InputTag("UCT2015Producer","MHTUnpacked")
)




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

process.l1extraParticlesReEmulOldTaus = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("gtDigis"),
    etTotalSource = cms.InputTag("gctReEmulDigis"),
    nonIsolatedEmSource = cms.InputTag("gctReEmulDigis","nonIsoEm"),
    etMissSource = cms.InputTag("gctReEmulDigis"),
    htMissSource = cms.InputTag("gctReEmulDigis"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("gctReEmulDigis","forJets"),
    centralJetSource = cms.InputTag("gctReEmulDigis","cenJets"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("gctReEmulDigisOldTaus","tauJets"),
    isolatedEmSource = cms.InputTag("gctReEmulDigis","isoEm"),
    etHadSource = cms.InputTag("gctReEmulDigis"),
    hfRingEtSumsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    hfRingBitCountsSource = cms.InputTag("gctReEmulDigis"), # these are empty
    centralBxOnly = cms.bool(True),
    ignoreHtMiss = cms.bool(False)
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

process.l1ExtraTreeProducerReEmulOldTaus = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus"),
   metLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticlesReEmulOldTaus"),
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
    fileName = cms.string('L1TreeReRun_recoveringElectrons_brokenTaus.root')
)

from L1Trigger.GlobalTrigger.gtDigis_cfi import gtDigis

process.gtReEmulDigis   = gtDigis.clone()

process.gtReEmulDigis.GmtInputTag  = cms.InputTag("gtDigis") # this is original GMT info from DATA (GMT is read out by GT FED)
process.gtReEmulDigis.GctInputTag  = cms.InputTag("gctReEmulDigis")
process.gtReEmulDigis.EmulateBxInEvent = cms.int32(1)

process.gtReEmulDigisOldTaus   = gtDigis.clone()

process.gtReEmulDigisOldTaus.GmtInputTag  = cms.InputTag("gtDigis") # this is original GMT info from DATA (GMT is read out by GT FED)
process.gtReEmulDigisOldTaus.GctInputTag  = cms.InputTag("gctReEmulDigisOldTaus")
process.gtReEmulDigisOldTaus.EmulateBxInEvent = cms.int32(1)


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
    *process.gctReEmulDigisOldTaus
    *process.gtReEmulDigis
    *process.gtReEmulDigisOldTaus
    *process.scalersRawToDigi*
    process.l1extraParticles*
    process.l1extraParticlesReEmul*  
    process.l1extraParticlesReEmulOldTaus*
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerReEmul*
    process.l1ExtraTreeProducerReEmulOldTaus

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


