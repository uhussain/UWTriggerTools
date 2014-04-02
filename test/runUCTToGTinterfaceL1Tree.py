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
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/6AF2C1E2-DF7F-E311-B452-003048679162.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/9643B5C2-0280-E311-A66E-00261894394D.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/A0B9E0A6-C87F-E311-92E4-0026189438EF.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/4AEC3425-0580-E311-B5A2-0025905A611E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20002/9673477E-D57A-E311-ACD1-002618943983.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/341F0B18-157B-E311-A26C-0025905A6088.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/EAC88847-C27F-E311-8824-0026189438F5.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/0A53FC3E-0780-E311-9D84-0025905A60AA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/B08574C9-E57F-E311-99A5-003048D15DDA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/4A4B2971-C07F-E311-843E-0025905A6094.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/1EFB29C6-BE7F-E311-83C7-003048FF9AA6.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/6AD89D7E-A57F-E311-B867-0025905A611E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/AA99D298-B47F-E311-AB06-002618943901.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/9EF8E9D8-AA7F-E311-8589-0025905A60CA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/BAB14DA6-CF7C-E311-BDBA-0026189438A7.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20002/5068EC78-B47A-E311-8CEB-003048678FAE.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20001/306835C8-A679-E311-9A07-00261894387B.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/0CB7AE7E-A57F-E311-A1FF-00261894382A.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/BA95FC23-A87F-E311-BA73-001A928116D2.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/683C7D19-C97F-E311-AF6E-003048FF86CA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/049861F6-AB7F-E311-AF20-0025905A6134.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/CC19FEF5-BB7F-E311-B3BC-0025905A6064.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20001/B885C2C2-9779-E311-B93F-001A92971BC8.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/00995EB2-F57F-E311-8C2A-002618943945.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/1C99ABBB-EC7F-E311-8DE9-003048678B18.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/301A9AFB-0080-E311-9586-003048FF86CA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/F06FC13B-B17F-E311-98C8-0025905A60B0.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/ACAF37A2-BA7F-E311-8FC3-002618943927.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/94C3E755-B07F-E311-BA9D-003048FFD76E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/DCE31B33-AF7F-E311-BD2C-0025905A60A6.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/805BEA08-AE7F-E311-8553-0025905A607E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/78834016-B27F-E311-81D2-002590596490.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/10817AFD-A57F-E311-9910-0025905A6060.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/7A763396-AC7C-E311-BA84-002618943949.root",
                             )   
                             )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20000)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V2::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
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
    fileName = cms.string('L1TreeReRun_recoveringElectrons_long.root')
)

process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra *
    process.l1ExtraTreeProducer*
    process.l1ExtraTreeProducerUCT
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
