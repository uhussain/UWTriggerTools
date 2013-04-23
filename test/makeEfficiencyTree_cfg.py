#flake8: noqa
'''

Generate trees for measuring and comparing L1 and UCT efficiencies with
respect to RECO objects.

Usage:

    ./makeEfficiencyTree_cfg.py

Optional arguments:

    inputFiles=myFile.root outputFile=outputFile.root maxEvents=-1

Authors: L. Dodd, N. Woods, I. Ojalvo, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

# Get command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# Set useful defaults
#options.inputFiles = '/store/user/tapas/ETauSkim/skim_12_1_erV.root'
options.outputFile = "uct_efficiency_tree.root"
options.register(
    'eicIsolationThreshold',
    3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "EIC Isolation threshold")
options.register(
    "stage1B",
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "If 1, enable production of Stage1B trees"
)
options.register(
    "stage1",
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "If 0, *disable* production of Stage1 trees"
)
options.register(
    'ecalCalib',
    'CALIB_V4',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Can be CALIB_V1, CALIB_V3, or CALIB_V4')
options.register(
    'eicCardHcalOnly',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'If 1, turn off the ECAL for the stage1 EGTau path.')
options.register(
    'isMC',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.')

options.parseArguments()

process = cms.Process("L1UCTEfficiency")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# Load the correct global tag, based on the release
if 'CMSSW_6' in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = 'POSTLS161_V12::All'
    print "Using global tag for upgrade MC: %s" % process.GlobalTag.globaltag
    if not options.isMC:
        raise ValueError("There is no data in CMSSW 6, you must mean isMC=1")
else:
    if not options.isMC:
        # CMSSW 5 data
        process.GlobalTag.globaltag = 'GR_R_53_V21::All'
    else:
        # CMSSW 5 MC
        process.GlobalTag.globaltag = 'START53_V7B::All'
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    print "Using global tag for 52X data: %s" % process.GlobalTag.globaltag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

# Load emulation and RECO sequences
if not options.isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

process.load("L1Trigger.UCT2015.recoObjects_cfi")

# Determine which calibration to use
from L1Trigger.UCT2015.emulation_cfi import \
        eg_calib_v1, eg_calib_v3, eg_calib_v4

calib_map = {
    'CALIB_V1': eg_calib_v1,
    'CALIB_V3': eg_calib_v3,
    'CALIB_V4': eg_calib_v4
}

ecal_calibration = calib_map[options.ecalCalib]
process.RCTConfigProducers.eGammaECalScaleFactors = ecal_calibration
process.RCTConfigProducers.jetMETECalScaleFactors = ecal_calibration
process.UCT2015EClusterProducer.ecalCalibration = ecal_calibration

if options.eicCardHcalOnly:
    print "Disabling ECAL in Stage1 EGTau path"
    # Set all scale factors to 0.
    process.RCTConfigProducers.eGammaECalScaleFactors = [
        0.0 for _ in eg_calib_v1]

# Common branches to add to the ntuple
common_ntuple_branches = cms.PSet(
    index = cms.string("index"), # Index of reco object in the event
    nRecoObjects = cms.string("nTotalObjects"), # Number of reco objects in the event
    nPVs = cms.string("nPVs"), # number of reco'ed vertices in the event

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

    # For tuning isolation and PU subtraction
    l1gPU = cms.string("? l1gMatch ? l1g.getFloat('puLevel', -4) : -2"),
    l1gPUUIC = cms.string("? l1gMatch ? l1g.getFloat('puLevelUIC', -4) : -2"),
    l1gRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedRegionEt', -4) : -2"),

    l1gEtaCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnEta') : 0", "I"),
    l1gPhiCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnPhi') : 0", "I"),

    l1gDPhi = cms.string("? l1gMatch ? deltaPhi(l1g.phi, reco.phi) : -1"),
    l1gDEta = cms.string("? l1gMatch ? l1g.eta - reco.eta : -10"),
    l1gDR = cms.string("? l1gMatch ? deltaR(l1g.eta, l1g.phi, reco.eta, reco.phi) : -1"),
)

# Specific to EG tau objects
egtau_branches = cms.PSet(
    l1gJetPt = cms.string("? l1gMatch ? l1g.getFloat('associatedJetPt', -4) : -2"),
    l1gEllIso = cms.string("? l1gMatch ? l1g.getInt('ellIsolation', -4) : -2"),
    l1gTauVeto = cms.string("? l1gMatch ? l1g.getInt('tauVeto', -4) : -2"),
    l1gMIP = cms.string("? l1gMatch ? l1g.getInt('mipBit', -4) : -2"),
)

stage1b_branches = cms.PSet(
    l1gPUEM = cms.string("? l1gMatch ? l1g.getFloat('puLevelEM', -4) : -2"),
    l1gPUUICEM = cms.string("? l1gMatch ? l1g.getFloat('puLevelUICEM', -4) : -2"),
    l1gEffArea = cms.string("? l1gMatch ? l1g.getFloat('effArea', -4) : -2"),
    l1gRegionEtEM = cms.string("? l1gMatch ? l1g.getFloat('associatedRegionEtEM', -4) : -2"),
    l1gJetPtEM = cms.string("? l1gMatch ? l1g.getFloat('associatedJetPtEM', -4) : -2"),
    # only defined for taus, EG objects are EM clusters by defintion.
    l1gEmClusterEt = cms.string("? l1gMatch ? l1g.getFloat('emClusterEt', -4) : -2"),
    l1gEmClusterCenterEt = cms.string("? l1gMatch ? l1g.getFloat('emClusterCenterEt', -4) : -2"),
    l1gEmClusterStripEt = cms.string("? l1gMatch ? l1g.getFloat('emClusterStripEt', -4) : -2"),
    l1gEmClusterFG = cms.string("? l1gMatch ? l1g.getInt('emClusterCenterFG', -4) : -2"),

    l1gHighestCenter2x1 = cms.string("? l1gMatch ? l1g.getFloat('highestCenter2x1Et', -4) : -2"),
    l1gHighestNeighbor2x1 = cms.string("? l1gMatch ? l1g.getFloat('highestNeighbor2x1Et', -4) : -2"),

    l1g4RegionEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).totalEt : -2"),
    l1g4RegionEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).totalEtEcal : -2"),
    l1g4RegionPattern = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).patternPass : -2"),
    l1g4RegionMips = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).numberOfMips : -2"),
    l1g4RegionLowestEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).lowestRegionEt : -2"),
    l1g4RegionLowestEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).lowestRegionEtEcal : -2"),
    l1g4RegionNumRegions = cms.string("? l1gMatch ? l1g.regionDiscriminant(4).numberOfRegions : -2"),

    l1g3RegionEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).totalEt : -2"),
    l1g3RegionEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).totalEtEcal : -2"),
    l1g3RegionPattern = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).patternPass : -2"),
    l1g3RegionMips = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).numberOfMips : -2"),
    l1g3RegionLowestEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).lowestRegionEt : -2"),
    l1g3RegionLowestEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).lowestRegionEtEcal : -2"),
    l1g3RegionNumRegions = cms.string("? l1gMatch ? l1g.regionDiscriminant(3).numberOfRegions : -2"),

    l1g2RegionEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).totalEt : -2"),
    l1g2RegionEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).totalEtEcal : -2"),
    l1g2RegionPattern = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).patternPass : -2"),
    l1g2RegionMips = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).numberOfMips : -2"),
    l1g2RegionLowestEt = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).lowestRegionEt : -2"),
    l1g2RegionLowestEtEcal = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).lowestRegionEtEcal : -2"),
    l1g2RegionNumRegions = cms.string("? l1gMatch ? l1g.regionDiscriminant(2).numberOfRegions : -2"),
)

# Keep track of electron isolation values
electron_branches = cms.PSet(
    dr03TkSumPt  = cms.string("reco.dr03TkSumPt"),
    dr03EcalRecHitSumEt  = cms.string("reco.dr03EcalRecHitSumEt"),
    dr03HcalTowerSumEt  = cms.string("reco.dr03HcalTowerSumEt"),
    dr03CombinedEt  = cms.string("reco.dr03TkSumPt + reco.dr03EcalRecHitSumEt + reco.dr03HcalTowerSumEt"),
)

# Keep track of information about the ECAL/HCAL composition of taus
tau_branches = cms.PSet(
    emFraction = cms.string("reco.emFraction"),
    decayMode = cms.string("reco.decayMode"),
    # EK - as far as I can tell, this does not use the lead track at all
    hcal = cms.string("reco.hcalTotOverPLead"),
)

# Define the tree producers
process.isoTauEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)
# Define the tree producers
process.rlxTauEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)

# Define the tree producers
process.rlxTauPlusJetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Central"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)

# Note that the input electron collection is not isolated, this needs to be done
# at the ntuple level.
process.isoEGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Isolated")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

process.rlxEGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "NonIsolated"),
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

# So we can compare relaxed UCT + ntuple isolation cuts versus stock L1 IsoEG
process.rlxUCTisoL1EGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

# Package all of the lepton efficiencies into one sequence
process.leptonEfficiencies = cms.Sequence(
    #process.isoTauEfficiency *
    process.rlxTauEfficiency *
    #process.rlxTauPlusJetEfficiency *
    #process.isoEGEfficiency *
    process.rlxEGEfficiency *
    process.rlxUCTisoL1EGEfficiency
)

process.jetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "JetUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)

process.corrjetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "CorrJetUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)

process.highPtPF = cms.EDFilter(
    "GenericPFCandidateSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("charge != 0 && pt > 10")
)

process.printPF = cms.EDAnalyzer(
    "CandInfoPrinter",
    src = cms.InputTag("highPtPF"),
    printHeader=cms.bool(True),
    pt = cms.string("pt"),
    eta = cms.string("eta"),
    pdgId = cms.string("pdgId")
)


process.printTaus = cms.EDAnalyzer(
    "CandInfoPrinter",
    src = cms.InputTag("hpsPFTauProducer"),
    printHeader=cms.bool(True),
    pt = cms.string("pt"),
    ncands = cms.string("signalPFCands.size"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    dm = cms.string("decayMode"),
)

process.printTPGs = cms.EDFilter(
    "TPGDebugger",
    ecalSrc = cms.InputTag("ecalDigis:EcalTriggerPrimitives"),
    #hcalSrc = cms.InputTag("hackHCALMIPs"),
    hcalSrc = cms.InputTag("simHcalTriggerPrimitiveDigis"),
    toPrint = cms.VPSet(
        # everything in all events
        cms.PSet(
            run = cms.uint32(0),
            minIEta = cms.int32(-1000),
            maxIEta = cms.int32(1000),
            maxIPhi = cms.int32(1000),
            minIPhi = cms.int32(-1000),
        )
    )
)

#process.hackHCALMIPs.src = "simHcalTriggerPrimitiveDigis"
process.pionEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("highPtPF"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Tau"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCTStage1BProducer", "RelaxedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        stage1b_branches
    )
)


process.p1 = cms.Path(
    process.recoObjects *
    process.emulationSequence *
    #process.printTaus *
    #process.highPtPF *
    #process.printPF *
    #process.printTPGs *
    #process.dump *
    #process.pionEfficiency *
    process.jetEfficiency *
    process.corrjetEfficiency
)

if options.stage1:
    print "Building Stage1 trees"
    process.p1 += process.leptonEfficiencies

if options.stage1B:
    print "Building Stage1B trees"
    # Make a copy of the lepton efficiency trees using stage 1B inputs.
    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
    process.leptonEfficienciesStage1B = cloneProcessingSnippet(
        process, process.leptonEfficiencies, 'Stage1B')
    # Update input tags to the stage 1B producer
    for stage1BTreeMaker in [#process.isoTauEfficiencyStage1B,
                             process.rlxTauEfficiencyStage1B,
                             #process.isoEGEfficiencyStage1B,
                             process.rlxEGEfficiencyStage1B,
                             process.rlxUCTisoL1EGEfficiencyStage1B]:
        stage1BTreeMaker.l1GSrc[0].setModuleLabel("UCTStage1BProducer")
        for branch in stage1b_branches.parameterNames_():
            setattr(stage1BTreeMaker.ntuple, branch,
                    getattr(stage1b_branches, branch))
    # add the computation of stage1b trees
    process.p1 += process.leptonEfficienciesStage1B



################################################################################
###  Semileptonic ttbar skim for sums ###########################################
################################################################################

process.oneMuon = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tightMuons"),
    minNumber = cms.uint32(1),
)

process.jetsPt30 = cms.EDFilter(
    "PFJetSelector",
    src = cms.InputTag("ak5PFJetsNOMuons"),
    filter = cms.bool(True),
    cut = cms.string("pt > 30")
)

process.atLeastThreeJets = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("jetsPt30"),
    minNumber = cms.uint32(3),
)

# Computing CaloMET significance
process.load("RecoMET.METProducers.CaloMETSignif_cfi")
process.metsignificance.noHF = True

# Computing RECO level Sum ET and Sum HT
process.load("L1Trigger.UCT2015.PFSumET_cfi")

process.l1SumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    l1MHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticles", "MET"),
    # Evan said change l1METSigSrc to match recoMETSigSrc
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    #l1METSigSrc = cms.InputTag("metsignificance"),
    # fixme
    l1SHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticles", "MET"),
    recoMHTSrc = cms.InputTag("htMetAK5"),
    recoMETSrc = cms.InputTag("metNoHF"), # calomet
    recoMETSigSrc = cms.InputTag("metsignificance"), # calomet
    recoSHTSrc = cms.InputTag("pfSumET", "sht"),
    recoSETSrc = cms.InputTag("pfSumET", "set"),
)
process.uctSumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    l1MHTSrc = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    l1METSrc = cms.InputTag("UCT2015Producer", "METUnpacked"),
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("UCT2015Producer", "SHTUnpacked"),
    l1SETSrc = cms.InputTag("UCT2015Producer", "SETUnpacked"),

    recoMHTSrc = cms.InputTag("htMetAK5"),
    recoMETSrc = cms.InputTag("metNoHF"), # calomet
    recoMETSigSrc = cms.InputTag("metsignificance"), # calomet
    recoSHTSrc = cms.InputTag("pfSumET", "sht"),
    recoSETSrc = cms.InputTag("pfSumET", "set"),
)

# Make a version of UCT without PU corrections.
process.UCT2015ProducerNoPU = process.UCT2015Producer.clone(
    puCorrect = False
)
process.uctSumsNoPUEfficiency = process.uctSumsEfficiency.clone(
    l1MHTSrc = cms.InputTag("UCT2015ProducerNoPU", "MHTUnpacked"),
    l1METSrc = cms.InputTag("UCT2015ProducerNoPU", "METUnpacked"),
    l1METSigSrc = cms.InputTag("UCT2015ProducerNoPU", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("UCT2015ProducerNoPU", "SHTUnpacked"),
    l1SETSrc = cms.InputTag("UCT2015ProducerNoPU", "SETUnpacked"),
)

process.semileptonicTTBarPath = cms.Path(
    process.cleanJets *
    process.oneMuon *
    process.jetsPt30 *
    process.atLeastThreeJets *
    process.pfSumET *
    process.metsignificance *
    process.l1SumsEfficiency *
    process.uctSumsEfficiency
    # w/o PU corrections
    #process.UCT2015ProducerNoPU *
    #process.uctSumsNoPUEfficiency
)

process.schedule = cms.Schedule(
    process.p1,
    process.semileptonicTTBarPath
)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

eic = options.eicIsolationThreshold
print "Setting EIC threshold to %i" % eic
process.RCTConfigProducers.eicIsolationThreshold = eic
