import FWCore.ParameterSet.Config as cms

#########################################################################
##-------- Find Isolated Muons and Electrons 

tightMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt>10 && isGlobalMuon && isTrackerMuon && abs(eta)<2.5"
                     " && globalTrack().normalizedChi2<10"
                     " && globalTrack().hitPattern().numberOfValidMuonHits>0"
                     " && globalTrack().hitPattern().numberOfValidPixelHits>0"
                     " && numberOfMatchedStations>1"
                     " && globalTrack().hitPattern().trackerLayersWithMeasurement>5"	
		     " && (pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2"

	)
)

tightElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( "gsfElectrons" ),
    cut = cms.string(
    " (et>20.0)"
    " && (gsfTrack.trackerExpectedHitsInner.numberOfHits==0 && !(-0.02<convDist<0.02 && -0.02<convDcot<0.02))"
    " && ( dr03TkSumPt +  dr03EcalRecHitSumEt+  dr03HcalTowerSumEt)/pt  < 0.15 "+
    " && (  (isEB"
    " && abs(sigmaIetaIeta)<0.01 &&  abs(deltaPhiSuperClusterTrackAtVtx)<0.03 && abs(deltaEtaSuperClusterTrackAtVtx)<0.004 )"
    " || (isEE"
    " && abs(sigmaIetaIeta)<0.03 &&   abs(deltaPhiSuperClusterTrackAtVtx)<0.02 &&abs(deltaEtaSuperClusterTrackAtVtx)<0.005 )"
    ")"
    )
)


##########################################################################
##-------- Remove electrons and muons from jet collection ----------------------
ak5PFJetsNOMuons = cms.EDProducer("PFJetCleaner",
    srcJets = cms.InputTag("ak5PFJets"),
    module_label = cms.string(""),
    idLevel = cms.int32(0), # No Jet Id Required
    etaMin  =  cms.double(0.0),
    etaMax  =  cms.double(2.7),
    ptMin   =  cms.double(0.0),                                 
    srcObjects = cms.VInputTag(cms.InputTag("tightMuons")),
    deltaRMin = cms.double(0.3)
)

ak5PFJetsNOMuonsNOElectrons = cms.EDProducer("PFJetCleaner",
    srcJets = cms.InputTag("ak5PFJetsNOMuons"),
    module_label = cms.string(""),
    idLevel = cms.int32(0),  # No Jet Id required
    etaMin  =  cms.double(0.0),
    etaMax  =  cms.double(100),
    ptMin   =  cms.double(0.0),                                 
    srcObjects = cms.VInputTag(cms.InputTag("tightElectrons")),
    deltaRMin = cms.double(0.5)
)


cleanJets=cms.Sequence(tightMuons+tightElectrons+ak5PFJetsNOMuons+ak5PFJetsNOMuonsNOElectrons)



