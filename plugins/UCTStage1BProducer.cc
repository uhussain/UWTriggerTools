// -*- C++ -*-
//
// Package:    UCTStage1BProducer
// Class:      UCTStage1BProducer
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Thu Jun  7 13:29:52 CDT 2012
//


// system include files
#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "L1Trigger/UCT2015/interface/UCTCandidate.h"
#include "L1Trigger/UCT2015/interface/helpers.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

class UCTStage1BProducer : public edm::EDProducer {
public:

  // Concrete collection of L1Gobjects (with extra tuning information)
  typedef vector<UCTCandidate> UCTCandidateCollection;
  typedef std::auto_ptr<UCTCandidateCollection> UCTCandidateCollectionPtr;

  explicit UCTStage1BProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // Helper methods

  void makeEGs();
  void makeTaus();
  void puSubtraction();

  // ----------member data ---------------------------

  bool puCorrect;

  unsigned int egSeed;
  unsigned int tauSeed;

  double regionalHoECut;
  double egRelativeRgnIsolationCut;
  double egRelativeJetIsolationCut;
  double egRelativeEMRgnIsolationCut;
  double egRelativeEMJetIsolationCut;
  double tauRelativeRgnIsolationCut;
  double tauRelativeJetIsolationCut;
  double tauRelativeEMRgnIsolationCut;
  double tauRelativeEMJetIsolationCut;

  double egLSB_;
  double tauLSB_;
  double regionLSB_;

  double puLevelUIC;
  double puLevel;

  Handle<L1CaloRegionCollection> newRegions;
  Handle<L1CaloRegionCollection> newEMRegions;
  Handle<L1CaloEmCollection> tauCands;
  Handle<UCTCandidateCollection> emClusters;

  list<UCTCandidate> rlxEGList;
  list<UCTCandidate> isoEGList;

  list<UCTCandidate> rlxTauList;
  list<UCTCandidate> isoTauList;

};

UCTStage1BProducer::UCTStage1BProducer(const edm::ParameterSet& iConfig) :
  puCorrect(iConfig.getParameter<bool>("puCorrect")),
  egSeed(iConfig.getParameter<unsigned int>("egSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  regionalHoECut(iConfig.getParameter<double>("regionalHoECut")),
  egRelativeRgnIsolationCut(iConfig.getParameter<double>("egRelativeRgnIsolationCut")),
  egRelativeJetIsolationCut(iConfig.getParameter<double>("egRelativeJetIsolationCut")),
  egRelativeEMRgnIsolationCut(iConfig.getParameter<double>("egRelativeEMRgnIsolationCut")),
  egRelativeEMJetIsolationCut(iConfig.getParameter<double>("egRelativeEMJetIsolationCut")),
  tauRelativeRgnIsolationCut(iConfig.getParameter<double>("tauRelativeRgnIsolationCut")),
  tauRelativeJetIsolationCut(iConfig.getParameter<double>("tauRelativeJetIsolationCut")),
  tauRelativeEMRgnIsolationCut(iConfig.getParameter<double>("tauRelativeEMRgnIsolationCut")),
  tauRelativeEMJetIsolationCut(iConfig.getParameter<double>("tauRelativeEMJetIsolationCut")),
  egLSB_(iConfig.getParameter<double>("egLSB")),
  tauLSB_(iConfig.getParameter<double>("tauLSB")),
  regionLSB_(iConfig.getParameter<double>("regionLSB"))
{
  // Also declare we produce unpacked collections (which have more info)
  produces<UCTCandidateCollection>( "RelaxedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "RelaxedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedTauUnpacked" ) ;
}


// ------------ method called for each event  ------------
void
UCTStage1BProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iEvent.getByLabel("uctDigis", newRegions);
  iEvent.getByLabel("uctDigis", tauCands);
  iEvent.getByLabel("UCT2015EClusterProducer", "EClustersUnpacked", emClusters);
  iEvent.getByLabel("UCT2015EClusterProducer", "ERegions", newEMRegions);

  makeEGs();
  makeTaus();

  UCTCandidateCollectionPtr unpackedRlxTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedRlxEGs(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoEGs(new UCTCandidateCollection);

  for(list<UCTCandidate>::iterator rlxTau = rlxTauList.begin();
      rlxTau != rlxTauList.end();
      rlxTau++) {
    unpackedRlxTaus->push_back(*rlxTau);
  }
  for(list<UCTCandidate>::iterator isoTau = isoTauList.begin();
      isoTau != isoTauList.end();
      isoTau++) {
    unpackedIsoTaus->push_back(*isoTau);
  }
  for(list<UCTCandidate>::iterator rlxEG = rlxEGList.begin();
      rlxEG != rlxEGList.end();
      rlxEG++) {
    unpackedRlxEGs->push_back(*rlxEG);
  }
  for(list<UCTCandidate>::iterator isoEG = isoEGList.begin();
      isoEG != isoEGList.end();
      isoEG++) {
    unpackedIsoEGs->push_back(*isoEG);
  }

  iEvent.put(unpackedRlxTaus, "RelaxedTauUnpacked");
  iEvent.put(unpackedIsoTaus, "IsolatedTauUnpacked");
  iEvent.put(unpackedRlxEGs, "RelaxedEGUnpacked");
  iEvent.put(unpackedIsoEGs, "IsolatedEGUnpacked");

}

void UCTStage1BProducer::puSubtraction()
{
  puLevel = 0;
  puLevelUIC = 0;
  double r_puLevelUIC=0.0;

  int puCount = 0;
  double Rarea=0.0;

  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++){
    double regionEt = newRegion->et()*regionLSB_;
    double puETMax = 10;
    if(regionEt <= puETMax) {
      puLevel += regionEt;
      puCount++;
      r_puLevelUIC += regionEt;
      Rarea += getRegionArea(newRegion->gctEta());
    }
  }
  puLevel = puLevel / puCount;
  r_puLevelUIC = r_puLevelUIC / Rarea;
}

void UCTStage1BProducer::makeEGs() {
  rlxEGList.clear();
  isoEGList.clear();
  for(vector<UCTCandidate>::const_iterator emCluster = emClusters->begin();
      emCluster != emClusters->end();
      emCluster++) {
    if(emCluster->et() > egSeed) {
      double emClusterEt = emCluster->et();
      unsigned emClusterRegionIPhi = emCluster->getInt("rgnPhi");
      unsigned emClusterRegionIEta = emCluster->getInt("rgnEta");
      double C = 0;
      double N = 0;
      double S = 0;
      double E = 0;
      double W = 0;
      double NE = 0;
      double SE = 0;
      double NW = 0;
      double SW = 0;
      for(L1CaloRegionCollection::const_iterator region = newRegions->begin();
	  region != newRegions->end(); region++) {

        // in the physical scale.
        double regionEt = region->et()*regionLSB_;

	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare emCluster to the total E+H within the region to define rgnIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      double rgnIsolation = std::min(0., C - emClusterEt);

      double jetIsolation = C + N + S + E + W + NE + NW + SE + SW - emClusterEt;
      // Now compute EM only isolations
      C = 0;
      N = 0;
      S = 0;
      E = 0;
      W = 0;
      NE = 0;
      SE = 0;
      NW = 0;
      SW = 0;
      for(L1CaloRegionCollection::const_iterator region = newEMRegions->begin();
	  region != newEMRegions->end(); region++) {
        // FIXME - double check this is the correct scale.
        double regionEt = region->et()*regionLSB_;
	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare emCluster to the total E within the region to define rgnEMIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      double rgnEMIsolation = std::min(0., C - emClusterEt);
      double jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - emCluster->pt();
      // Every emCluster which passes HoECut is an egObject
      double h = C - emCluster->pt();
      double e = emClusterEt;
      if((h / e) < regionalHoECut) {
        // make a copy of the candidate
        UCTCandidate egCand = *emCluster;
	// egObject that passes all isolations is an isolated egObject
	double relativeRgnEMIsolation = (double) rgnEMIsolation / (double) emCluster->pt();
	double relativeJetEMIsolation = (double) jetEMIsolation / (double) emCluster->pt();
	double relativeRgnIsolation = ((double) rgnIsolation) / ((double) emCluster->pt());
	double relativeJetIsolation = ((double) jetIsolation) / ((double) emCluster->pt());
        egCand.setFloat("rgnEMIsolation", rgnEMIsolation);
        egCand.setFloat("jetEMIsolation", jetEMIsolation);
        egCand.setFloat("rgnIsolation", rgnIsolation);
        egCand.setFloat("jetIsolation", jetIsolation);
        egCand.setFloat("h", h);
        egCand.setFloat("e", e);
	rlxEGList.push_back(egCand);

	if(relativeRgnIsolation < egRelativeRgnIsolationCut &&
	   relativeJetIsolation < egRelativeJetIsolationCut &&
	   relativeRgnEMIsolation < egRelativeEMRgnIsolationCut &&
	   relativeJetEMIsolation < egRelativeEMJetIsolationCut) {
	  isoEGList.push_back(rlxEGList.back());
	}
      }
    }
  }
  rlxEGList.sort();
  isoEGList.sort();
  rlxEGList.reverse();
  isoEGList.reverse();
}

void UCTStage1BProducer::makeTaus() {
  rlxTauList.clear();
  isoTauList.clear();
  for(L1CaloEmCollection::const_iterator tauCand = tauCands->begin();
      tauCand != tauCands->end(); tauCand++){
    double tauCandEt = tauCand->rank() * tauLSB_;
    if(tauCandEt > tauSeed) {
      unsigned tauCandRegionIPhi = tauCand->regionId().iphi();
      unsigned tauCandRegionIEta = tauCand->regionId().ieta();
      double C = 0;
      double N = 0;
      double S = 0;
      double E = 0;
      double W = 0;
      double NE = 0;
      double SE = 0;
      double NW = 0;
      double SW = 0;
      for(L1CaloRegionCollection::const_iterator region = newRegions->begin();
	  region != newRegions->end(); region++) {
        double regionEt = region->et()*regionLSB_;
	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare tauCand to the total E+H within the region to define rgnIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      double rgnIsolation = std::min(0., C - tauCandEt);
      double jetIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEt;
      // Now compute EM only isolations
      C = 0;
      N = 0;
      S = 0;
      E = 0;
      W = 0;
      NE = 0;
      SE = 0;
      NW = 0;
      SW = 0;
      // FIXME - make sure the LSB here is correct.
      for(L1CaloRegionCollection::const_iterator region = newEMRegions->begin();
	  region != newEMRegions->end(); region++) {
        double regionEt = region->et()*regionLSB_;
	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare tauCand to the total E within the region to define rgnEMIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      double rgnEMIsolation = std::min(0., C - tauCandEt);
      unsigned jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEt;
      // Use finer grain position resolution if possible using the emClusters
      // with default value being center of the region
      // phi: 0-17 becomes 0-71; eta: 0-13 becomes 0-55;
      double tauCandEta = convertRegionEta(tauCandRegionIEta);
      double tauCandPhi = convertRegionPhi(tauCandRegionIPhi);

      double matchedEmClusterEt = -1;

      for(vector<UCTCandidate>::const_iterator emCluster = emClusters->begin();
	  emCluster != emClusters->end();
	  emCluster++) {
	if(emCluster->getInt("rgnPhi") == (int)tauCandRegionIPhi && emCluster->getInt("rgnEta") == (int)tauCandRegionIEta) {
	  tauCandEta = emCluster->eta();
	  tauCandPhi = emCluster->phi();
          matchedEmClusterEt = emCluster->pt();
          break;
	}
      }
      UCTCandidate theTau(tauCandEt, tauCandEta, tauCandPhi);
      // tauObject that passes all isolations is an isolated tauObject
      double relativeRgnEMIsolation = ((double) rgnEMIsolation) / ((double) tauCandEt);
      double relativeJetEMIsolation = ((double) jetEMIsolation) / ((double) tauCandEt);
      double relativeRgnIsolation = ((double) rgnIsolation) / ((double) tauCandEt);
      double relativeJetIsolation = ((double) jetIsolation) / ((double) tauCandEt);

      theTau.setFloat("rgnEMIsolation", rgnEMIsolation);
      theTau.setFloat("jetEMIsolation", jetEMIsolation);
      theTau.setFloat("rgnIsolation", rgnIsolation);
      theTau.setFloat("jetIsolation", jetIsolation);
      theTau.setFloat("emClusterEt", matchedEmClusterEt);
      theTau.setInt("rgnEta", tauCandRegionIEta);
      theTau.setInt("rgnPhi", tauCandRegionIPhi);

      rlxTauList.push_back(theTau);

      if(relativeRgnIsolation < tauRelativeRgnIsolationCut &&
	 relativeJetIsolation < tauRelativeJetIsolationCut &&
	 relativeRgnEMIsolation < tauRelativeEMRgnIsolationCut &&
	 relativeJetEMIsolation < tauRelativeEMJetIsolationCut) {
	isoTauList.push_back(rlxTauList.back());
      }
    }
  }
  rlxTauList.sort();
  isoTauList.sort();
  rlxTauList.reverse();
  isoTauList.reverse();
}

//define this as a plug-in
DEFINE_FWK_MODULE(UCTStage1BProducer);
