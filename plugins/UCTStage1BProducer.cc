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

#include "L1Trigger/UCT2015/src/L1GObject.h"
#include "L1Trigger/UCT2015/interface/helpers.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

class UCTStage1BProducer : public edm::EDProducer {
public:

  // Concrete collection of L1Gobjects (with extra tuning information)
  typedef vector<L1GObject> L1GObjectCollection;
  typedef std::auto_ptr<L1GObjectCollection> L1GObjectCollectionPtr;

  explicit UCTStage1BProducer(const edm::ParameterSet&);
  ~UCTStage1BProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // Helper methods

  void makeEGs();
  void makeTaus();

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

  Handle<L1CaloRegionCollection> newRegions;
  Handle<L1CaloRegionCollection> newEMRegions;
  Handle<L1CaloEmCollection> tauCands;
  Handle<L1GObjectCollection> emClusters;

  list<L1GObject> rlxEGList;
  list<L1GObject> isoEGList;

  list<L1GObject> rlxTauList;
  list<L1GObject> isoTauList;

};

//
// constructors and destructor
//
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
  produces<L1GObjectCollection>( "RelaxedEGUnpacked" ) ;
  produces<L1GObjectCollection>( "IsolatedEGUnpacked" ) ;
  produces<L1GObjectCollection>( "RelaxedTauUnpacked" ) ;
  produces<L1GObjectCollection>( "IsolatedTauUnpacked" ) ;

}


UCTStage1BProducer::~UCTStage1BProducer()
{
}


//
// member functions
//

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

  L1GObjectCollectionPtr unpackedRlxTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedRlxEGs(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoEGs(new L1GObjectCollection);

  for(list<L1GObject>::iterator rlxTau = rlxTauList.begin();
      rlxTau != rlxTauList.end();
      rlxTau++) {
    unpackedRlxTaus->push_back(*rlxTau);
  }
  for(list<L1GObject>::iterator isoTau = isoTauList.begin();
      isoTau != isoTauList.end();
      isoTau++) {
    unpackedIsoTaus->push_back(*isoTau);
  }
  for(list<L1GObject>::iterator rlxEG = rlxEGList.begin();
      rlxEG != rlxEGList.end();
      rlxEG++) {
    unpackedRlxEGs->push_back(*rlxEG);
  }
  for(list<L1GObject>::iterator isoEG = isoEGList.begin();
      isoEG != isoEGList.end();
      isoEG++) {
    unpackedIsoEGs->push_back(*isoEG);
  }

  iEvent.put(unpackedRlxTaus, "RelaxedTauUnpacked");
  iEvent.put(unpackedIsoTaus, "IsolatedTauUnpacked");
  iEvent.put(unpackedRlxEGs, "RelaxedEGUnpacked");
  iEvent.put(unpackedIsoEGs, "IsolatedEGUnpacked");

}


// ------------ method called once each job just before starting event loop  ------------
void
UCTStage1BProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
UCTStage1BProducer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
UCTStage1BProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
UCTStage1BProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
UCTStage1BProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
UCTStage1BProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void UCTStage1BProducer::makeEGs() {
  rlxEGList.clear();
  isoEGList.clear();
  for(vector<L1GObject>::const_iterator emCluster = emClusters->begin();
      emCluster != emClusters->end();
      emCluster++) {
    if(emCluster->et() > egSeed) {
      double emClusterEt = emCluster->ptValue();
      unsigned emClusterEtInRegionScale = emClusterEt / regionLSB_;
      unsigned emClusterRegionIPhi = emCluster->phiIndex();
      unsigned emClusterRegionIEta = emCluster->etaIndex();
      unsigned C = 0;
      unsigned N = 0;
      unsigned S = 0;
      unsigned E = 0;
      unsigned W = 0;
      unsigned NE = 0;
      unsigned SE = 0;
      unsigned NW = 0;
      unsigned SW = 0;
      for(L1CaloRegionCollection::const_iterator region = newRegions->begin();
	  region != newRegions->end(); region++) {
	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = region->et();
	}
      }
      // Compare emCluster to the total E+H within the region to define rgnIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      unsigned rgnIsolation = 0;
      if(C > emClusterEtInRegionScale) rgnIsolation = C - emClusterEtInRegionScale;
      unsigned jetIsolation = C + N + S + E + W + NE + NW + SE + SW - emClusterEtInRegionScale;
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
	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = region->et();
	}
      }
      // Compare emCluster to the total E within the region to define rgnEMIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      unsigned rgnEMIsolation = 0;
      if(C > emCluster->ptCode()) rgnEMIsolation = C - emCluster->ptCode();
      unsigned jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - emCluster->ptCode();
      // Every emCluster which passes HoECut is an egObject
      double h = ((double) C * regionLSB_) - ((double) emCluster->ptCode() * egLSB_);
      double e = ((double) emCluster->ptCode() * egLSB_);
      if((h / e) < regionalHoECut) {
	rlxEGList.push_back(L1GObject(emCluster->ptCode(), emCluster->etaTwrIndex(), emCluster->phiTwrIndex(), "RelaxedEG", true));
	// egObject that passes all isolations is an isolated egObject
	double relativeRgnEMIsolation = (double) rgnEMIsolation / (double) emCluster->ptCode();
	double relativeJetEMIsolation = (double) jetEMIsolation / (double) emCluster->ptCode();
	double relativeRgnIsolation = ((double) rgnIsolation * regionLSB_) / ((double) emCluster->ptCode() * egLSB_);
	double relativeJetIsolation = ((double) jetIsolation * regionLSB_) / ((double) emCluster->ptCode() * egLSB_);
	if(relativeRgnIsolation < egRelativeRgnIsolationCut &&
	   relativeJetIsolation < egRelativeJetIsolationCut &&
	   relativeRgnEMIsolation < egRelativeEMRgnIsolationCut &&
	   relativeJetEMIsolation < egRelativeEMJetIsolationCut) {
	  isoEGList.push_back(L1GObject(emCluster->ptCode(), emCluster->etaTwrIndex(), emCluster->phiTwrIndex(), "IsolatedEG", true));
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
    if(tauCand->rank() > tauSeed) {
      double tauCandEt = tauCand->rank() * tauLSB_;
      unsigned tauCandEtInRegionScale = tauCandEt / regionLSB_;
      unsigned tauCandEtCode = tauCandEtInRegionScale;
      unsigned tauCandRegionIPhi = tauCand->regionId().iphi();
      unsigned tauCandRegionIEta = tauCand->regionId().ieta();
      unsigned C = 0;
      unsigned N = 0;
      unsigned S = 0;
      unsigned E = 0;
      unsigned W = 0;
      unsigned NE = 0;
      unsigned SE = 0;
      unsigned NW = 0;
      unsigned SW = 0;
      for(L1CaloRegionCollection::const_iterator region = newRegions->begin();
	  region != newRegions->end(); region++) {
	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = region->et();
	  // If the ET in the region is larger than the addition to the 2x1
	  // we can use the region ET instead as the tauCandEt
	  if(C > tauCandEtCode) tauCandEtCode = C;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = region->et();
	}
      }
      // Compare tauCand to the total E+H within the region to define rgnIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      unsigned rgnIsolation = 0;
      if(C > tauCandEtCode) rgnIsolation = C - tauCandEtCode;
      unsigned jetIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEtCode;
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
	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = region->et();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = region->et();
	}
      }
      // Compare tauCand to the total E within the region to define rgnEMIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      unsigned rgnEMIsolation = 0;
      if(C > tauCandEtInRegionScale) rgnEMIsolation = C - tauCandEtCode;
      unsigned jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEtCode;
      // Use finer grain position resolution if possible using the emClusters
      // with default value being center of the region
      // phi: 0-17 becomes 0-71; eta: 0-13 becomes 0-55;
      unsigned tauCandIEta = (tauCandRegionIEta - 4) * 4;
      unsigned tauCandIPhi = tauCandRegionIPhi * 4 + 1;
      // FIXME make this next for loop work.
      for(vector<L1GObject>::const_iterator emCluster = emClusters->begin();
	  emCluster != emClusters->end();
	  emCluster++) {
	if(emCluster->phiIndex() == tauCandIPhi && emCluster->etaIndex() == tauCandIEta) {
	  tauCandIPhi = emCluster->phiTwrIndex();
	  tauCandIEta = emCluster->etaTwrIndex();
	}
      }
      rlxTauList.push_back(L1GObject(tauCandEtCode, tauCandIEta, tauCandIPhi, "RelaxedTau", true));
      // tauObject that passes all isolations is an isolated tauObject
      double relativeRgnEMIsolation = ((double) rgnEMIsolation * egLSB_) / ((double) tauCandEtCode * regionLSB_);
      double relativeJetEMIsolation = ((double) jetEMIsolation * egLSB_) / ((double) tauCandEtCode * regionLSB_);
      double relativeRgnIsolation = ((double) rgnIsolation) / ((double) tauCandEtCode);
      double relativeJetIsolation = ((double) jetIsolation) / ((double) tauCandEtCode);
      if(relativeRgnIsolation < tauRelativeRgnIsolationCut &&
	 relativeJetIsolation < tauRelativeJetIsolationCut &&
	 relativeRgnEMIsolation < tauRelativeEMRgnIsolationCut &&
	 relativeJetEMIsolation < tauRelativeEMJetIsolationCut) {
	isoTauList.push_back(L1GObject(tauCandEtCode, tauCandIEta, tauCandIPhi, "IsolatedTau", true));
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
