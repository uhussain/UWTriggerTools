// -*- C++ -*-
//
// Package:    UCT2015EClusterProducer
// Class:      UCT2015EClusterProducer
//
/**\class UCT2015EClusterProducer UCT2015EClusterProducer.cc L1Trigger/UCT2015/src/UCT2015EClusterProducer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Thu Jun  7 13:29:52 CDT 2012
// $Id: UCT2015EClusterProducer.cc,v 1.7 2013/02/22 17:08:09 friis Exp $
//
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
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"

#include "L1Trigger/UCT2015/src/L1GObject.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

//
// class declaration
//
//

class UCT2015EClusterProducer : public edm::EDProducer {
public:

  static const unsigned int N_TOWER_PHI;
  static const unsigned int N_TOWER_ETA;

  // A packed uint = pt, eta, and phi packed into 32 bits
  typedef vector<unsigned int> PackedUIntCollection;
  typedef std::auto_ptr<PackedUIntCollection> PackedUIntCollectionPtr;

  // Concrete collection of output objects (with extra tuning information)
  typedef vector<L1GObject> L1GObjectCollection;
  typedef std::auto_ptr<L1GObjectCollection> L1GObjectCollectionPtr;

  explicit UCT2015EClusterProducer(const edm::ParameterSet&);
  ~UCT2015EClusterProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // Helper methods

  void puSubtraction();
  void makeERegions();
  void makeEClusters();
  void printEClusters();

  // ----------member data ---------------------------

  bool debug;
  bool puCorrect;
  unsigned int puETMax;
  unsigned int eClusterSeed;
  double eLSB_;

  std::vector<edm::InputTag> ecalDigis;

  vector<vector <unsigned int> > eTowerETCode;
  vector<vector <bool> > eTowerFGVeto;

  unsigned int puLevel;

  list<L1GObject> eClusterList;
  L1CaloRegionCollection eRegionList;

};

//
// constants, enums and typedefs
//

// static data member definitions
//

unsigned int const UCT2015EClusterProducer::N_TOWER_PHI = 72;
unsigned int const UCT2015EClusterProducer::N_TOWER_ETA = 56;

//
// constructors and destructor
//
UCT2015EClusterProducer::UCT2015EClusterProducer(const edm::ParameterSet& iConfig) :
  debug(iConfig.getParameter<bool>("debug")),
  puCorrect(iConfig.getParameter<bool>("puCorrect")),
  puETMax(iConfig.getParameter<unsigned int>("puETMax")),
  eClusterSeed(iConfig.getParameter<unsigned int>("eClusterSeed")),
  eLSB_(iConfig.getParameter<double>("ecalLSB")),
  ecalDigis(iConfig.getParameter<std::vector<edm::InputTag> >("ecalDigis")),
  eTowerETCode(N_TOWER_PHI, vector<unsigned int>(N_TOWER_ETA)),
  eTowerFGVeto(N_TOWER_PHI, vector<bool>(N_TOWER_ETA))
{
  produces<PackedUIntCollection>( "PULevel" ) ;
  produces<PackedUIntCollection>( "EClusters" ) ;
  produces<L1GObjectCollection>( "EClustersUnpacked" ) ;
  produces<L1CaloRegionCollection>( "ERegions" );
}


UCT2015EClusterProducer::~UCT2015EClusterProducer()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
UCT2015EClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  puLevel = 0;

  PackedUIntCollectionPtr packedPULevels(new PackedUIntCollection);
  PackedUIntCollectionPtr packedEClusters(new PackedUIntCollection);
  L1GObjectCollectionPtr unpackedEClusters(new L1GObjectCollection);
  std::auto_ptr<L1CaloRegionCollection> ERegions (new L1CaloRegionCollection);

  edm::Handle<EcalTrigPrimDigiCollection> ecal;
  iEvent.getByLabel(ecalDigis[0], ecal);
  EcalTrigPrimDigiCollection ecalCollection;
  if (ecal.isValid()) { ecalCollection = *ecal;}
  else {return;}
  unsigned int nEcalDigi = ecalCollection.size();
  if (nEcalDigi>4032) {nEcalDigi=4032;}
  for (unsigned int i = 0; i < nEcalDigi; i++){
    int ieta = ecalCollection[i].id().ieta();
    // Note TPG absIeta counts from 1-28 (not 0-27)
    int cal_iphi = ecalCollection[i].id().iphi();
    unsigned int iEta;
    // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
    // And we have mapped [-28, -1], [1, 28] onto [0, 55]
    if(ieta < 0)
      iEta = ieta + 28;
    else if(ieta > 0)
      iEta = ieta + 27;
    else {
      std::cout << "UCT2015EClusterProducer: Illegal value for ecalCollection->id().ieta() -- should never be zero!" << std::endl;
      return;
    }

    // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
    unsigned int iPhi = (cal_iphi-1);
    eTowerETCode[iPhi][iEta] = ecalCollection[i].compressedEt();
    eTowerFGVeto[iPhi][iEta] = (bool) ecalCollection[i].fineGrain();     // 0 or 1
  }

  if(puCorrect) puSubtraction();
  packedPULevels->push_back(puLevel);
  makeERegions();
  makeEClusters();
  if(debug) printEClusters();

  for(list<L1GObject>::iterator eCluster = eClusterList.begin();
      eCluster != eClusterList.end();
      eCluster++) {
    packedEClusters->push_back(eCluster->packedObject());
    unpackedEClusters->push_back(*eCluster);
  }

  iEvent.put(packedPULevels, "PULevel");
  iEvent.put(packedEClusters, "EClusters");
  iEvent.put(unpackedEClusters, "EClustersUnpacked");
  iEvent.put(ERegions, "ERegions");
}


// ------------ method called once each job just before starting event loop  ------------
void
UCT2015EClusterProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
UCT2015EClusterProducer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
UCT2015EClusterProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
UCT2015EClusterProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
UCT2015EClusterProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
UCT2015EClusterProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// NB PU is not in the physical scale!!  Needs to be multiplied by regionLSB
void UCT2015EClusterProducer::puSubtraction()
{
  puLevel = 0;
  unsigned int puCount = 0;
  for(unsigned int iPhi = 0; iPhi < N_TOWER_PHI; iPhi++) {
    for(unsigned int iEta = 0; iEta < N_TOWER_ETA; iEta++) {
      if(eTowerETCode[iPhi][iEta] <= puETMax) {
	puLevel += eTowerETCode[iPhi][iEta]; puCount++;
      }
    }
  }
  // Add a factor of 12x12, so it corresponds to a jet.  Reduces roundoff error.
  puLevel *= 144;
  puLevel = puLevel / puCount;
}

void UCT2015EClusterProducer::makeERegions() {
  eRegionList.clear();
  for(unsigned int iPhi = 0; iPhi < N_TOWER_PHI/4; iPhi++) {
    for(unsigned int iEta = 0; iEta < N_TOWER_ETA/4; iEta++) {
      unsigned int et = 0;
      bool overflow = false;
      bool finegrain = false;
      bool mip = false;
      bool quiet = false;
      for(unsigned int i = 0; i < 4; i++) {
	for(unsigned int j = 0; j < 4; j++) {
	  et += eTowerETCode[iPhi * 4 + i][iEta * 4 + j];
	}
      }
      if(et > 0xFFFF) overflow = true;
      if(et > 0x0001 && et < 0x0004) mip = true;
      if(et < 0x0002) quiet = true;
      eRegionList.push_back(L1CaloRegion(et, finegrain, overflow, mip, quiet, iEta, iPhi));
    }
  }
}

void UCT2015EClusterProducer::makeEClusters() {
  eClusterList.clear();
  for(unsigned int iPhi = 0; iPhi < N_TOWER_PHI; iPhi++) {
    for(unsigned int iEta = 0; iEta < N_TOWER_ETA; iEta++) {
      if(eTowerETCode[iPhi][iEta] > eClusterSeed) {
	unsigned int center_et = eTowerETCode[iPhi][iEta];
	unsigned int neighborN_et = 0;
	unsigned int neighborS_et = 0;
	unsigned int neighborE_et = 0;
	unsigned int neighborW_et = 0;
	unsigned int neighborNE_et = 0;
	unsigned int neighborNW_et = 0;
	unsigned int neighborSE_et = 0;
	unsigned int neighborSW_et = 0;
	unsigned int N;
	unsigned int S;
	if(iPhi > 0 && iPhi < (N_TOWER_PHI - 1)) {
	  N = iPhi + 1;
	  S = iPhi - 1;
	}
	else if(iPhi == (N_TOWER_PHI - 1)) {
	  N = 0;
	  S = iPhi - 1;
	}
	else if(iPhi == 0) {
	  N = iPhi + 1;
	  S = N_TOWER_PHI - 1;
	}
	if(iEta == 0) {
	  unsigned int E = iEta + 1;
	  neighborN_et = eTowerETCode[N][iEta];
	  neighborS_et = eTowerETCode[S][iEta];
	  neighborE_et = eTowerETCode[iPhi][E];
	  neighborNE_et = eTowerETCode[N][E];
	  neighborSE_et = eTowerETCode[S][E];
	}
	else if(iEta == N_TOWER_ETA - 1) {
	  unsigned int W = iEta - 1;
	  neighborN_et = eTowerETCode[N][iEta];
	  neighborS_et = eTowerETCode[S][iEta];
	  neighborW_et = eTowerETCode[iPhi][W];
	  neighborNW_et = eTowerETCode[N][W];
	  neighborSW_et = eTowerETCode[S][W];
	}
	else {
	  unsigned int E = iEta + 1;
	  unsigned int W = iEta - 1;
	  neighborN_et = eTowerETCode[N][iEta];
	  neighborS_et = eTowerETCode[S][iEta];
	  neighborE_et = eTowerETCode[iPhi][E];
	  neighborW_et = eTowerETCode[iPhi][W];
	  neighborNE_et = eTowerETCode[N][E];
	  neighborNW_et = eTowerETCode[N][W];
	  neighborSE_et = eTowerETCode[S][E];
	  neighborSW_et = eTowerETCode[S][W];
	}
	if(center_et > neighborN_et &&
	   center_et > neighborNW_et &&
	   center_et > neighborW_et &&
	   center_et > neighborSW_et &&
	   center_et >= neighborNE_et &&
	   center_et >= neighborE_et &&
	   center_et >= neighborSE_et &&
	   center_et >= neighborS_et) {
	  unsigned int eClusterET = center_et +
	    neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	    neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;
	  // Temporarily use the tower (iPhi, iEta) -- todo: convert to half-tower resolution
	  unsigned int eClusterPhi = iPhi;
	  unsigned int eClusterEta = iEta;
	  eClusterList.push_back(L1GObject(eClusterET, eClusterEta, eClusterPhi, "ECluster", true));
	}
      }
    }
  }
  eClusterList.sort();
  eClusterList.reverse();
}

void UCT2015EClusterProducer::printEClusters() {
  std::cout << "eClusterList.size() = " << eClusterList.size() << std::endl;
  for(list<L1GObject>::iterator eCluster = eClusterList.begin();
      eCluster != eClusterList.end();
      eCluster++) {
    std::cout << *eCluster << std::endl;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(UCT2015EClusterProducer);
