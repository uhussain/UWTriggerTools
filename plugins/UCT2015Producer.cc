// -*- C++ -*-
//
// Package:    UCT2015Producer
// Class:      UCT2015Producer
//
/**\class UCT2015Producer UCT2015Producer.cc L1Trigger/UCT2015/src/UCT2015Producer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Thu Jun  7 13:29:52 CDT 2012
// $Id: UCT2015Producer.cc,v 1.20 2013/01/17 11:57:51 friis Exp $
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
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "L1Trigger/UCT2015/src/L1GObject.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

//
// class declaration
//
//

// mapping GCT eta to physical eta
const double _etaValues[11] = {
 0.174, 0.522, 0.870, 1.218, 1.566, 1.956, 2.586, 3.250, 3.750, 4.250, 4.750
};

double etaValue(int myEta) {
  if(myEta < 11) {
    return -_etaValues[-(myEta - 10)]; // 0-10 are negative eta values
  }
  else if (myEta < 22) {
    return _etaValues[myEta - 11];     // 11-21 are positive eta values
  }
  return 999.;
}

class UCT2015Producer : public edm::EDProducer {
public:

  static const unsigned N_JET_PHI;
  static const unsigned N_JET_ETA;

  // A packed uint = pt, eta, and phi packed into 32 bits
  typedef vector<unsigned int> PackedUIntCollection;
  typedef std::auto_ptr<PackedUIntCollection> PackedUIntCollectionPtr;

  // Concrete collection of L1Gobjects (with extra tuning information)
  typedef vector<L1GObject> L1GObjectCollection;
  typedef std::auto_ptr<L1GObjectCollection> L1GObjectCollectionPtr;

  explicit UCT2015Producer(const edm::ParameterSet&);
  ~UCT2015Producer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  double egPhysicalEt(const L1CaloEmCand& cand) {
    return egLSB_*cand.rank();
  }

  double regionPhysicalEt(const L1CaloRegion& cand) {
    return regionLSB_*cand.et();
  }

  double jetcorr(const double ptraw, const int ieta, const uint pu) {

    double calib[] = {  // calib_pudefault
      0.961508, -0.961508, 15.947,
      1.36807, -1.36807, 6.68028,
      1.44409, -1.44409, 8.29767,
      1.51892, -1.51892, 13.5637,
      1.24288, -1.24288, 23.8283,
      1.19104, -1.19104, 29.0307,
      1.22199, -1.22199, 30.0691,
      1.22154, -1.22154, 30.5915,
      1.20026, -1.20026, 28.6067,
      1.15149, -1.15149, 27.4127,
      1.16534, -1.16534, 25.4566,
      1.16308, -1.16308, 25.8007,
      1.19049, -1.19049, 27.2719,
      1.21544, -1.21544, 28.4664,
      1.24526, -1.24526, 30.216,
      1.20326, -1.20326, 30.2056,
      1.17776, -1.17776, 29.7033,
      1.23713, -1.23713, 23.7761,
      1.49303, -1.49303, 14.181,
      1.46926, -1.46926, 7.82764,
      1.32482, -1.32482, 7.18464,
      0.974094, -0.974094, 15.9478
    };
    //==============================================

    double alpha = calib[3*ieta + 0];
    double beta  = calib[3*ieta + 1];
    double gamma = calib[3*ieta + 2];

    double pt =  alpha * ptraw + beta * pu + gamma;
    return pt;
  }

  double jetcorrUIC(const double ptraw, const int ieta, const uint pu) {

    double calib[] = {  // calib_puUIC
      1.04423, -1.09017, 13.5378,
      1.36887, -2.14366, 8.37627,
      1.43954, -2.25431, 10.1702,
      1.48093, -2.82625, 18.5013,
      1.2272, -2.25491, 27.179,
      1.18273, -2.08923, 31.7074,
      1.22839, -1.44659, 29.0857,
      1.22846, -1.33895, 29.088,
      1.20654, -1.31506, 27.1161,
      1.15725, -1.26133, 25.9607,
      1.17086, -1.27616, 23.9877,
      1.16862, -1.27373, 24.3519,
      1.19648, -1.30409, 25.7887,
      1.22124, -1.33107, 26.9854,
      1.25245, -1.36509, 28.6541,
      1.20829, -1.42292, 29.2688,
      1.16759, -2.06248, 32.3936,
      1.21936, -2.24049, 27.1477,
      1.44564, -2.75891, 19.1133,
      1.46288, -2.29087, 9.76789,
      1.3305, -2.08356, 8.76516,
      1.04901, -1.09517, 13.6107
    };
    //==============================================

    // std::cout << " JetCorr: " << ptraw << " " << ieta << " " << pu << std::endl;
    double alpha = calib[3*ieta + 0];
    double beta  = calib[3*ieta + 1];
    double gamma = calib[3*ieta + 2];

    double pt =  alpha * ptraw + beta * pu + gamma;
    return pt;
  }

  // Helper methods

  void puSubtraction();

  void makeSums();
  void makeJets();
  void makeEGTaus();

  list<L1GObject> correctJets(list<L1GObject>,string);

  // ----------member data ---------------------------

  bool puCorrect;
  bool useUICrho; // which PU denstity to use for energy correction determination

  unsigned int puETMax;
  unsigned int puLevel;
  //double puLevelUIC; // puLevel divided by puCount*Area, not multiply by 9.0
  unsigned int  puLevelUIC; // puLevel divided by puCount*Area, not multiply by 9.0

  unsigned int sumET;
  int sumEx;
  int sumEy;
  unsigned int MET;

  unsigned int regionETCutForHT;
  unsigned int regionETCutForMET;
  unsigned int minGctEtaForSums;
  unsigned int maxGctEtaForSums;
  unsigned int sumHT;
  int sumHx;
  int sumHy;
  unsigned int MHT;

  unsigned int sumExtraET;
  unsigned int extraMET;
  unsigned int sumExtraHT;
  unsigned int extraMHT;

  L1GObject METObject;
  L1GObject METSIGObject;
  L1GObject MHTObject;
  L1GObject SETObject;
  L1GObject SHTObject;

  unsigned int jetSeed;
  list<L1GObject> jetList, corrJetList;

  unsigned int egtSeed;
  double relativeIsolationCut;
  double relativeJetIsolationCut;
  list<L1GObject> rlxTauList, corrRlxTauList;
  list<L1GObject> rlxEGList;
  list<L1GObject> isoTauList, corrIsoTauList;
  list<L1GObject> isoEGList;

  Handle<L1CaloRegionCollection> newRegions;
  Handle<L1CaloEmCollection> newEMCands;

  vector<double> sinPhi;
  vector<double> cosPhi;

  double egLSB_;
  double regionLSB_;

};

//
// constants, enums and typedefs
//

// static data member definitions
//

unsigned const UCT2015Producer::N_JET_PHI = L1CaloRegionDetId::N_PHI * 4;
unsigned const UCT2015Producer::N_JET_ETA = L1CaloRegionDetId::N_ETA * 4;

//
// constructors and destructor
//
UCT2015Producer::UCT2015Producer(const edm::ParameterSet& iConfig) :
  puCorrect(iConfig.getParameter<bool>("puCorrect")),
  useUICrho(iConfig.getParameter<bool>("useUICrho")),
  puETMax(iConfig.getParameter<unsigned int>("puETMax")),
  regionETCutForHT(iConfig.getParameter<unsigned int>("regionETCutForHT")),
  regionETCutForMET(iConfig.getParameter<unsigned int>("regionETCutForMET")),
  minGctEtaForSums(iConfig.getParameter<unsigned int>("minGctEtaForSums")),
  maxGctEtaForSums(iConfig.getParameter<unsigned int>("maxGctEtaForSums")),
  jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
  egtSeed(iConfig.getParameter<unsigned int>("egtSeed")),
  relativeIsolationCut(iConfig.getParameter<double>("relativeIsolationCut")),
  relativeJetIsolationCut(iConfig.getParameter<double>("relativeJetIsolationCut")),
  egLSB_(iConfig.getParameter<double>("egammaLSB")),
  regionLSB_(iConfig.getParameter<double>("regionLSB"))
{
  puLevel = 0;
  puLevelUIC = 0.0;
  produces<PackedUIntCollection>( "PULevel" ) ;
  produces<PackedUIntCollection>( "PULevelUIC" ) ;
  produces<PackedUIntCollection>( "MET" ) ;
  produces<PackedUIntCollection>( "MHT" ) ;
  produces<PackedUIntCollection>( "SET" ) ;
  produces<PackedUIntCollection>( "SHT" ) ;
  produces<PackedUIntCollection>( "Jet" ) ;
  produces<PackedUIntCollection>( "CorrJet" ) ;
  produces<PackedUIntCollection>( "RelaxedEG" ) ;
  produces<PackedUIntCollection>( "IsolatedEG" ) ;
  produces<PackedUIntCollection>( "RelaxedTau" ) ;
  produces<PackedUIntCollection>( "IsolatedTau" ) ;
  produces<PackedUIntCollection>( "CorrRelaxedTau" ) ;
  produces<PackedUIntCollection>( "CorrIsolatedTau" ) ;

  // Also declare we produce unpacked collections (which have more info)
  produces<L1GObjectCollection>( "JetUnpacked" ) ;
  produces<L1GObjectCollection>( "CorrJetUnpacked" ) ;
  produces<L1GObjectCollection>( "RelaxedEGUnpacked" ) ;
  produces<L1GObjectCollection>( "IsolatedEGUnpacked" ) ;
  produces<L1GObjectCollection>( "RelaxedTauUnpacked" ) ;
  produces<L1GObjectCollection>( "IsolatedTauUnpacked" ) ;
  produces<L1GObjectCollection>( "CorrRelaxedTauUnpacked" ) ;
  produces<L1GObjectCollection>( "CorrIsolatedTauUnpacked" ) ;

  produces<L1GObjectCollection>( "PULevelUnpacked" ) ;
  produces<L1GObjectCollection>( "PULevelUICUnpacked" ) ;
  produces<L1GObjectCollection>( "METUnpacked" ) ;
  produces<L1GObjectCollection>( "METSIGUnpacked" ) ;
  produces<L1GObjectCollection>( "MHTUnpacked" ) ;
  produces<L1GObjectCollection>( "SETUnpacked" ) ;
  produces<L1GObjectCollection>( "SHTUnpacked" ) ;

  //now do what ever initialization is needed
  for(unsigned int i = 0; i < L1CaloRegionDetId::N_PHI; i++) {
    sinPhi.push_back(sin(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
    cosPhi.push_back(cos(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
  }
}


UCT2015Producer::~UCT2015Producer()
{
}


//
// member functions
//

// For the single objects, like MET/MHT, etc, convert them into a
// std::auto_ptr<L1GObjectCollection> suitable for putting into the edm::Event
// The "collection" contains only 1 object.
UCT2015Producer::L1GObjectCollectionPtr collectionize(const L1GObject& obj) {
  return UCT2015Producer::L1GObjectCollectionPtr(
      new UCT2015Producer::L1GObjectCollection(1, obj));
}

// ------------ method called for each event  ------------
void
UCT2015Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iEvent.getByLabel("uctDigis", newRegions);
  iEvent.getByLabel("uctDigis", newEMCands);

  if(puCorrect) puSubtraction();

  makeSums();
  makeJets();
  makeEGTaus();


  PackedUIntCollectionPtr packedPULevels(new PackedUIntCollection);
  PackedUIntCollectionPtr packedPULevelsUIC(new PackedUIntCollection);
  PackedUIntCollectionPtr packedMETs(new PackedUIntCollection);
  PackedUIntCollectionPtr packedMHTs(new PackedUIntCollection);
  PackedUIntCollectionPtr packedSETs(new PackedUIntCollection);
  PackedUIntCollectionPtr packedSHTs(new PackedUIntCollection);

  PackedUIntCollectionPtr packedJets(new PackedUIntCollection);
  PackedUIntCollectionPtr packedRlxTaus(new PackedUIntCollection);
  PackedUIntCollectionPtr packedIsoTaus(new PackedUIntCollection);
  PackedUIntCollectionPtr packedCorrJets(new PackedUIntCollection);
  PackedUIntCollectionPtr packedCorrRlxTaus(new PackedUIntCollection);
  PackedUIntCollectionPtr packedCorrIsoTaus(new PackedUIntCollection);
  PackedUIntCollectionPtr packedRlxEGs(new PackedUIntCollection);
  PackedUIntCollectionPtr packedIsoEGs(new PackedUIntCollection);

  L1GObjectCollectionPtr unpackedJets(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedRlxTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrJets(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrRlxTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrIsoTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedRlxEGs(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoEGs(new L1GObjectCollection);

  packedPULevels->push_back(puLevel);
  packedPULevelsUIC->push_back(puLevelUIC);
  packedMETs->push_back(METObject.packedObject());
  packedMHTs->push_back(MHTObject.packedObject());
  packedSETs->push_back(SETObject.packedObject());
  packedSHTs->push_back(SHTObject.packedObject());

  //uncorrected Jet and Tau collections
  for(list<L1GObject>::iterator jet = jetList.begin();
      jet != jetList.end();
      jet++) {
    packedJets->push_back(jet->packedObject());
    unpackedJets->push_back(*jet);
  }
  for(list<L1GObject>::iterator rlxTau = rlxTauList.begin();
      rlxTau != rlxTauList.end();
      rlxTau++) {
    packedRlxTaus->push_back(rlxTau->packedObject());
    unpackedRlxTaus->push_back(*rlxTau);
  }
  for(list<L1GObject>::iterator isoTau = isoTauList.begin();
      isoTau != isoTauList.end();
      isoTau++) {
    packedIsoTaus->push_back(isoTau->packedObject());
    unpackedIsoTaus->push_back(*isoTau);
  }
  //corrected Jet and Tau collections
  corrJetList=correctJets(jetList,"Jets");
  corrRlxTauList=correctJets(rlxTauList,"Tau");
  corrIsoTauList=correctJets(isoTauList,"IsoTau");

  for(list<L1GObject>::iterator jet = corrJetList.begin();
      jet != corrJetList.end();
      jet++) {
    packedCorrJets->push_back(jet->packedObject());
    unpackedCorrJets->push_back(*jet);
  }
  for(list<L1GObject>::iterator rlxTau = corrRlxTauList.begin();
      rlxTau != corrRlxTauList.end();
      rlxTau++) {
    packedCorrRlxTaus->push_back(rlxTau->packedObject());
    unpackedCorrRlxTaus->push_back(*rlxTau);
  }
  for(list<L1GObject>::iterator isoTau = corrIsoTauList.begin();
      isoTau != corrIsoTauList.end();
      isoTau++) {
    packedCorrIsoTaus->push_back(isoTau->packedObject());
    unpackedCorrIsoTaus->push_back(*isoTau);
  }

  // egamma collections
  for(list<L1GObject>::iterator rlxEG = rlxEGList.begin();
      rlxEG != rlxEGList.end();
      rlxEG++) {
    packedRlxEGs->push_back(rlxEG->packedObject());
    unpackedRlxEGs->push_back(*rlxEG);
  }
  for(list<L1GObject>::iterator isoEG = isoEGList.begin();
      isoEG != isoEGList.end();
      isoEG++) {
    packedIsoEGs->push_back(isoEG->packedObject());
    unpackedIsoEGs->push_back(*isoEG);
  }

  iEvent.put(collectionize(puLevel), "PULevelUnpacked");
  iEvent.put(collectionize(puLevelUIC), "PULevelUICUnpacked");
  iEvent.put(collectionize(METObject), "METUnpacked");
  iEvent.put(collectionize(MHTObject), "MHTUnpacked");
  iEvent.put(collectionize(SETObject), "SETUnpacked");
  iEvent.put(collectionize(SHTObject), "SHTUnpacked");
  iEvent.put(collectionize(METSIGObject), "METSIGUnpacked");

  iEvent.put(packedJets, "Jet");
  iEvent.put(packedRlxTaus, "RelaxedTau");
  iEvent.put(packedIsoTaus, "IsolatedTau");
  iEvent.put(packedCorrJets, "CorrJet");
  iEvent.put(packedCorrRlxTaus, "CorrRelaxedTau");
  iEvent.put(packedCorrIsoTaus, "CorrIsolatedTau");
  iEvent.put(packedRlxEGs, "RelaxedEG");
  iEvent.put(packedIsoEGs, "IsolatedEG");

  iEvent.put(unpackedJets, "JetUnpacked");
  iEvent.put(unpackedRlxTaus, "RelaxedTauUnpacked");
  iEvent.put(unpackedIsoTaus, "IsolatedTauUnpacked");
  iEvent.put(unpackedCorrJets, "CorrJetUnpacked");
  iEvent.put(unpackedCorrRlxTaus, "CorrRelaxedTauUnpacked");
  iEvent.put(unpackedCorrIsoTaus, "CorrIsolatedTauUnpacked");
  iEvent.put(unpackedRlxEGs, "RelaxedEGUnpacked");
  iEvent.put(unpackedIsoEGs, "IsolatedEGUnpacked");

  iEvent.put(packedPULevels, "PULevel");
  iEvent.put(packedPULevelsUIC, "PULevelUIC");
  iEvent.put(packedMETs, "MET");
  iEvent.put(packedMHTs, "MHT");
  iEvent.put(packedSETs, "SET");
  iEvent.put(packedSHTs, "SHT");
}


// ------------ method called once each job just before starting event loop  ------------
void
UCT2015Producer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
UCT2015Producer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
UCT2015Producer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
UCT2015Producer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
UCT2015Producer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
UCT2015Producer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}




// NB PU is not in the physical scale!!  Needs to be multiplied by regionLSB
void UCT2015Producer::puSubtraction()
{
  puLevel = 0;
  puLevelUIC = 0;
  double r_puLevelUIC=0.0;

  int puCount = 0;
  double Rarea=0.0;
  for(L1CaloRegionCollection::const_iterator newRegion =
	newRegions->begin();
      newRegion != newRegions->end(); newRegion++){
    if(regionPhysicalEt(*newRegion) <= puETMax) {
      puLevel += newRegion->et(); puCount++;
      r_puLevelUIC += newRegion->et();
      if (newRegion->gctEta()==0) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==1) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==2) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==3) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==4) Rarea += 0.828*0.348;
      if (newRegion->gctEta()==5) Rarea += 0.432*0.348;
      if (newRegion->gctEta()==6) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==7) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==8) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==9) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==10) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==21) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==20) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==19) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==18) Rarea += 0.5*0.348;
      if (newRegion->gctEta()==17) Rarea += 0.828*0.348;
      if (newRegion->gctEta()==16) Rarea += 0.432*0.348;
      if (newRegion->gctEta()==15) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==14) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==13) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==12) Rarea += 0.348*0.348;
      if (newRegion->gctEta()==11) Rarea += 0.348*0.348;
    }
  }
  // Add a factor of 9, so it corresponds to a jet.  Reduces roundoff error.
  puLevel *= 9;
  // std::cout << puLevel << " : " << puCount << " : " << puLevel/puCount << std::endl;
  puLevel = puLevel / puCount;
  r_puLevelUIC = r_puLevelUIC / Rarea;
  puLevelUIC=0;
  if (r_puLevelUIC > 0.) puLevelUIC = floor (r_puLevelUIC + 0.5);

}

// Resolutions taken from:
// http://cmslxr.fnal.gov/lxr/source/RecoMET/METProducers/python/CaloMETSignif_cfi.py

namespace uct_metsig {

double metSigEcalBarrelResolution(double et) {
  const double par[3] = {0.2,0.03,0.005};
  return et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
}

double metSigEcalEndCapResolution(double et) {
  const double par[3] = {0.2,0.03,0.005};
  return et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
}

double metSigHcalBarrelResolution(double et) {
  const double par[3] = {0.,1.22,0.05};
  return et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
}

double metSigHcalEndCapResolution(double et) {
  const double par[3] = {0.,1.3,0.05};
  return et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
}

double getEtResolution(bool isECAL, double et, int ieta) {
  if (et == 0)
    return 0;
  if (isECAL) {
    // ECAL like, use ECAL resolution
    if (std::abs(etaValue(ieta)) > 1.5) {
      return metSigEcalEndCapResolution(et);
    } else {
      return metSigEcalBarrelResolution(et);
    }
  } else {
    if (std::abs(etaValue(ieta)) > 1.5) {
      return metSigHcalEndCapResolution(et);
    } else {
      return metSigHcalBarrelResolution(et);
    }
  }
}

} // end uct_metsig (met significance functions).

void UCT2015Producer::makeSums()
{
  sumET = 0;
  sumEx = 0;
  sumEy = 0;
  sumHT = 0;
  sumHx = 0;
  sumHy = 0;

  const double METSIG_LSB = 1./128.;
  // The MET covariance matrix.
  int cov00 = 0;
  int cov01 = 0;
  int cov10 = 0;
  int cov11 = 0;

  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++){
    // Remove forward stuff
    if (newRegion->gctEta() < minGctEtaForSums || newRegion->gctEta() > maxGctEtaForSums) {
      continue;
    }
    //unsigned int regionET = newRegion->et() - puLevel;
    double regionET = std::max(regionPhysicalEt(*newRegion) - puLevel*regionLSB_/9., 0.);
    // Decide if it is ECAL-like or HCAL-like
    bool isECAL = !newRegion->tauVeto() && !newRegion->mip();
    double sigma_et = uct_metsig::getEtResolution(isECAL, regionET, newRegion->gctEta());
    double sigmaPhi = 0.35;

    double sigma0_2=sigma_et*sigma_et;
    double sigma1_2=sigmaPhi*sigmaPhi*regionET*regionET;

    double cosphi = cosPhi[newRegion->gctPhi()];
    double sinphi = sinPhi[newRegion->gctPhi()];

//    std::cout << "et: " << regionET << " sigma: " << sigma_et
//      << " sigma02: " << sigma0_2 << " sigma12: " << sigma1_2
//      << " cosphi: " << cosphi << " sinphi: " << sinphi << std::endl;

    // Not sure how this will translate to a FPGA.  Lets give it a LSB of
    // 0.125 GeV - these numbers should be small.
    cov00 += (int) ((sigma0_2*cosphi*cosphi + sigma1_2*sinphi*sinphi) / METSIG_LSB);
    cov01 += (int) ((cosphi*sinphi*(sigma0_2 - sigma1_2)) / METSIG_LSB);
    cov10 += (int) ((cosphi*sinphi*(sigma0_2 - sigma1_2)) / METSIG_LSB);
    cov11 += (int) ((sigma1_2*cosphi*cosphi + sigma0_2*sinphi*sinphi) / METSIG_LSB);
    if(regionET >= regionETCutForMET){
    	sumET += regionET;
    	sumEx += (int) (((double) regionET) * cosPhi[newRegion->gctPhi()]);
    	sumEy += (int) (((double) regionET) * sinPhi[newRegion->gctPhi()]);
    }
    if(regionET >= regionETCutForHT) {
      sumHT += regionET;
      sumHx += (int) (((double) regionET) * cosPhi[newRegion->gctPhi()]);
      sumHy += (int) (((double) regionET) * sinPhi[newRegion->gctPhi()]);
    }
  }
  MET = ((unsigned int) sqrt(sumEx * sumEx + sumEy * sumEy));
  MHT = ((unsigned int) sqrt(sumHx * sumHx + sumHy * sumHy));
  unsigned int iPhi = L1CaloRegionDetId::N_PHI * (atan2(sumEy, sumEx) + 3.1415927) / (2 * 3.1415927);
  METObject = L1GObject(MET, 0, iPhi, "MET");
  iPhi = L1CaloRegionDetId::N_PHI * (atan2(sumHy, sumHx) + 3.1415927) / (2 * 3.1415927);
  MHTObject = L1GObject(MHT, 0, iPhi, "MHT");
  SETObject = L1GObject(sumET, 0, 0, "SumET");
  SHTObject = L1GObject(sumHT, 0, 0, "SumHT");

  // Significance of MET = Transpose(met).Inverse(metCov).met
  // One factor of METSIG_LSB cancels out.
  int significance = (1./METSIG_LSB)*(
      -(cov11*sumEx*sumEx) + sumEy*((cov01 + cov10)*sumEx - cov00*sumEy)
      )/(cov01*cov10 - cov00*cov11);

//  std::cout << "covariance:" << std::endl;
//  std::cout << cov00 << " " << cov01 << std::endl;
//  std::cout << cov10 << " " << cov11 << std::endl;
//  std::cout << "ex: " << sumEx << " ey: " << sumEy << std::endl;
//  std::cout << "sig: " << significance << std::endl;

  METSIGObject = L1GObject(significance, 0, iPhi, "MET");
}

int deltaPhi18(int phi1, int phi2) {
  // Compute the difference in phi between two towers, wrapping at phi = 18
  int difference = phi1 - phi2;
  if (std::abs(phi1 - phi2) == 17) {
    difference = -difference/std::abs(difference);
  }
  return difference;
}

int deltaGctPhi(const L1CaloRegion& r1, const L1CaloRegion& r2) {
  return deltaPhi18(r1.gctPhi(), r2.gctPhi());
}

void debug(const std::string& dir, const L1CaloRegion& r1, const L1CaloRegion& r2) {
  //std::cout << "Found " << dir << " neighbor of (" << r1.gctPhi() << "," << r1.gctEta() << ") at (" << r2.gctPhi() << "," << r2.gctEta() << ")" << std::endl;
}

void UCT2015Producer::makeJets() {
  jetList.clear();
  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++) {
    double regionET = regionPhysicalEt(*newRegion);
    if(regionET > jetSeed) {
      double neighborN_et = 0;
      double neighborS_et = 0;
      double neighborE_et = 0;
      double neighborW_et = 0;
      double neighborNE_et = 0;
      double neighborSW_et = 0;
      double neighborNW_et = 0;
      double neighborSE_et = 0;
      unsigned int nNeighbors = 0;
      //std::cout << "Looking for seed @ " << newRegion->gctPhi() << " " << newRegion->gctEta() << std::endl;
      for(L1CaloRegionCollection::const_iterator neighbor = newRegions->begin();
	  neighbor != newRegions->end(); neighbor++) {
	double neighborET = regionPhysicalEt(*neighbor);
	if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	   (newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborN_et = neighborET;
          nNeighbors++;
          debug("N", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborS_et = neighborET;
          nNeighbors++;
          debug("S", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborE_et = neighborET;
          nNeighbors++;
          debug("E", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborW_et = neighborET;
          nNeighbors++;
          debug("W", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborNE_et = neighborET;
          nNeighbors++;
          debug("NE", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborSW_et = neighborET;
          nNeighbors++;
          debug("SW", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborNW_et = neighborET;
          nNeighbors++;
          debug("NW", *newRegion, *neighbor);
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborSE_et = neighborET;
          nNeighbors++;
          debug("SE", *newRegion, *neighbor);
	  continue;
	}
      }
      if(regionET > neighborN_et &&
	 regionET > neighborNW_et &&
	 regionET > neighborW_et &&
	 regionET > neighborSW_et &&
	 regionET >= neighborNE_et &&
	 regionET >= neighborE_et &&
	 regionET >= neighborSE_et &&
	 regionET >= neighborS_et) {
	unsigned int jetET = regionET +
	  neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	  neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;
	/*
	  int jetPhi = newRegion->gctPhi() * 4 +
	  ( - 2 * (neighborS_et + neighborSE_et + neighborSW_et)
	  + 2 * (neighborN_et + neighborNE_et + neighborNW_et) ) / jetET;
	  if(jetPhi < 0) {

	  }
	  else if(jetPhi >= ((int) N_JET_PHI)) {
	  jetPhi -= N_JET_PHI;
	  }
	  int jetEta = newRegion->gctEta() * 4 +
	  ( - 2 * (neighborW_et + neighborNW_et + neighborSW_et)
	  + 2 * (neighborE_et + neighborNE_et + neighborSE_et) ) / jetET;
	  if(jetEta < 0) jetEta = 0;
	  if(jetEta >= ((int) N_JET_ETA)) jetEta = N_JET_ETA - 1;
	*/
	// Temporarily use the region granularity -- we will try to improve as above when code is debugged
	int jetPhi = newRegion->gctPhi();
	int jetEta = newRegion->gctEta();

        bool neighborCheck = (nNeighbors == 8);
        // On the eta edge we only expect 5 neighbors
        if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
          neighborCheck = true;

        if (!neighborCheck) {
          std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
          assert(false);
        }
	jetList.push_back(L1GObject(jetET, jetEta, jetPhi, "Jet"));
        // Embed the puLevel information in the jet object for later tuning
        jetList.back().puLevel_ = puLevel;
       jetList.back().puLevelUIC_ = puLevelUIC;
        // Store information about the "core" PT of the jet (central region)
        jetList.back().associatedRegionEt_ = regionET;
      }
    }
  }
  jetList.sort();
  jetList.reverse();
}

list<L1GObject> UCT2015Producer::correctJets(list<L1GObject> jets, string L1ObjName) {

  list<L1GObject> corrlist;
  if ( ! puCorrect ) return corrlist;  // jet corrections only valid if PU density has been calculated

  corrlist.clear();

  for(list<L1GObject>::iterator jet = jets.begin(); jet != jets.end(); jet++) {

    double jetET=jet->ptValue();
    double jpt;
    unsigned int corjetET=0;

    //ccla  apply Michael's jet correction function
    if (useUICrho){
      jpt=jetcorrUIC(jetET,jet->etaIndex(),jet->puLevelUIC());
    }else{
      jpt=jetcorr(jetET,jet->etaIndex(),jet->puLevel());
    }
    if (jpt>0) corjetET = floor (jpt + 0.5);

    // std::cout << L1ObjName << " " << jetET << " " << jpt << " " << corjetET << std::endl;
    corrlist.push_back(L1GObject(corjetET, jet->etaIndex(), jet->phiIndex(), L1ObjName));

    if (useUICrho){ // only store the PU density used for to obtain the JEC. Can be used to let user easily know which correction was used
      corrlist.back().puLevel_ = 0;
      corrlist.back().puLevelUIC_ = jet->puLevelUIC();
    }else{
      corrlist.back().puLevel_ = jet->puLevel();
      corrlist.back().puLevelUIC_ = 0;
    }
    // Store information about the "core" PT of the jet (central region)
    corrlist.back().associatedRegionEt_ = jet->associatedRegionEt();

  }

  corrlist.sort();
  corrlist.reverse();

  return corrlist;

}

void UCT2015Producer::makeEGTaus() {
  rlxTauList.clear();
  isoTauList.clear();
  rlxEGList.clear();
  isoEGList.clear();
  for(L1CaloEmCollection::const_iterator egtCand =
	newEMCands->begin();
      egtCand != newEMCands->end(); egtCand++){
    double et = egPhysicalEt(*egtCand);
    if(et > egtSeed) {
      for(L1CaloRegionCollection::const_iterator region = newRegions->begin();
	  region != newRegions->end(); region++) {
	if(egtCand->regionId().iphi() == region->gctPhi() &&
	   egtCand->regionId().ieta() == region->gctEta())
	  {
            double regionEt = regionPhysicalEt(*region);
            // Debugging
            if (false && egtCand->rank() > regionEt) {
              std::cout << "Mismatch!" << std::endl;
              std::cout << "egPhi = " << egtCand->regionId().iphi() << std::endl;
              std::cout << "egEta = " << egtCand->regionId().iphi() << std::endl;
              std::cout << "egRank = " << egtCand->rank() << std::endl;
              std::cout << "regionEt = " << regionEt << std::endl;
              std::cout << "ratio = " << egtCand->rank()*1./regionEt << std::endl;
            }
	    // A 2x1 and 1x2 cluster above egtSeed is always in tau list
	    rlxTauList.push_back(L1GObject(et, egtCand->regionId().ieta(), egtCand->regionId().iphi(), "Tau"));

            // Embed the isolation information in the L1GObject for later
            // tuning
            rlxTauList.back().associatedJetPt_ = -3; // we haven't found the jet yet.
            rlxTauList.back().puLevel_= puLevel;
	    rlxTauList.back().puLevelUIC_ = puLevelUIC;
            rlxTauList.back().associatedRegionEt_= regionEt;
            rlxTauList.back().ellIsolation_= egtCand->isolated();

	    // Note tauVeto now refers to emActivity pattern veto; Good patterns are from EG candidates
            // EKF - temporarily remove selection, do it at ntuple level.
	    //if(!region->tauVeto() && !region->mip()) {
            rlxEGList.push_back(L1GObject(et, egtCand->regionId().ieta(), egtCand->regionId().iphi(), "EG"));
            rlxEGList.back().associatedJetPt_ = -3; // we haven't found the jet yet.
            rlxEGList.back().puLevel_ = puLevel;
            rlxEGList.back().associatedRegionEt_ = regionEt;
            rlxEGList.back().ellIsolation_ = egtCand->isolated();
            rlxEGList.back().tauVeto_ = region->tauVeto();
            rlxEGList.back().mipBit_ = region->mip();
	    //}

	    // Look for overlapping jet and require that isolation be passed
	    for(list<L1GObject>::iterator jet = jetList.begin(); jet != jetList.end(); jet++) {
	      if(egtCand->regionId().iphi() == jet->phiIndex() &&
		 egtCand->regionId().ieta() == jet->etaIndex()) {
                // Embed tuning parameters into the relaxed objects
                rlxTauList.back().associatedJetPt_ = jet->pt();
                if (!region->tauVeto() && !region->mip())
                  rlxEGList.back().associatedJetPt_ = jet->pt();

		double isolation = regionEt - (regionLSB_*puLevel/9.) - et;   // Core isolation (could go less than zero)
		double relativeIsolation = isolation / et;
		double jetIsolation = jet->ptValue() - regionLSB_*puLevel - et;        // Jet isolation
		double relativeJetIsolation = jetIsolation / et;
		// A 2x1 and 1x2 cluster above egtSeed passing relative isolation will be in tau list
		if(relativeIsolation < relativeIsolationCut && relativeJetIsolation < relativeJetIsolationCut && egtCand->isolated()) {
		  isoTauList.push_back(L1GObject(et, egtCand->regionId().ieta(), egtCand->regionId().iphi(), "IsoTau"));
		  isoTauList.back().puLevel_ = puLevel;
		  isoTauList.back().puLevelUIC_ = puLevelUIC;
		  // Good patterns of EG candidate + relative isolation makes it to IsoEG
		  if(!region->tauVeto() && !region->mip()) {
		    isoEGList.push_back(L1GObject(et, egtCand->regionId().ieta(), egtCand->regionId().iphi(), "IsoEG"));
		  }
		}
		break;
	      }
	    }
	    break;
	  }
      }
    }
  }
  rlxEGList.sort();
  rlxTauList.sort();
  isoEGList.sort();
  isoTauList.sort();
  rlxEGList.reverse();
  rlxTauList.reverse();
  isoEGList.reverse();
  isoTauList.reverse();

}

//define this as a plug-in
DEFINE_FWK_MODULE(UCT2015Producer);
