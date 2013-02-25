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
// $Id: UCT2015Producer.cc,v 1.22 2013/02/21 13:07:51 friis Exp $
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
#include "L1Trigger/UCT2015/interface/helpers.h"
#include "L1Trigger/UCT2015/interface/jetcorrections.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;


class UCT2015Producer : public edm::EDProducer {
public:

  static const unsigned N_JET_PHI;
  static const unsigned N_JET_ETA;

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

  double egPhysicalEt(const L1CaloEmCand& cand) const {
    return egLSB_*cand.rank();
  }

  double regionPhysicalEt(const L1CaloRegion& cand) const {
    return regionLSB_*cand.et();
  }

  // Find information about observables in the annulus.  We define the annulus
  // as all regions around the central region, with the exception of the second
  // highest in ET, as this could be sharing the 2x1.
  // MIPS in annulus refers to number of regions in the annulus which have
  // their MIP bit set.
  // egFlags is the number where (!tauVeto && !mip)
  void findAnnulusInfo(int ieta, int iphi,
      const L1CaloRegionCollection& regions,
      double* associatedSecondRegionEt,
      unsigned int* mipsInAnnulus,
      unsigned int* egFlagsInAnnulus) const;

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

  L1GObjectCollectionPtr unpackedJets(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedRlxTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrJets(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrRlxTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedCorrIsoTaus(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedRlxEGs(new L1GObjectCollection);
  L1GObjectCollectionPtr unpackedIsoEGs(new L1GObjectCollection);

  //uncorrected Jet and Tau collections
  for(list<L1GObject>::iterator jet = jetList.begin();
      jet != jetList.end();
      jet++) {
    unpackedJets->push_back(*jet);
  }
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
  //corrected Jet and Tau collections
  corrJetList=correctJets(jetList,"Jets");
  corrRlxTauList=correctJets(rlxTauList,"Tau");
  corrIsoTauList=correctJets(isoTauList,"IsoTau");

  for(list<L1GObject>::iterator jet = corrJetList.begin();
      jet != corrJetList.end();
      jet++) {
    unpackedCorrJets->push_back(*jet);
  }
  for(list<L1GObject>::iterator rlxTau = corrRlxTauList.begin();
      rlxTau != corrRlxTauList.end();
      rlxTau++) {
    unpackedCorrRlxTaus->push_back(*rlxTau);
  }
  for(list<L1GObject>::iterator isoTau = corrIsoTauList.begin();
      isoTau != corrIsoTauList.end();
      isoTau++) {
    unpackedCorrIsoTaus->push_back(*isoTau);
  }

  // egamma collections
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

  iEvent.put(collectionize(puLevel), "PULevelUnpacked");
  iEvent.put(collectionize(puLevelUIC), "PULevelUICUnpacked");
  iEvent.put(collectionize(METObject), "METUnpacked");
  iEvent.put(collectionize(MHTObject), "MHTUnpacked");
  iEvent.put(collectionize(SETObject), "SETUnpacked");
  iEvent.put(collectionize(SHTObject), "SHTUnpacked");

  iEvent.put(unpackedJets, "JetUnpacked");
  iEvent.put(unpackedRlxTaus, "RelaxedTauUnpacked");
  iEvent.put(unpackedIsoTaus, "IsolatedTauUnpacked");
  iEvent.put(unpackedCorrJets, "CorrJetUnpacked");
  iEvent.put(unpackedCorrRlxTaus, "CorrRelaxedTauUnpacked");
  iEvent.put(unpackedCorrIsoTaus, "CorrIsolatedTauUnpacked");
  iEvent.put(unpackedRlxEGs, "RelaxedEGUnpacked");
  iEvent.put(unpackedIsoEGs, "IsolatedEGUnpacked");

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
      Rarea += getRegionArea(newRegion->gctEta());
    }
  }
  // Add a factor of 9, so it corresponds to a jet.  Reduces roundoff error.
  puLevel *= 9;
  puLevel = puLevel / puCount;
  r_puLevelUIC = r_puLevelUIC / Rarea;
  puLevelUIC=0;
  if (r_puLevelUIC > 0.) puLevelUIC = floor (r_puLevelUIC + 0.5);

}

void UCT2015Producer::makeSums()
{
  sumET = 0;
  sumEx = 0;
  sumEy = 0;
  sumHT = 0;
  sumHx = 0;
  sumHy = 0;

  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++){
    // Remove forward stuff
    if (newRegion->gctEta() < minGctEtaForSums || newRegion->gctEta() > maxGctEtaForSums) {
      continue;
    }
    //unsigned int regionET = newRegion->et() - puLevel;
    double regionET = std::max(regionPhysicalEt(*newRegion) - puLevel*regionLSB_/9., 0.);
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
}

int deltaGctPhi(const L1CaloRegion& r1, const L1CaloRegion& r2) {
  return deltaPhiWrapAtN(18, r1.gctPhi(), r2.gctPhi());
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
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborS_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborE_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborNE_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborSW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborNW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborSE_et = neighborET;
          nNeighbors++;
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

// Given a region at iphi/ieta, find the highest region in the surrounding
// regions.
void UCT2015Producer::findAnnulusInfo(int ieta, int iphi,
    const L1CaloRegionCollection& regions,
    double* associatedSecondRegionEt,
    unsigned int* mipsInAnnulus,
    unsigned int* egFlagsInAnnulus) const {

  unsigned int neighborsFound = 0;
  unsigned int mipsCount = 0;
  unsigned int egFlagCount = 0;
  double highestNeighborEt = 0;
  // We don't want to count the contribution of the highest neighbor, this allows
  // us to subtract off the highest neighbor at the end, so we only loop once.
  bool highestNeighborHasMip = false;
  bool highestNeighborHasEGFlag = false;

  for(L1CaloRegionCollection::const_iterator region = regions.begin();
      region != regions.end(); region++) {
    int regionPhi = region->gctPhi();
    int regionEta = region->gctEta();
    unsigned int deltaPhi = std::abs(deltaPhiWrapAtN(18, iphi, regionPhi));
    unsigned int deltaEta = std::abs(ieta - regionEta);
    if ((deltaPhi + deltaEta) > 0 && deltaPhi < 2 && deltaEta < 2) {
      double regionET = regionPhysicalEt(*region);
      if (regionET > highestNeighborEt) {
        highestNeighborEt = regionET;
        // Keep track of what flags the highest neighbor has
        highestNeighborHasMip = region->mip();
        highestNeighborHasEGFlag = !region->mip() && !region->tauVeto();
      }

      // count how many neighbors pass the flags.
      if (region->mip()) {
        mipsCount++;
      }
      if (!region->mip() && !region->tauVeto()) {
        egFlagCount++;
      }

      // If we already found all 8 neighbors, we don't need to keep looping
      // over the regions.
      neighborsFound++;
      if (neighborsFound == 8) {
        break;
      }
    }
  }
  // check if we need to remove the highest neighbor from the flag count.
  if (highestNeighborHasMip)
    mipsCount--;
  if (highestNeighborHasEGFlag)
    egFlagCount--;

  // set output
  *associatedSecondRegionEt = highestNeighborEt;
  *mipsInAnnulus = mipsCount;
  *egFlagsInAnnulus = egFlagCount;
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

            // Find the highest region in the 3x3 annulus around the center
            // region.
            double associatedSecondRegionEt = 0;
            unsigned int mipsInAnnulus = 0;
            unsigned int egFlagsInAnnulus = 0;
            findAnnulusInfo(
                egtCand->regionId().ieta(), egtCand->regionId().iphi(),
                *newRegions,
                &associatedSecondRegionEt, &mipsInAnnulus, &egFlagsInAnnulus);

            // Embed the isolation information in the L1GObject for later
            // tuning
            rlxTauList.back().associatedJetPt_ = -3; // we haven't found the jet yet.
            rlxTauList.back().puLevel_= puLevel;
	    rlxTauList.back().puLevelUIC_ = puLevelUIC;
            rlxTauList.back().associatedRegionEt_= regionEt;
            rlxTauList.back().associatedSecondRegionEt_= associatedSecondRegionEt;
            rlxTauList.back().mipsInAnnulus_= mipsInAnnulus;
            rlxTauList.back().egFlagsInAnnulus_= egFlagsInAnnulus;
            rlxTauList.back().ellIsolation_= egtCand->isolated();

	    // Note tauVeto now refers to emActivity pattern veto; Good patterns are from EG candidates
            // EKF - temporarily remove selection, do it at ntuple level.
	    //if(!region->tauVeto() && !region->mip()) {
            rlxEGList.push_back(L1GObject(et, egtCand->regionId().ieta(), egtCand->regionId().iphi(), "EG"));
            rlxEGList.back().associatedJetPt_ = -3; // we haven't found the jet yet.
            rlxEGList.back().puLevel_ = puLevel;
            rlxEGList.back().associatedRegionEt_ = regionEt;
            rlxEGList.back().associatedSecondRegionEt_= associatedSecondRegionEt;
            rlxEGList.back().mipsInAnnulus_= mipsInAnnulus;
            rlxEGList.back().egFlagsInAnnulus_= egFlagsInAnnulus;
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
