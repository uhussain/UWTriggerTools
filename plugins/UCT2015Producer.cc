//
// Package:    UCT2015Producer
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
  typedef vector<UCTCandidate> UCTCandidateCollection;
  typedef std::auto_ptr<UCTCandidateCollection> UCTCandidateCollectionPtr;

  explicit UCT2015Producer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

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
      unsigned int* egFlagsInAnnulus,
      unsigned int* mipInSecondRegion) const;

  // Helper methods

  void puSubtraction();

  void makeSums();
  void makeJets();
  void makeEGTaus();

  list<UCTCandidate> correctJets(const list<UCTCandidate>&);

  // ----------member data ---------------------------

  bool puCorrect;
  bool useUICrho; // which PU denstity to use for energy correction determination
  bool useHI; // do HI-style background subtraction

  unsigned int puETMax;
  unsigned int puLevel;
  //double puLevelUIC; // puLevel divided by puCount*Area, not multiply by 9.0
  unsigned int  puLevelUIC; // puLevel divided by puCount*Area, not multiply by 9.0
  vector<int> puLevelHI;

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

  UCTCandidate METObject;
  UCTCandidate MHTObject;
  UCTCandidate SETObject;
  UCTCandidate SHTObject;

  unsigned int jetSeed;
  list<UCTCandidate> jetList, corrJetList;

  unsigned int egtSeed;
  double relativeIsolationCut;
  double relativeJetIsolationCut;
  list<UCTCandidate> rlxTauList, corrRlxTauList;
  list<UCTCandidate> rlxEGList;
  list<UCTCandidate> isoTauList, corrIsoTauList;
  list<UCTCandidate> isoEGList;

  Handle<L1CaloRegionCollection> newRegions;
  Handle<L1CaloEmCollection> newEMCands;

  vector<double> sinPhi;
  vector<double> cosPhi;

  double egLSB_;
  double regionLSB_;

};

unsigned const UCT2015Producer::N_JET_PHI = L1CaloRegionDetId::N_PHI * 4;
unsigned const UCT2015Producer::N_JET_ETA = L1CaloRegionDetId::N_ETA * 4;

//
// constructors and destructor
//
UCT2015Producer::UCT2015Producer(const edm::ParameterSet& iConfig) :
  puCorrect(iConfig.getParameter<bool>("puCorrect")),
  useUICrho(iConfig.getParameter<bool>("useUICrho")),
  useHI(iConfig.getParameter<bool>("useHI")),
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
  puLevelHI.resize(L1CaloRegionDetId::N_ETA);
  for(unsigned i = 0; i < L1CaloRegionDetId::N_ETA; ++i)
    puLevelHI[i] = 0;

  // Also declare we produce unpacked collections (which have more info)
  produces<UCTCandidateCollection>( "JetUnpacked" ) ;
  produces<UCTCandidateCollection>( "CorrJetUnpacked" ) ;
  produces<UCTCandidateCollection>( "RelaxedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "RelaxedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "CorrRelaxedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "CorrIsolatedTauUnpacked" ) ;

  produces<UCTCandidateCollection>( "PULevelUnpacked" ) ;
  produces<UCTCandidateCollection>( "PULevelUICUnpacked" ) ;
  produces<UCTCandidateCollection>( "METUnpacked" ) ;
  produces<UCTCandidateCollection>( "MHTUnpacked" ) ;
  produces<UCTCandidateCollection>( "SETUnpacked" ) ;
  produces<UCTCandidateCollection>( "SHTUnpacked" ) ;

  //now do what ever initialization is needed
  for(unsigned int i = 0; i < L1CaloRegionDetId::N_PHI; i++) {
    sinPhi.push_back(sin(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
    cosPhi.push_back(cos(2. * 3.1415927 * i * 1.0 / L1CaloRegionDetId::N_PHI));
  }
}


// For the single objects, like MET/MHT, etc, convert them into a
// std::auto_ptr<UCTCandidateCollection> suitable for putting into the edm::Event
// The "collection" contains only 1 object.
UCT2015Producer::UCTCandidateCollectionPtr collectionize(const UCTCandidate& obj) {
  return UCT2015Producer::UCTCandidateCollectionPtr(
      new UCT2015Producer::UCTCandidateCollection(1, obj));
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

  UCTCandidateCollectionPtr unpackedJets(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedRlxTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedCorrJets(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedCorrRlxTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedCorrIsoTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedRlxEGs(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoEGs(new UCTCandidateCollection);

  //uncorrected Jet and Tau collections
  for(list<UCTCandidate>::iterator jet = jetList.begin();
      jet != jetList.end();
      jet++) {
    unpackedJets->push_back(*jet);
  }
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
  //corrected Jet and Tau collections
  corrJetList = correctJets(jetList);
  corrRlxTauList = correctJets(rlxTauList);
  corrIsoTauList = correctJets(isoTauList);

  for(list<UCTCandidate>::iterator jet = corrJetList.begin();
      jet != corrJetList.end();
      jet++) {
    unpackedCorrJets->push_back(*jet);
  }
  for(list<UCTCandidate>::iterator rlxTau = corrRlxTauList.begin();
      rlxTau != corrRlxTauList.end();
      rlxTau++) {
    unpackedCorrRlxTaus->push_back(*rlxTau);
  }
  for(list<UCTCandidate>::iterator isoTau = corrIsoTauList.begin();
      isoTau != corrIsoTauList.end();
      isoTau++) {
    unpackedCorrIsoTaus->push_back(*isoTau);
  }

  // egamma collections
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

  // Just store these as cands to make life easier.
  UCTCandidate puLevelAsCand(puLevel, 0, 0);
  UCTCandidate puLevelUICAsCand(puLevel, 0, 0);

  iEvent.put(collectionize(puLevelAsCand), "PULevelUnpacked");
  iEvent.put(collectionize(puLevelUICAsCand), "PULevelUICUnpacked");
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

// NB PU is not in the physical scale!!  Needs to be multiplied by regionLSB
void UCT2015Producer::puSubtraction()
{
  puLevel = 0;
  puLevelUIC = 0;
  double r_puLevelUIC=0.0;
  double r_puLevelHI[L1CaloRegionDetId::N_ETA];
  int etaCount[L1CaloRegionDetId::N_ETA];
  for(unsigned i = 0; i < L1CaloRegionDetId::N_ETA; ++i)
  {
    puLevelHI[i] = 0;
    r_puLevelHI[i] = 0.0;
    etaCount[i] = 0;
  }

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
    r_puLevelHI[newRegion->gctEta()] += newRegion->et();
    etaCount[newRegion->gctEta()]++;
  }
  // Add a factor of 9, so it corresponds to a jet.  Reduces roundoff error.
  puLevel *= 9;
  if(puCount != 0) puLevel = puLevel / puCount;
  r_puLevelUIC = r_puLevelUIC / Rarea;
  puLevelUIC=0;
  if (r_puLevelUIC > 0.) puLevelUIC = floor (r_puLevelUIC + 0.5);

  for(unsigned i = 0; i < L1CaloRegionDetId::N_ETA; ++i)
  {
    puLevelHI[i] = floor(r_puLevelHI[i]/etaCount[i] + 0.5);
  }
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

  double physicalPhi = atan2(sumEy, sumEx) + 3.1415927;
  unsigned int iPhi = L1CaloRegionDetId::N_PHI * physicalPhi / (2 * 3.1415927);
  METObject = UCTCandidate(MET, 0, physicalPhi);
  METObject.setInt("rgnPhi", iPhi);

  double physicalPhiHT = atan2(sumHy, sumHx) + 3.1415927;
  iPhi = L1CaloRegionDetId::N_PHI * (physicalPhiHT) / (2 * 3.1415927);
  MHTObject = UCTCandidate(MHT, 0, physicalPhiHT);
  MHTObject.setInt("rgnPhi", iPhi);

  SETObject = UCTCandidate(sumET, 0, 0);
  SHTObject = UCTCandidate(sumHT, 0, 0);
}

void UCT2015Producer::makeJets() {
  jetList.clear();
  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++) {
    double regionET = regionPhysicalEt(*newRegion);
    if(puCorrect && useHI)
      regionET = std::max(0.,regionET -
			  (puLevelHI[newRegion->gctEta()]*regionLSB_));
    if((regionET > jetSeed) || (puCorrect && useHI)) {
      double neighborN_et = 0;
      double neighborS_et = 0;
      double neighborE_et = 0;
      double neighborW_et = 0;
      double neighborNE_et = 0;
      double neighborSW_et = 0;
      double neighborNW_et = 0;
      double neighborSE_et = 0;
      unsigned int nNeighbors = 0;
      for(L1CaloRegionCollection::const_iterator neighbor = newRegions->begin();
	  neighbor != newRegions->end(); neighbor++) {
	double neighborET = regionPhysicalEt(*neighbor);
	if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	   (newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborN_et = neighborET;
	  if(puCorrect && useHI)
	    neighborN_et = std::max(0.,neighborET -
				    (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborS_et = neighborET;
	  if(puCorrect && useHI)
	    neighborS_et = std::max(0.,neighborET -
				    (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborE_et = neighborET;
	  if(puCorrect && useHI)
	    neighborE_et = std::max(0.,neighborET -
				    (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborW_et = neighborET;
	  if(puCorrect && useHI)
	    neighborW_et = std::max(0.,neighborET -
				    (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborNE_et = neighborET;
	  if(puCorrect && useHI)
	    neighborNE_et = std::max(0.,neighborET -
				     (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborSW_et = neighborET;
	  if(puCorrect && useHI)
	    neighborSW_et = std::max(0.,neighborET -
				     (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborNW_et = neighborET;
	  if(puCorrect && useHI)
	    neighborNW_et = std::max(0.,neighborET -
				     (puLevelHI[neighbor->gctEta()]*regionLSB_));
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborSE_et = neighborET;
	  if(puCorrect && useHI)
	    neighborSE_et = std::max(0.,neighborET -
				     (puLevelHI[neighbor->gctEta()]*regionLSB_));
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
        UCTCandidate theJet(jetET, convertRegionEta(jetEta), convertRegionPhi(jetPhi));
        theJet.setInt("rgnEta", jetEta);
        theJet.setInt("rgnPhi", jetPhi);
        theJet.setInt("rctEta",  newRegion->rctEta());
        theJet.setInt("rctPhi", newRegion->rctPhi());

        // Embed the puLevel information in the jet object for later tuning
        theJet.setFloat("puLevel", puLevel);
        theJet.setFloat("puLevelUIC", puLevelUIC);
        // Store information about the "core" PT of the jet (central region)
        theJet.setFloat("associatedRegionEt", regionET);
        jetList.push_back(theJet);
      }
    }
  }
  jetList.sort();
  jetList.reverse();
}

list<UCTCandidate>
UCT2015Producer::correctJets(const list<UCTCandidate>& jets) {
  // jet corrections only valid if PU density has been calculated
  list<UCTCandidate> corrlist;
  if (!puCorrect) return corrlist;

  corrlist.clear();

  for(list<UCTCandidate>::const_iterator jet = jets.begin(); jet != jets.end(); jet++) {

    const double jetET=jet->pt();
    double jpt = 0;
    unsigned int corjetET = 0;

    //apply Michael's jet correction function
    if (useUICrho){
      jpt = jetcorrUIC(jetET, jet->getInt("rgnEta"), jet->getFloat("puLevelUIC"));
    }else if (useHI) {
      jpt = jetET;
    }else{
      jpt = jetcorr(jetET, jet->getInt("rgnEta"), jet->getFloat("puLevel"));
    }
    if (jpt>0)
      corjetET = floor(jpt + 0.5);

    UCTCandidate newJet = *jet;
    newJet.setP4(reco::LeafCandidate::PolarLorentzVector(
          corjetET, jet->eta(), jet->phi(), jet->mass()));
    newJet.setFloat("uncorrectedPt", jetET);

    corrlist.push_back(newJet);
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
    unsigned int* egFlagsInAnnulus,
    unsigned int* mipInSecondRegion) const {

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
  *mipInSecondRegion = highestNeighborHasMip;
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

            // Find the highest region in the 3x3 annulus around the center
            // region.
            double associatedSecondRegionEt = 0;
            unsigned int mipsInAnnulus = 0;
            unsigned int egFlagsInAnnulus = 0;
            unsigned int mipInSecondRegion = 0;
            findAnnulusInfo(
                egtCand->regionId().ieta(), egtCand->regionId().iphi(),
                *newRegions,
                &associatedSecondRegionEt, &mipsInAnnulus, &egFlagsInAnnulus,
                &mipInSecondRegion);

            UCTCandidate egtauCand(
                et,
                convertRegionEta(egtCand->regionId().ieta()),
                convertRegionPhi(egtCand->regionId().iphi()));

            // Add extra information to the candidate
            egtauCand.setInt("rgnEta", egtCand->regionId().ieta());
            egtauCand.setInt("rgnPhi", egtCand->regionId().iphi());
            egtauCand.setInt("rctEta", egtCand->regionId().rctEta());
            egtauCand.setInt("rctPhi", egtCand->regionId().rctPhi());
            egtauCand.setInt("rank", egtCand->rank());
            egtauCand.setFloat("test", et);
            egtauCand.setFloat("associatedJetPt", -3);
            egtauCand.setFloat("associatedRegionEt", regionEt);
            egtauCand.setFloat("associatedSecondRegionEt", associatedSecondRegionEt);
            egtauCand.setInt("associatedSecondRegionMIP", mipInSecondRegion);
            egtauCand.setFloat("puLevel", puLevel);
            egtauCand.setFloat("puLevelUIC", puLevelUIC);
            egtauCand.setInt("ellIsolation", egtCand->isolated());
            egtauCand.setInt("tauVeto", region->tauVeto());
            egtauCand.setInt("mipBit", region->mip());

	    // A 2x1 and 1x2 cluster above egtSeed is always in tau list
            rlxTauList.push_back(egtauCand);

	    // Note tauVeto now refers to emActivity pattern veto;
            // Good patterns are from EG candidates
            // EKF - temporarily remove selection, do it at ntuple level.

	    //if(!region->tauVeto() && !region->mip()) {
            rlxEGList.push_back(egtauCand);
	    //}

	    // Look for overlapping jet and require that isolation be passed
	    for(list<UCTCandidate>::iterator jet = jetList.begin(); jet != jetList.end(); jet++) {

	      if((int)egtCand->regionId().iphi() == jet->getInt("rgnPhi") &&
		 (int)egtCand->regionId().ieta() == jet->getInt("rgnEta")) {
                // Embed tuning parameters into the relaxed objects
                rlxTauList.back().setFloat("associatedJetPt", jet->pt());
                // EG ID disabled - EKF
                //if (!region->tauVeto() && !region->mip())
                  rlxEGList.back().setFloat("associatedJetPt", jet->pt());

		double isolation = regionEt - (regionLSB_*puLevel/9.) - et;   // Core isolation (could go less than zero)
		double relativeIsolation = isolation / et;
		double jetIsolation = jet->pt() - regionLSB_*puLevel - et;        // Jet isolation
		double relativeJetIsolation = jetIsolation / et;
		// A 2x1 and 1x2 cluster above egtSeed passing relative isolation will be in tau list
		if(relativeIsolation < relativeIsolationCut && relativeJetIsolation < relativeJetIsolationCut && egtCand->isolated()) {
                  isoTauList.push_back(rlxTauList.back());
		  // Good patterns of EG candidate + relative isolation makes it to IsoEG
		  if(!region->tauVeto() && !region->mip()) {
		    isoEGList.push_back(rlxEGList.back());
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
