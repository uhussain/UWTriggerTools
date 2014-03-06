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
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloMipQuietRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"

#include "L1Trigger/UCT2015/interface/UCTCandidate.h"
#include "L1Trigger/UCT2015/interface/helpers.h"
#include "L1Trigger/UCT2015/interface/pileupmult_corrections.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

class Laura : public edm::EDProducer {
	public:

		// Concrete collection of output objects (with extra tuning information)
		typedef vector<UCTCandidate> UCTCandidateCollection;
		typedef std::auto_ptr<UCTCandidateCollection> UCTCandidateCollectionPtr;

		explicit Laura(const edm::ParameterSet&);

	private:
		virtual void produce(edm::Event&, const edm::EventSetup&);

		double regionPhysicalEt(const L1CaloRegion& cand) const {
			return regionLSB_*cand.et();
		}

		// Helper methods

		// ----------member data ---------------------------

		bool debug;

		unsigned int puMult;
		bool puMultCorrect;

		double regionLSB_;

		L1CaloRegionCollection LauraRegionList;

};


Laura::Laura(const edm::ParameterSet& iConfig) :
	puMultCorrect(iConfig.getParameter<bool>("puMultCorrect")),
	regionLSB_(iConfig.getParameter<double>("regionLSB"))
{
	produces<L1CaloRegionCollection>("LauraRegions");
}


	void
Laura::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if(!puMultCorrect) return;
	std::auto_ptr<L1CaloRegionCollection> LauraRegions(new L1CaloRegionCollection);

	Handle<L1CaloRegionCollection> notCorrectedRegions;

	iEvent.getByLabel("uctDigis", notCorrectedRegions);

	//-------- does something with the notCorrectedRegions
	puMult = 0;

	for(L1CaloRegionCollection::const_iterator notCorrectedRegion =
			notCorrectedRegions->begin();
			notCorrectedRegion != notCorrectedRegions->end(); notCorrectedRegion++){
		double regionET =  regionPhysicalEt(*notCorrectedRegion);
		// cout << "regionET: " << regionET <<endl; 
		if (regionET > 0) {puMult++;}
	}


	LauraRegionList.clear();
	for(L1CaloRegionCollection::const_iterator notCorrectedRegion =
			notCorrectedRegions->begin();
			notCorrectedRegion != notCorrectedRegions->end(); notCorrectedRegion++){
		double regionET =  regionPhysicalEt(*notCorrectedRegion);
		unsigned int regionEta = notCorrectedRegion->gctEta();
		//the divide by regionLSB to get back to gct Digis

		double regionEtCorr = (pumcorr(regionET, regionEta ,puMult))/regionLSB_;	

		if(regionEta<18 && regionEta>3) //if !hf
		{		

			bool overflow=notCorrectedRegion->overFlow();
			bool tau=notCorrectedRegion->tauVeto();
			bool mip=notCorrectedRegion->mip();
			bool quiet= notCorrectedRegion->quiet();
			unsigned crate=notCorrectedRegion->rctCrate();
			unsigned card=notCorrectedRegion->rctCard();
			unsigned rgn=notCorrectedRegion->rctRegionIndex();
			LauraRegionList.push_back(L1CaloRegion(regionEtCorr, overflow, tau,mip,quiet,crate,card,rgn));
		}
		else //if hf
		{
			bool fineGrain= notCorrectedRegion->fineGrain();
			unsigned crate= notCorrectedRegion->rctCrate();
			unsigned hfRgn=notCorrectedRegion->rctRegionIndex();
			LauraRegionList.push_back(L1CaloRegion(regionEtCorr,fineGrain,crate, hfRgn));
		}
	}
	for(L1CaloRegionCollection::const_iterator LaurasNewRegion = LauraRegionList.begin();
			LaurasNewRegion != LauraRegionList.end(); ++LaurasNewRegion) {
		LauraRegions->push_back(*LaurasNewRegion);
	}


	iEvent.put(LauraRegions, "LauraRegions");
}
DEFINE_FWK_MODULE(Laura);
