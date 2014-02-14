/*
 * =====================================================================================
 *
 *       Filename:  puMultipTreeRegionEt.cc
 *
 *    Description:  Produce a tree for computing rates.
 *
 *         Author:  Laura Dodd, laura.dodd@cern.ch Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 * =====================================================================================
 */
#include "L1Trigger/UCT2015/interface/ExpressionNtuple.h"
#include "L1Trigger/UCT2015/interface/L1RecoMatch.h"


#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TTree.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "L1Trigger/UCT2015/interface/UCTCandidate.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//typedef std::vector<edm::InputTag> VInputTag;
//typedef std::vector<unsigned int> PackedUIntCollection;


using namespace std;
using namespace edm;


class puMultipTreeRegionEt : public edm::EDAnalyzer {
	public:
		explicit puMultipTreeRegionEt(const edm::ParameterSet& pset);

	private:
		virtual void analyze(const edm::Event& evt, const edm::EventSetup& es);
		double egPhysicalEt(const L1CaloEmCand& cand) const {
			return egLSB_*cand.rank();
		}

		double regionPhysicalEt(const L1CaloRegion& cand) const {
			return regionLSB_*cand.et();
		}

		TTree* tree;
		unsigned int run_;
		unsigned int lumi_;
		unsigned int puMult0_;
		unsigned int puMult1_;
		unsigned int puMult2_;
		unsigned int puMult3_;
		unsigned int puMult4_;
		unsigned int puMult5_;
		unsigned long int event_;
		InputTag scalerSrc_;
		InputTag uctDigis_;
		InputTag pvSrc_;
		float instLumi_;
		unsigned int npvs_;
		Handle<L1CaloRegionCollection> newRegions;
		Handle<LumiScalersCollection> lumiScalers;
		Handle<reco::VertexCollection> vertices;


                vector<float> regionPt_;
                vector<int> etas_;

		vector<double> sinPhi;
		vector<double> cosPhi;

		double egLSB_;
		double regionLSB_;
		// Add UCT-only variables
		bool isUCT_;
		bool useHF_;
		unsigned int minGctEtaForSums;
		unsigned int maxGctEtaForSums;
};



puMultipTreeRegionEt::puMultipTreeRegionEt(const edm::ParameterSet& pset) 
{
	// Initialize the ntuple builder
	edm::Service<TFileService> fs;
	tree = fs->make<TTree>("Ntuple", "Ntuple");
        tree->Branch("regionPt", "std::vector<float>", &regionPt_);
        tree->Branch("regionEta", "std::vector<int>", &etas_);

	tree->Branch("run", &run_, "run/i");
	tree->Branch("lumi", &lumi_, "lumi/i");
	tree->Branch("evt", &event_, "evt/l");
	tree->Branch("npvs", &npvs_, "npvs/i");
	tree->Branch("instlumi", &instLumi_, "instlumi/F");
	tree->Branch("puMult0", &puMult0_, "puMult0/i");
	tree->Branch("puMult1", &puMult1_, "puMult1/i");
	tree->Branch("puMult2", &puMult2_, "puMult2/i");
	tree->Branch("puMult3", &puMult3_, "puMult3/i");
	tree->Branch("puMult4", &puMult4_, "puMult4/i");
	tree->Branch("puMult5", &puMult5_, "puMult5/i");
	scalerSrc_ = pset.exists("scalerSrc") ? pset.getParameter<InputTag>("scalerSrc") : InputTag("scalersRawToDigi");
	// UCT variables
	isUCT_ = pset.getParameter<bool>("isUCT");
	useHF_ = pset.getParameter<bool>("useHF");
	pvSrc_ = pset.exists("pvSrc") ? pset.getParameter<InputTag>("pvSrc") : InputTag("offlinePrimaryVertices");
	regionLSB_ = pset.getParameter<double>("regionLSB");
	minGctEtaForSums = pset.getParameter<unsigned int>("minGctEtaForSums");
	maxGctEtaForSums = pset.getParameter<unsigned int>("maxGctEtaForSums");
}


void puMultipTreeRegionEt::analyze(const edm::Event& evt, const edm::EventSetup& es) {

	// Setup meta info
	run_ = evt.id().run();
	lumi_ = evt.id().luminosityBlock();
	event_ = evt.id().event();

        etas_.clear();
        regionPt_.clear();

	// Get instantaneous lumi from the scalers
	// thx to Carlo Battilana
	//Handle<LumiScalersCollection> lumiScalers;
	//Handle<L1CaloRegionCollection> newRegions;
	//edm::DetSetVector<L1CaloRegionCollection> newRegion;
	//edm::Handle<L1CaloRegionCollection>::const_iterator newRegion;

	evt.getByLabel(scalerSrc_, lumiScalers);
	evt.getByLabel("uctDigis", newRegions);
	evt.getByLabel(pvSrc_, vertices);


	instLumi_ = -1;
	npvs_ = 0;
	puMult0_ = 0;
	puMult1_ = 0;
	puMult2_ = 0;
	puMult3_ = 0;
	puMult4_ = 0;
	puMult5_ = 0;

	npvs_ = vertices->size();

	if (lumiScalers->size())
		instLumi_ = lumiScalers->begin()->instantLumi();


	for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin(); newRegion != newRegions->end(); newRegion++)
	{
		// Remove forward stuff
		if (!useHF_)
		{
			if (newRegion->gctEta() < minGctEtaForSums || newRegion->gctEta() > maxGctEtaForSums) {
				continue;
			}
		}
		double regionET =  regionPhysicalEt(*newRegion);
		unsigned int regionEta = newRegion->gctEta(); 
		// cout << "regionET: " << regionET <<endl; 
		regionPt_.push_back(regionET);
		etas_.push_back(regionEta);
		if (regionET > 0) {puMult0_++;}
		if (regionET < 1) {puMult1_++;}
		if (regionET < 2) {puMult2_++;}
		if (regionET < 3) {puMult3_++;}
		if (regionET < 4) {puMult4_++;}
		if (regionET < 5) {puMult5_++;}
	}

	tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(puMultipTreeRegionEt);
