/*
 * =====================================================================================
 *
 *       Filename:  RateTree.cc
 *
 *    Description:  Produce a tree for computing rates.
 *
 *         Author:  Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "L1Trigger/UCT2015/src/L1GObject.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class RateTree : public edm::EDAnalyzer {
  public:
    RateTree(const edm::ParameterSet& pset);
    virtual ~RateTree();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    VInputTag src_;
    TTree* tree;
    std::vector<Float_t>* pts_;
    std::vector<Float_t>* etas_;
    std::vector<Float_t>* phis_;
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;
    edm::InputTag scalerSrc_;
    Float_t instLumi_;
    // Add UCT-only variables
    bool isUCT_;
    std::vector<Float_t>* jetPt_;
    std::vector<Float_t>* regionPt_;
    std::vector<Float_t>* secondRegionPt_;
    std::vector<Int_t>* ellIso_;
    std::vector<Int_t>* pu_;
    std::vector<Int_t>* puUIC_;
    std::vector<bool>* taus_;
    std::vector<bool>* mips_;

};

RateTree::RateTree(const edm::ParameterSet& pset) {
  // Initialize the ntuple builder
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Ntuple", "Ntuple");
  pts_ = new std::vector<Float_t>();
  etas_ = new std::vector<Float_t>();
  phis_ = new std::vector<Float_t>();
  taus_ = new std::vector<bool>();
  mips_ = new std::vector<bool>();


  tree->Branch("pt", "std::vector<float>", &pts_);
  tree->Branch("eta", "std::vector<float>", &etas_);
  tree->Branch("phi", "std::vector<float>", &phis_);
  tree->Branch("run", &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt", &event_, "evt/l");
  tree->Branch("instlumi", &instLumi_, "instlumi/F");

  // UCT variables
  isUCT_ = pset.getParameter<bool>("isUCT");
  jetPt_ = new std::vector<Float_t>();
  regionPt_ = new std::vector<Float_t>();
  secondRegionPt_ = new std::vector<Float_t>();
  ellIso_ = new std::vector<Int_t>();
  pu_ = new std::vector<Int_t>();
  puUIC_ = new std::vector<Int_t>();

  if (isUCT_) {
    tree->Branch("jetPt", "std::vector<float>", &jetPt_);
    tree->Branch("regionPt", "std::vector<float>", &regionPt_);
    tree->Branch("secondRegionPt", "std::vector<float>", &secondRegionPt_);
    tree->Branch("ellIso", "std::vector<int>", &ellIso_);
    tree->Branch("pu", "std::vector<int>", &pu_);
    tree->Branch("puUIC", "std::vector<int>", &puUIC_);
    tree->Branch("tauVeto", "std::vector<bool>", &taus_);
    tree->Branch("mipBit", "std::vector<bool>", &mips_);
  }

  src_ = pset.getParameter<VInputTag>("src");
  scalerSrc_ = pset.exists("scalerSrc") ?
    pset.getParameter<edm::InputTag>("scalerSrc") : edm::InputTag("scalersRawToDigi");
}

RateTree::~RateTree() {
  delete pts_;
  delete etas_;
  delete phis_;
  delete jetPt_;
  delete regionPt_;
  delete ellIso_;
  delete pu_;
  delete puUIC_;
  delete mips_;
  delete taus_;
}


namespace {

  // Predicate to sort candidates by descending pt
  class CandPtSorter {
    public:
      bool operator()(const reco::Candidate* candA, const reco::Candidate* candB)
        const {
          return candA->pt() > candB->pt();
        }
  };

  // Turn a set of InputTags into a colleciton of candidate pointers.
  std::vector<const reco::Candidate*> getCollections(
      const edm::Event& evt, const VInputTag& collections) {
    std::vector<const reco::Candidate*> output;
    // Loop over collections
    for (size_t i = 0; i < collections.size(); ++i) {
      edm::Handle<edm::View<reco::Candidate> > handle;
      evt.getByLabel(collections[i], handle);
      // Loop over objects in current collection
      for (size_t j = 0; j < handle->size(); ++j) {
        const reco::Candidate& object = handle->at(j);
        output.push_back(&object);
      }
    }
    return output;
  }

}

void RateTree::analyze(const edm::Event& evt, const edm::EventSetup& es) {

  // Get the objects.
  std::vector<const reco::Candidate*> objects = getCollections(
      evt, src_);

  std::sort(objects.begin(), objects.end(), CandPtSorter());


  // Clear previous event's objects
  pts_->clear();
  etas_->clear();
  phis_->clear();
  jetPt_->clear();
  regionPt_->clear();
  secondRegionPt_->clear();
  ellIso_->clear();
  pu_->clear();
  puUIC_->clear();
  taus_->clear();
  mips_->clear();

  // Setup meta info
  run_ = evt.id().run();
  lumi_ = evt.id().luminosityBlock();
  event_ = evt.id().event();

  // Get instantaneous lumi from the scalers
  // thx to Carlo Battilana
  edm::Handle<LumiScalersCollection> lumiScalers;
  evt.getByLabel(scalerSrc_, lumiScalers);
  instLumi_ = -1;
  if (lumiScalers->size())
    instLumi_ = lumiScalers->begin()->instantLumi();

  for (size_t i = 0; i < objects.size(); ++i) {
    pts_->push_back(objects[i]->pt());
    etas_->push_back(objects[i]->eta());
    phis_->push_back(objects[i]->phi());
    if (isUCT_) {
      const L1GObject* l1g = dynamic_cast<const L1GObject*>(objects[i]);
      assert(l1g);
      jetPt_->push_back(l1g->associatedJetPt());
      regionPt_->push_back(l1g->associatedRegionEt());
      secondRegionPt_->push_back(l1g->associatedSecondRegionEt());
      ellIso_->push_back(l1g->ellIsolation());
      pu_->push_back(l1g->puLevel());
      puUIC_->push_back(l1g->puLevelUIC());
      mips_->push_back(l1g->mipBit());
      taus_->push_back(l1g->tauVeto());
    }
  }

  // pad everything, to work around the MaxIf bug.
  pts_->push_back(-9999);
  etas_->push_back(-9999);
  phis_->push_back(-9999);
  jetPt_->push_back(-9999);
  regionPt_->push_back(-9999);
  secondRegionPt_->push_back(-9999);
  ellIso_->push_back(-9999);
  pu_->push_back(-9999);
  puUIC_->push_back(-9999);
  mips_->push_back(-9999);
  taus_->push_back(-9999);

  tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RateTree);
