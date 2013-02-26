/*
 * =====================================================================================
 *
 *       Filename:  SumsRateTree.cc
 *
 *    Description:  Produces a tree of UCT event-level sums objects -
 *                  like MET, MT, MHT, etc.
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
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class SumsRateTree : public edm::EDAnalyzer {
  public:
    SumsRateTree(const edm::ParameterSet& pset);
    virtual ~SumsRateTree();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::InputTag l1MHTSrc_;
    edm::InputTag l1METSrc_;
    edm::InputTag l1SHTSrc_;
    edm::InputTag l1SETSrc_;

    TTree* tree;

    Float_t l1MHT_;
    Float_t l1MET_;
    Float_t l1SHT_;
    Float_t l1SET_;

    Int_t pu_;
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;

    edm::InputTag scalerSrc_;
    Float_t instLumi_;
};

SumsRateTree::SumsRateTree(const edm::ParameterSet& pset) {
  // Initialize the ntuple builder
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Ntuple", "Ntuple");

  tree->Branch("l1MHT", &l1MHT_, "l1MHT/F");
  tree->Branch("l1MET", &l1MET_, "l1MET/F");
  tree->Branch("l1SHT", &l1SHT_, "l1SHT/F");
  tree->Branch("l1SET", &l1SET_, "l1SET/F");

  tree->Branch("pu", &run_, "pu/i");
  tree->Branch("run", &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt", &event_, "evt/l");
  tree->Branch("instlumi", &instLumi_, "instlumi/F");

  l1MHTSrc_ = pset.getParameter<edm::InputTag>("l1MHTSrc");
  l1METSrc_ = pset.getParameter<edm::InputTag>("l1METSrc");
  l1SHTSrc_ = pset.getParameter<edm::InputTag>("l1SHTSrc");
  l1SETSrc_ = pset.getParameter<edm::InputTag>("l1SETSrc");

  scalerSrc_ = pset.exists("scalerSrc") ?
    pset.getParameter<edm::InputTag>("scalerSrc") : edm::InputTag("scalersRawToDigi");
}

SumsRateTree::~SumsRateTree() {}

namespace {

double getValue(const edm::Event& evt, const edm::InputTag& tag) {
  edm::Handle<edm::View<reco::Candidate> > handle;
  evt.getByLabel(tag, handle);
  return handle->at(0).et();
}

}

void SumsRateTree::analyze(const edm::Event& evt, const edm::EventSetup& es) {

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

  l1MHT_ = getValue(evt, l1MHTSrc_);
  l1MET_ = getValue(evt, l1METSrc_);
  l1SHT_ = getValue(evt, l1SHTSrc_);
  l1SET_ = getValue(evt, l1SETSrc_);

  tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SumsRateTree);
