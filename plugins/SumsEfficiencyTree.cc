/*
 * =====================================================================================
 *
 *       Filename:  SumsEfficiencyTree.cc
 *
 *    Description:  Produces a tree of UCT event-level sums objects -
 *                  like MET, MT, MHT, etc.
 *
 *         Author:  Evan Friis, evan.friis@cern.ch with contributions from
 *                  C. Veelken (LLR)
 *
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "L1Trigger/UCT2015/src/L1GObject.h"

#include "DataFormats/METReco/interface/MET.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class SumsEfficiencyTree : public edm::EDAnalyzer {
  public:
    SumsEfficiencyTree(const edm::ParameterSet& pset);
    virtual ~SumsEfficiencyTree();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    edm::InputTag l1MHTSrc_;
    edm::InputTag l1METSrc_;
    edm::InputTag l1SHTSrc_;
    edm::InputTag l1SETSrc_;
    edm::InputTag recoMHTSrc_;
    edm::InputTag recoMETSrc_;
    edm::InputTag recoSHTSrc_;
    edm::InputTag recoSETSrc_;
    edm::InputTag pvSrc_;

    TTree* tree;

    Float_t l1MHT_;
    Float_t l1MET_;
    Float_t l1MHTPhi_;
    Float_t l1METPhi_;
    Float_t l1SHT_;
    Float_t l1SET_;

    Float_t recoMHT_;
    Float_t recoMET_;
    Float_t recoMHTPhi_;
    Float_t recoMETPhi_;
    Float_t recoSHT_;
    Float_t recoSET_;

    Int_t pvs_;
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;

    Int_t L1SingleJet36_;
    Int_t L1ETM50_;
};

SumsEfficiencyTree::SumsEfficiencyTree(const edm::ParameterSet& pset) {
  // Initialize the ntuple builder
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Ntuple", "Ntuple");

  tree->Branch("l1MHT", &l1MHT_, "l1MHT/F");
  tree->Branch("l1MET", &l1MET_, "l1MET/F");
  tree->Branch("l1MHTPhi", &l1MHTPhi_, "l1MHTPhi/F");
  tree->Branch("l1METPhi", &l1METPhi_, "l1METPhi/F");

  tree->Branch("l1SHT", &l1SHT_, "l1SHT/F");
  tree->Branch("l1SET", &l1SET_, "l1SET/F");

  tree->Branch("recoMHT", &recoMHT_, "recoMHT/F");
  tree->Branch("recoMET", &recoMET_, "recoMET/F");
  tree->Branch("recoMHTPhi", &recoMHTPhi_, "recoMHTPhi/F");
  tree->Branch("recoMETPhi", &recoMETPhi_, "recoMETPhi/F");

  tree->Branch("recoSHT", &recoSHT_, "recoSHT/F");
  tree->Branch("recoSET", &recoSET_, "recoSET/F");

  tree->Branch("pvs", &pvs_, "pvs/i");
  tree->Branch("run", &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt", &event_, "evt/l");

  // Include the actual L1 bits, to cross check the L1extra particles.
  tree->Branch("L1SingleJet36", &L1SingleJet36_, "L1SingleJet36/i");
  tree->Branch("L1ETM50", &L1ETM50_, "L1ETM50/i");

  l1MHTSrc_ = pset.getParameter<edm::InputTag>("l1MHTSrc");
  l1METSrc_ = pset.getParameter<edm::InputTag>("l1METSrc");
  l1SHTSrc_ = pset.getParameter<edm::InputTag>("l1SHTSrc");
  l1SETSrc_ = pset.getParameter<edm::InputTag>("l1SETSrc");

  recoMHTSrc_ = pset.getParameter<edm::InputTag>("recoMHTSrc");
  recoMETSrc_ = pset.getParameter<edm::InputTag>("recoMETSrc");
  recoSHTSrc_ = pset.getParameter<edm::InputTag>("recoSHTSrc");
  recoSETSrc_ = pset.getParameter<edm::InputTag>("recoSETSrc");

  pvSrc_ = pset.exists("pvSrc") ? pset.getParameter<edm::InputTag>("pvSrc") : edm::InputTag("offlinePrimaryVertices");
}

SumsEfficiencyTree::~SumsEfficiencyTree() {}

namespace {

void getValue(const edm::Event& evt, const edm::InputTag& tag, Float_t& et, Float_t& phi) {
  edm::Handle<edm::View<reco::Candidate> > handle;
  evt.getByLabel(tag, handle);
  et = handle->at(0).et();
  phi = handle->at(0).phi();
}


// Taken from C. Veelken
void checkL1bit(const DecisionWord& l1GtDecision, const L1GtTriggerMenu& l1GtTriggerMenu,
    const std::string& l1bitName, bool& l1bit_found, bool& l1bit_passed)
{
  l1bit_found = false;
  l1bit_passed = false;

  const AlgorithmMap& l1GtAlgorithms = l1GtTriggerMenu.gtAlgorithmMap();
  for ( AlgorithmMap::const_iterator l1GtAlgorithm = l1GtAlgorithms.begin();
      l1GtAlgorithm != l1GtAlgorithms.end(); ++l1GtAlgorithm ) {
    std::string l1BitName_idx = l1GtAlgorithm->second.algoName();
    int idx = l1GtAlgorithm->second.algoBitNumber();
    if ( l1BitName_idx == l1bitName ) {
      l1bit_found = true;
      l1bit_passed = l1GtDecision[idx];
      break;
    }
  }
}

void printL1bits(const DecisionWord& l1GtDecision, const L1GtTriggerMenu& l1GtTriggerMenu)
{
  std::cerr << "Existing L1 bits:" << std::endl;
  const AlgorithmMap& l1GtAlgorithms = l1GtTriggerMenu.gtAlgorithmMap();
  for ( AlgorithmMap::const_iterator l1GtAlgorithm = l1GtAlgorithms.begin();
      l1GtAlgorithm != l1GtAlgorithms.end(); ++l1GtAlgorithm ) {
    std::string l1BitName = l1GtAlgorithm->second.algoName();
    int index = l1GtAlgorithm->second.algoBitNumber();
    std::string triggerDecision = ( l1GtDecision[index] ) ? "passed" : "failed";
    std::cout << " L1 bit = " << l1BitName << ": " << triggerDecision << std::endl;
  }
}

}

void SumsEfficiencyTree::analyze(const edm::Event& evt, const edm::EventSetup& es) {

  // Setup meta info
  run_ = evt.id().run();
  lumi_ = evt.id().luminosityBlock();
  event_ = evt.id().event();

  // Get PV collection
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(pvSrc_, vertices);
  pvs_ = vertices->size();

  getValue(evt, recoMHTSrc_, recoMHT_, recoMHTPhi_);
  getValue(evt, recoMETSrc_, recoMET_, recoMETPhi_);

  Float_t dummyPhi;

  getValue(evt, recoSHTSrc_, recoSHT_, dummyPhi);
  getValue(evt, recoSETSrc_, recoSET_, dummyPhi);

  getValue(evt, l1MHTSrc_, l1MHT_, l1MHTPhi_);
  getValue(evt, l1METSrc_, l1MET_, l1METPhi_);
  getValue(evt, l1SHTSrc_, l1SHT_, dummyPhi);
  getValue(evt, l1SETSrc_, l1SET_, dummyPhi);

  // Fill the L1 bit branches
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  evt.getByLabel(edm::InputTag("gtDigis", "", "RECO"), l1GtReadoutRecord);

  edm::ESHandle<L1GtTriggerMenu> l1GtTriggerMenu;
  es.get<L1GtTriggerMenuRcd>().get(l1GtTriggerMenu);

  DecisionWord l1GtDecision = l1GtReadoutRecord->decisionWord();

  bool monitorL1Bit_found = false;
  bool monitorL1Bit_passed = false;
  checkL1bit(l1GtDecision, *l1GtTriggerMenu, "L1_SingleJet36", monitorL1Bit_found, monitorL1Bit_passed);
  if (false && !monitorL1Bit_found) {
    std::cerr << "Couldn't find the L1_SingleJet36 bit!" << std::endl;
    printL1bits(l1GtDecision, *l1GtTriggerMenu);
  }
  if (!monitorL1Bit_found)
    L1SingleJet36_ = -1;
  else
    L1SingleJet36_ = monitorL1Bit_passed;

  checkL1bit(l1GtDecision, *l1GtTriggerMenu, "L1_ETM50", monitorL1Bit_found, monitorL1Bit_passed);
  if (false && !monitorL1Bit_found)
    std::cerr << "Couldn't find the L1_ETM50 bit!" << std::endl;
  if (!monitorL1Bit_found)
    L1ETM50_ = -1;
  else
    L1ETM50_ = monitorL1Bit_passed;

  tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SumsEfficiencyTree);
