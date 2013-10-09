#ifndef L1UCTGTTRANSLATOR
#define L1UCTGTTRANSLATOR

// system includes


// EDM includes
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"
#include "L1Trigger/UCT2015/interface/UCTCandidate.h"


class UCT2015GctCandsProducer : public edm::EDProducer {
 public:

  /// typedefs
  // Concrete collection of L1Gobjects (with extra tuning information)
  typedef std::vector<UCTCandidate> UCTCandidateCollection;
  typedef std::auto_ptr<UCTCandidateCollection> UCTCandidateCollectionPtr;


  /// constructor
  explicit UCT2015GctCandsProducer(const edm::ParameterSet& ps);

  /// destructor
  ~UCT2015GctCandsProducer();

 private:
  void beginJob() ;
  void produce(edm::Event& e, const edm::EventSetup& c);
  void endJob() ;

  // untracked parameters
  edm::InputTag egSourceRlx_;
  edm::InputTag egSourceIso_;

  edm::InputTag tauSourceRlx_;
  edm::InputTag tauSourceIso_;

  edm::InputTag jetSource_;
  edm::InputTag setSource_;
  edm::InputTag shtSource_;
  edm::InputTag metSource_;
  edm::InputTag mhtSource_;

  unsigned int maxEGs_;
  unsigned int maxIsoEGs_;
  unsigned int maxTaus_;
  unsigned int maxIsoTaus_;
  unsigned int maxJets_;

  // tracked parameters

};

#endif
