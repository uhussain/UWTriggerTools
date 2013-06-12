#include "L1Trigger/UCT2015/interface/L1RecoMatch.h"
#include "DataFormats/TauReco/interface/PFTau.h"

L1RecoMatch::L1RecoMatch(const reco::Candidate* reco, const reco::Candidate* l1,
    const reco::Candidate* l1g, edm::EventID id,
    unsigned int index, unsigned int nTotalObjects, unsigned int nPVs):
  reco_(reco), l1extra_(l1), l1g_(l1g), id_(id),
  index_(index), nTotalObjects_(nTotalObjects), nPVs_(nPVs) {}

const reco::Candidate* L1RecoMatch::reco() const {
  return reco_;
}

const reco::Candidate* L1RecoMatch::l1() const {
  return l1extra_;
}

const reco::Candidate* L1RecoMatch::l1g() const {
  return l1g_;
}

bool L1RecoMatch::l1Match() const {
  return l1extra_ != NULL;
}

bool L1RecoMatch::l1gMatch() const {
  return l1g_ != NULL;
}

const edm::EventID& L1RecoMatch::id() const {
  return id_;
}

unsigned int L1RecoMatch::index() const {
  return index_;
}

unsigned int L1RecoMatch::nTotalObjects() const {
  return nTotalObjects_;
}

unsigned int L1RecoMatch::nPVs() const {
  return nPVs_;
}

double L1RecoMatch::ecalEnergy() const {
  double output = 0;
  const reco::PFTau* the_tau = dynamic_cast<const reco::PFTau*>(this->reco());
  if (!the_tau) {
    return -5;
  }
  for (size_t i = 0; i < the_tau->signalPFCands().size(); ++i) {
    const reco::PFCandidate& pfcand = *(the_tau->signalPFCands().at(i));
    output += pfcand.ecalEnergy() * sin(pfcand.theta());
  }
  return output;
}

double L1RecoMatch::hcalEnergy() const {
  double output = 0;
  const reco::PFTau* the_tau = dynamic_cast<const reco::PFTau*>(this->reco());
  if (!the_tau) {
    return -5;
  }
  for (size_t i = 0; i < the_tau->signalPFCands().size(); ++i) {
    const reco::PFCandidate& pfcand = *(the_tau->signalPFCands().at(i));
    output += pfcand.hcalEnergy() * sin(pfcand.theta());
  }
  return output;
}
