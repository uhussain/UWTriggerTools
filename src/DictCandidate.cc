#include "../interface/DictCandidate.h"
#include "FWCore/Utilities/interface/Exception.h"

DictCandidate::DictCandidate() : reco::LeafCandidate() {}

DictCandidate::DictCandidate(int pdgId, int charge, double pt, double eta,
    double phi, double mass) : reco::LeafCandidate(
      charge,
      reco::LeafCandidate::PolarLorentzVector(pt, eta, phi, mass),
      reco::LeafCandidate::Point(0, 0, 0), pdgId) {
}

// Helper to retrieve an item from a std::map.  Throws exception if the key
// doesn't exist.
template<typename T>
typename T::mapped_type getIfExists(const T& mapping, const std::string& key) {
  typename T::const_iterator result = mapping.find(key);
  if (result == mapping.end()) {
    throw cms::Exception("KeyError") << "Can't find the key "
      << key << " in the DictCandidate, sorry." << std::endl;
  }
  return result->second;
}

float DictCandidate::getFloat(const std::string& item) const {
  return getIfExists(floatData_, item);
}

int DictCandidate::getInt(const std::string& item) const {
  return getIfExists(intData_, item);
}

const std::string DictCandidate::getString(const std::string& item) const {
  return getIfExists(stringData_, item);
}

void DictCandidate::setFloat(const std::string& item, float value) {
  floatData_[item] = value;
}

void DictCandidate::setInt(const std::string& item, int value) {
  intData_[item] = value;
}

void DictCandidate::setString(const std::string& item,
    const std::string& value) {
  stringData_[item] = value;
}
