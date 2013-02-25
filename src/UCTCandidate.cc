#include "../interface/UCTCandidate.h"
#include "FWCore/Utilities/interface/Exception.h"

UCTCandidate::UCTCandidate() : reco::LeafCandidate() {}

UCTCandidate::UCTCandidate(double pt, double eta, double phi, double mass) :
  reco::LeafCandidate(
      0, reco::LeafCandidate::PolarLorentzVector(pt, eta, phi, mass),
      reco::LeafCandidate::Point(0, 0, 0), 0) {
}

// Helper to retrieve an item from a std::map.  Throws exception if the key
// doesn't exist.
template<typename T>
typename T::mapped_type getIfExists(const T& mapping, const std::string& key) {
  typename T::const_iterator result = mapping.find(key);
  if (result == mapping.end()) {
    throw cms::Exception("KeyError") << "Can't find the key "
      << key << " in the UCTCandidate, sorry." << std::endl;
  }
  return result->second;
}

float UCTCandidate::getFloat(const std::string& item) const {
  return getIfExists(floatData_, item);
}

int UCTCandidate::getInt(const std::string& item) const {
  return getIfExists(intData_, item);
}

const std::string UCTCandidate::getString(const std::string& item) const {
  return getIfExists(stringData_, item);
}

void UCTCandidate::setFloat(const std::string& item, float value) {
  floatData_[item] = value;
}

void UCTCandidate::setInt(const std::string& item, int value) {
  intData_[item] = value;
}

void UCTCandidate::setString(const std::string& item,
    const std::string& value) {
  stringData_[item] = value;
}
