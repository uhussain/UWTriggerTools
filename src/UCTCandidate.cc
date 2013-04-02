#include "../interface/UCTCandidate.h"
#include "L1Trigger/UCT2015/interface/RegionAlgos.h"
#include "FWCore/Utilities/interface/Exception.h"

UCTCandidate::UCTCandidate() : reco::LeafCandidate() {}

UCTCandidate::UCTCandidate(double pt, double eta, double phi, double mass,
    const std::vector<UCTRegion>& regions) :
  reco::LeafCandidate(
      0, reco::LeafCandidate::PolarLorentzVector(pt, eta, phi, mass),
      reco::LeafCandidate::Point(0, 0, 0), 0) {
    // copy over the region information.
    if (regions.size()) {
      hasRegions_ = true;
      std::copy(regions.begin(), regions.end(), std::back_inserter(regions_));
    } else {
      hasRegions_ = false;
    }
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

// version with default value.
template<typename T>
typename T::mapped_type getIfExists(const T& mapping,
    const std::string& key, const typename T::mapped_type& defaultVal) {
  typename T::const_iterator result = mapping.find(key);
  if (result == mapping.end()) {
    return defaultVal;
  }
  return result->second;
}

float UCTCandidate::getFloat(const std::string& item) const {
  return getIfExists(floatData_, item);
}

float UCTCandidate::getFloat(const std::string& item, float defaultVal) const {
  return getIfExists(floatData_, item, defaultVal);
}

int UCTCandidate::getInt(const std::string& item) const {
  return getIfExists(intData_, item);
}

int UCTCandidate::getInt(const std::string& item, int defaultVal) const {
  return getIfExists(intData_, item, defaultVal);
}

const std::string UCTCandidate::getString(const std::string& item) const {
  return getIfExists(stringData_, item);
}

const std::string UCTCandidate::getString(const std::string& item,
    const std::string& defaultVal) const {
  return getIfExists(stringData_, item, defaultVal);
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

bool UCTCandidate::operator < (const UCTCandidate& other) const {
  return this->pt() < other.pt();
}

// Get a region
const UCTRegion& UCTCandidate::getRegion(int etaPos, int phiPos) const {
  if (hasRegions_) {
    for (size_t i = 0; i < regions_.size(); ++i) {
      if (etaPos == regions_[i].etaPos && phiPos == regions_[i].phiPos)
        return regions_[i];
    }
  }
  throw cms::Exception("Region missing") << "The UCT candidate does not have"
    << " a region at (" << etaPos << ", " << phiPos << "), wtf";
}

// Get all regions
const std::vector<UCTRegion>& UCTCandidate::regions() const {
  return regions_;
}

// Set the regions
void UCTCandidate::setRegions(const std::vector<UCTRegion>& in) {
  regions_ = in;
  // sort in ascending PT
  std::sort(regions_.begin(), regions_.end());
}

RegionDiscriminantInfo UCTCandidate::regionDiscriminant(unsigned int N) const {
  return makeDiscriminant(N, regions_);
}

std::ostream& operator<<(std::ostream &os, const UCTCandidate& t) {
  os << "UCTCandidate(" << t.pt()
    << ", " << t.eta() << ", " << t.phi() << ")";
  return os;
}
