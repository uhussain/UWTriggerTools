#ifndef UCTREGION_HEV5YQ2S
#define UCTREGION_HEV5YQ2S

#include <vector>

// a simple data format to hold region about the associated regions.
struct UCTRegion {
  int etaPos;
  int phiPos;
  double ecalEt;
  double et;
  bool mip;
  bool tauVeto;
  // order by total ET
  bool operator<(const UCTRegion& other) const {
    return this->et < other.et;
  }
};

typedef std::vector<UCTRegion> RegionCollection;

#endif /* end of include guard: UCTREGION_HEV5YQ2S */
