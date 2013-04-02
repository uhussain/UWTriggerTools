#ifndef REGIONALGOS_ASQME1L9
#define REGIONALGOS_ASQME1L9
/*
 * Algorithms which extract information about regions associated
 *
 */
#include "L1Trigger/UCT2015/interface/UCTRegion.h"

struct RegionDiscriminantInfo {
  double totalEt;
  double totalEtEcal;
  double lowestRegionEt;
  double lowestRegionEtEcal;
  int numberOfMips;
  int numberOfRegions;
  bool patternPass;
};

// Compute composite discriminant info for N regions
RegionDiscriminantInfo makeDiscriminant(unsigned int NRegions,
    const RegionCollection& regions);

// Get the highest N regions by total energy.  We expect the input regions to
// be ordered by ascending pt.
//
RegionCollection getTopNRegions(unsigned int N, const RegionCollection& regions);

double totalEt(const RegionCollection&);

double totalEtEcal(const RegionCollection&);

// Check whether the regions match the corresponding tau pattern
bool matchesTauPattern(const RegionCollection&);

// Check whether the regions match the tau pattern, i.e. the bits are arranged
//
// XX  XX XX X X
// XX  X     X

// Number of bits which pass MIP
unsigned int numberOfMips(const RegionCollection&);

#endif /* end of include guard: REGIONALGOS_ASQME1L9 */
