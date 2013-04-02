#include "L1Trigger/UCT2015/interface/RegionAlgos.h"

RegionCollection getTopNRegions(unsigned int N, const RegionCollection& regions) {
  RegionCollection output;
  for (RegionCollection::const_reverse_iterator region = regions.rbegin();
      region != regions.rend(); ++region) {
    if (!region->et)
      break;
    if (output.size() < N)
      output.push_back(*region);
    else
      break;
  }
  return output;
}

RegionDiscriminantInfo makeDiscriminant(unsigned int NRegions,
    const RegionCollection& regions) {
  RegionCollection topRegions = getTopNRegions(NRegions, regions);
  RegionDiscriminantInfo output;
  output.totalEt = totalEt(topRegions);
  output.totalEtEcal = totalEtEcal(topRegions);
  output.patternPass = matchesTauPattern(topRegions);
  output.numberOfMips = numberOfMips(topRegions);
  output.numberOfRegions = topRegions.size();
  if (output.numberOfRegions) {
    output.lowestRegionEt = topRegions.back().et;
    output.lowestRegionEtEcal = topRegions.back().ecalEt;
  }
  return output;
}

double totalEt(const RegionCollection& regions) {
  double output = 0;
  for (unsigned int i = 0; i < regions.size(); ++i) {
    output += regions[i].et;
  }
  return output;
}

double totalEtEcal(const RegionCollection& regions) {
  double output = 0;
  for (unsigned int i = 0; i < regions.size(); ++i) {
    output += regions[i].ecalEt;
  }
  return output;
}

bool matchesTauPattern(const RegionCollection& regions) {
  if (regions.size() == 1) {
    return true;
  }
  int max_eta_pos = -5;
  int min_eta_pos = +5;
  int max_phi_pos = -5;
  int min_phi_pos = +5;
  for (unsigned int i = 0; i < regions.size(); ++i) {
    const UCTRegion& region = regions[i];
    max_eta_pos = std::max(region.etaPos, max_eta_pos);
    min_eta_pos = std::min(region.etaPos, min_eta_pos);
    max_phi_pos = std::max(region.phiPos, max_phi_pos);
    min_phi_pos = std::min(region.phiPos, min_phi_pos);
  }
  if (max_phi_pos - min_phi_pos > 1)
    return false;
  if (max_eta_pos - min_eta_pos > 1)
    return false;
  if (regions.size() == 2) {
    // if there are two regions, they have to be aligned.
    if (regions[0].etaPos != regions[1].etaPos &&
        regions[0].phiPos != regions[1].phiPos) {
      return false;
    }
  }
  return true;
}

unsigned int numberOfMips(const RegionCollection& regions) {
  unsigned int output = 0;
  for (unsigned int i = 0; i < regions.size(); ++i) {
    output += regions[i].mip;
  }
  return output;
}
