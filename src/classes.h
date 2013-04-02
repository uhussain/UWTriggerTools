/*
 * =====================================================================================
 *
 *    Description:  Define dataformats for the UCT.
 *
 *         Author:  Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#include "L1Trigger/UCT2015/interface/L1RecoMatch.h"
#include "L1Trigger/UCT2015/interface/UCTCandidate.h"
#include "L1Trigger/UCT2015/src/L1GObject.h"

namespace {

  L1RecoMatch dummyMatch;
  L1GObject dummyL1G;
  std::vector<L1GObject> dummyL1GCollection;
  edm::Wrapper<std::vector<L1GObject> > dummyL1GWrapper;

  UCTCandidate dummytDictC;
  std::vector<UCTCandidate> dummyUCTCandCollection;
  edm::Wrapper<std::vector<UCTCandidate> > dummyUCTCandCollectionWrapper;

  UCTRegion dummytRegion;
  std::vector<UCTRegion> dummyUCTRegionCollection;
  RegionDiscriminantInfo dummyDiscInfo;
  std::vector<RegionDiscriminantInfo> dummyDiscInfoVector;

  std::map<std::string, float> dummyMap;

}
