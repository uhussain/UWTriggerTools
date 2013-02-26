#include "L1Trigger/UCT2015/interface/helpers.h"
#include <cmath>
#include <math.h>

int deltaPhiWrapAtN(unsigned int N, int phi1, int phi2) {
  int difference = phi1 - phi2;
  if (std::abs(phi1 - phi2) == N-1) {
    difference = -difference/std::abs(difference);
  }
  return difference;
}

double convertRegionPhi(int iPhi) {
  if (iPhi < 10)
    return 2. * M_PI * iPhi / 18.;
  if (iPhi < 18)
    return -M_PI + 2. * M_PI * (iPhi - 9) / 18.;
  return -9;
}

double convertTPGPhi(int iPhi) {
  if (iPhi < 37)
    return 2. * M_PI * iPhi / 72.;
  if (iPhi < 72)
    return -M_PI + 2. * M_PI * (iPhi - 36) / 72.;
  return -9;
}

double convertRegionEta(int iEta) {
  const double rgnEtaValues[11] = {
     0.174, // HB and inner HE bins are 0.348 wide
     0.522,
     0.870,
     1.218,
     1.566,
     1.956, // Last two HE bins are 0.432 and 0.828 wide
     2.586,
     3.250, // HF bins are 0.5 wide
     3.750,
     4.250,
     4.750
  };
  if(iEta < 11) {
    return -rgnEtaValues[-(iEta - 10)]; // 0-10 are negative eta values
  }
  else if (iEta < 22) {
    return rgnEtaValues[iEta - 11];     // 11-21 are positive eta values
  }
  return -9;
}

double convertTPGEta(int iEta) {
  // 0-27 are negative eta values, 28-55 are positive eta values
  int correctedIndex = iEta < 28 ? -(iEta - 27) : iEta - 28;

  double etaValue = -9;
  if (correctedIndex < 20) {
    etaValue = 0.0435 + correctedIndex * 0.087;
  } else {
    const double endcapEtaValues[8] = {
       1.785,
       1.880,
       1.9865,
       2.1075,
       2.247,
       2.411,
       2.575,
       2.825
    };
    etaValue = endcapEtaValues[correctedIndex-20];
  }
  if (iEta < 28)
    return -etaValue;
  else
    return etaValue;
}

double getRegionArea(int gctEta) {
  switch (gctEta) {
    case 0: return 0.5*0.348;
    case 1: return 0.5*0.348;
    case 2: return 0.5*0.348;
    case 3: return 0.5*0.348;
    case 4: return 0.828*0.348;
    case 5: return 0.432*0.348;
    case 6: return 0.348*0.348;
    case 7: return 0.348*0.348;
    case 8: return 0.348*0.348;
    case 9: return 0.348*0.348;
    case 10: return 0.348*0.348;
    case 21: return 0.5*0.348;
    case 20: return 0.5*0.348;
    case 19: return 0.5*0.348;
    case 18: return 0.5*0.348;
    case 17: return 0.828*0.348;
    case 16: return 0.432*0.348;
    case 15: return 0.348*0.348;
    case 14: return 0.348*0.348;
    case 13: return 0.348*0.348;
    case 12: return 0.348*0.348;
    case 11: return 0.348*0.348;
  }
  return 0;
}

int twrPhi2RegionPhi(int iPhi) {
  unsigned int rgnIdx = (iPhi + 2)/4;
  // 70 and 71 are actually in GCT phi 0
  if (rgnIdx == 18)
    return 0;
  return rgnIdx;
}

int twrEta2RegionEta(int iEta) {
  // 4 towers per region.  Non-HF regions start at ieta=4 in RCT land.
  unsigned int rgnIdx = (iEta / 4) + 4;
  return rgnIdx;
}
