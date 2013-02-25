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
