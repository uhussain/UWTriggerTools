/*
 * =====================================================================================
 *
 *       Filename:  Helpers.h
 *
 *    Description:  Common UCT functions.
 *
 *         Author:  M. Cepeda, S. Dasu, E. Friis
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#ifndef HELPERS_W9QK6HND
#define HELPERS_W9QK6HND

// Compute the difference in phi between two towers, wrapping at phi = N
int deltaPhiWrapAtN(unsigned int N, int phi1, int phi2);

// Calo detector mapping.
//
// See: https://twiki.cern.ch/twiki/bin/view/CMS/RCTMap

// Get the physical phi for a given TPG index.
double convertTPGPhi(int iPhi);

// Get the physical eta for a given TPG index
double convertTPGEta(int iEta);

// Convert a region index into physical phi (at center of region)
double convertRegionPhi(int iPhi);

// Convert a region index into physical eta (at center of region)
double convertRegionEta(int iEta);

#endif /* end of include guard: HELPERS_W9QK6HND */
