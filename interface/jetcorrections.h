#ifndef JETCORRECTIONS_MSO6R8RU
#define JETCORRECTIONS_MSO6R8RU

/*
 * Jet corrections from Mike, Florent and friends.
 * The constants are in the cc file.
 *
 * Factored out by Evan.
 *
 */

// Version using normal PU subtraction algorithm
double jetcorr(double pt, int ieta, unsigned int pileup);

// Version using UIC PU eta-dependent correction method.
double jetcorrUIC(double pt, int ieta, unsigned int pileup);

#endif
