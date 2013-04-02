#ifndef DICTCANDIDATE_GNUPVDDA
#define DICTCANDIDATE_GNUPVDDA

/*
 * =====================================================================================
 *
 *       Filename:  UCTCandidate.h
 *
 *    Description:  A simple candidate which has arbitrary float/int/string
 *                  data attached to it.
 *
 *         Author:  Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include <map>
#include <string>
#include <iostream>

#include "L1Trigger/UCT2015/interface/UCTRegion.h"
#include "L1Trigger/UCT2015/interface/RegionAlgos.h"

class UCTCandidate : public reco::LeafCandidate {
  public:
    UCTCandidate();
    UCTCandidate(double pt, double eta, double phi, double mass=0,
        const std::vector<UCTRegion>& regions=std::vector<UCTRegion>());


    // Data attribute retrieval, throwing exceptions if attr is undefined.
    float getFloat(const std::string& item) const;
    int getInt(const std::string& item) const;
    const std::string getString(const std::string& item) const;

    // Data attribute retrieval with a default value if the attr doesn't exist.
    float getFloat(const std::string& item, float defaultVal) const;
    int getInt(const std::string& item, int defaultVal) const;
    const std::string getString(const std::string& item,
        const std::string& defaultVal) const;

    // Data attribute definition
    void setFloat(const std::string& item, float value);
    void setInt(const std::string& item, int value);
    void setString(const std::string& item, const std::string& value);

    // Sort by ascending PT per default.
    bool operator < (const UCTCandidate& other) const;

    // Get a region.
    const UCTRegion& getRegion(int etaPos, int phiPos) const;
    const std::vector<UCTRegion>& regions() const;
    // Set the regions
    void setRegions(const std::vector<UCTRegion>&  in);

    // Get discriminant into for patterns with N objects
    RegionDiscriminantInfo regionDiscriminant(unsigned int N) const;

    friend std::ostream& operator<<(std::ostream &os, const UCTCandidate& t);

  private:
    // data storage
    std::map<std::string, float> floatData_;
    std::map<std::string, int> intData_;
    std::map<std::string, std::string> stringData_;

    // an array of the associated regions.
    bool hasRegions_;
    std::vector<UCTRegion> regions_;
};

#endif /* end of include guard: DICTCANDIDATE_GNUPVDDA */
