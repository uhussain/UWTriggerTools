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

class UCTCandidate : public reco::LeafCandidate {
  public:
    UCTCandidate();
    UCTCandidate(double pt, double eta, double phi, double mass=0);

    // Data attribute retrieval
    float getFloat(const std::string& item) const;
    int getInt(const std::string& item) const;
    const std::string getString(const std::string& item) const;

    // Data attribute definition
    void setFloat(const std::string& item, float value);
    void setInt(const std::string& item, int value);
    void setString(const std::string& item, const std::string& value);

    // Sort by PT per default.
    bool operator < (const UCTCandidate& other) {
      return this->pt() < other.pt();
    }

  private:
    // data storage
    std::map<std::string, float> floatData_;
    std::map<std::string, int> intData_;
    std::map<std::string, std::string> stringData_;
};

#endif /* end of include guard: DICTCANDIDATE_GNUPVDDA */
