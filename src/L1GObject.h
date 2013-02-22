#ifndef L1GObject_h
#define L1GObject_h

#include <iostream>
using std::iostream;
using std::ostream;

#include <string>
using std::string;

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

/**
 * L1GObject represents a calorimeter global trigger object that
 * is made from global calorimeter quantities, total ET, missing ET.
 * These objects are created/filled while processing the calorimeter
 * trigger information at the card level.
 */

class L1GObject : public reco::LeafCandidate
{
public:

  // Constructors

  L1GObject() {}

  L1GObject(unsigned int et, unsigned int eta, unsigned int phi, string name = "L1GObject", bool twrGranularity = false)
    : myEt(et), myEta(eta), myPhi(phi), myName(name), myTwrGranularity(twrGranularity) {initialize();}
    
  L1GObject(unsigned int packedObject, string name = "L1GObject", bool twrGranularity = false) {
    myEt = (packedObject & 0xFFFF0000) >> 16;
    myEta = (packedObject & 0x0000FF00) >> 8;
    myPhi = (packedObject & 0x000000FF);
    myName = name;
    myTwrGranularity = twrGranularity;
    initialize();
  }

  unsigned int packedObject() {
    if(myEt > 0xFFFF) myEt = 0xFFFF;
    unsigned int etBits = (myEt << 16);
    if(myEta < 0xFF) {
      unsigned int etaBits = (myEta << 8);
      if(myPhi < 0xFF) {
	return (etBits + etaBits + myPhi);
      }
    }
    std::cerr << "L1GObject: Cannot pack content - fatal error: " << myEt << ", " << myEta << ", " << myPhi << std::endl;
    return (etBits);
  }

  L1GObject(const L1GObject& t)
    {
      myTwrGranularity = t.myTwrGranularity;
      myName = t.myName;
      myPhi = t.myPhi;
      myEta = t.myEta;
      myEt = t.myEt;
      associatedRegionEt_ = t.associatedRegionEt_;
      associatedJetPt_ = t.associatedJetPt_;
      ellIsolation_ = t.ellIsolation_;
      puLevel_ = t.puLevel_;
      puLevelUIC_ = t.puLevelUIC_;
      tauVeto_ = t.tauVeto_;
      mipBit_ = t.mipBit_;
      initialize();
    }

  L1GObject& operator=(const L1GObject& t)
    {
      if(this != &t)
	{
	  myTwrGranularity = t.myTwrGranularity;
	  myName = t.myName;
	  myPhi = t.myPhi;
	  myEta = t.myEta;
	  myEt = t.myEt;
          associatedRegionEt_ = t.associatedRegionEt_;
          associatedJetPt_ = t.associatedJetPt_;
          ellIsolation_ = t.ellIsolation_;
          puLevel_ = t.puLevel_;
          puLevelUIC_ = t.puLevelUIC_;
          tauVeto_ = t.tauVeto_;
          mipBit_ = t.mipBit_;
	  initialize();
	}
      return *this;
    }

  // Destructor

  virtual ~L1GObject() {}

  // Access functions

  string name() const {return myName;}

  bool empty() const {return false;}

  double ptValue() const {
    return myLSB * myEt;
  }

  double etaValue() const {
    if(myTwrGranularity) {
      if(myEta < 28) {
	return -twrEtaValues[-(myEta - 27)]; // 0-27 are negative eta values
      }
      else if(myEta < 56) {
	return twrEtaValues[myEta - 28];     // 28-55 are positive eta values
      }
    }
    else {
      if(myEta < 11) {
	return -rgnEtaValues[-(myEta - 10)]; // 0-10 are negative eta values
      }
      else if (myEta < 22) {
	return rgnEtaValues[myEta - 11];     // 11-21 are positive eta values
      }
    }
    return 999.;
  }

  double phiValue() const {
    if(myTwrGranularity) {
      if(myPhi < 72)
	return twrPhiValues[myPhi];
    }
    else {
      if(myPhi < 18)
	return rgnPhiValues[myPhi];
    }
    return 999.;
  }

  unsigned int ptCode() const {return myEt;}

  unsigned int etaIndex() const {
    if(myTwrGranularity) {
      return myEta / 4;
    }
    return myEta;
  }

  unsigned int phiIndex() const {
    if(myTwrGranularity) {
      return myPhi / 4;
    }
    return myPhi;
  }

  unsigned int etaTwrIndex() const {
    if(myTwrGranularity) {
      return myEta;
    }
    return myEta * 4 - 2;
  }

  unsigned int phiTwrIndex() const {
    if(myTwrGranularity) {
      return myPhi;
    }
    return myPhi * 4 - 2;
  }

  bool twrGranularity() const {return myTwrGranularity;}

  // Operators required for sorting lists of these objects

  bool operator==(const L1GObject& t) const
    {
      if(myEt == t.myEt) return true;
      else return false;
    }

  bool operator<(const L1GObject& t) const
    {
      if(myEt < t.myEt) return true;
      else return false;
    }

  bool operator>(const L1GObject& t) const
    {
      if(myEt > t.myEt) return true;
      else return false;
    }

  bool operator<=(const L1GObject& t) const
    {
      if(myEt <= t.myEt) return true;
      else return false;
    }

  bool operator>=(const L1GObject& t) const
    {
      if(myEt >= t.myEt) return true;
      else return false;
    }

  friend ostream& operator<<(ostream &os, const L1GObject& t)
    {
      os << "L1GObject : Name = " << t.name()
	 << "(Et, Eta, Phi) = ("
	 << t.myEt << ", "
	 << t.myEta << ", "
	 << t.myPhi << ") ("
	 << t.ptValue() << ", "
	 << t.etaValue() << ", "
	 << t.phiValue() << ")";
      return os;
    }

  void setEt(unsigned int et) {myEt = et;}
  void setEta(unsigned int eta) {myEta = eta;}
  void setPhi(unsigned int phi) {myPhi = phi;}
  void setName(string name) {myName = name;}
  void setLSB(double lsb) {myLSB = lsb;}
  void setTwrGranularity(bool twrGranularity) {myTwrGranularity = twrGranularity;}

  void initialize()
  {
    for(unsigned int i = 0; i < 10; i++) {
      rgnPhiValues[i] = 2. * 3.1415927 * i / 18;
    }
    for(unsigned int j = 10; j < 18; j++) {
      rgnPhiValues[j] = -3.1415927 + 2. * 3.1415927 * (j - 9) / 18;
    }
    for(unsigned int i = 0; i < 37; i++) {
      twrPhiValues[i] = 2. * 3.1415927 * i / 72;
    }
    for(unsigned int j = 37; j < 72; j++) {
      twrPhiValues[j] = -3.1415927 + 2. * 3.1415927 * (j - 36) / 72;
    }
    rgnEtaValues[ 0] = 0.174; // HB and inner HE bins are 0.348 wide
    rgnEtaValues[ 1] = 0.522;
    rgnEtaValues[ 2] = 0.870;
    rgnEtaValues[ 3] = 1.218;
    rgnEtaValues[ 4] = 1.566;
    rgnEtaValues[ 5] = 1.956; // Last two HE bins are 0.432 and 0.828 wide
    rgnEtaValues[ 6] = 2.586;
    rgnEtaValues[ 7] = 3.250; // HF bins are 0.5 wide
    rgnEtaValues[ 8] = 3.750;
    rgnEtaValues[ 9] = 4.250;
    rgnEtaValues[10] = 4.750;
    for(unsigned int i = 0; i < 20; i++) {
      twrEtaValues[i] = 0.0435 + i * 0.087;
    }
    twrEtaValues[20] = 1.785;
    twrEtaValues[21] = 1.880;
    twrEtaValues[22] = 1.9865;
    twrEtaValues[23] = 2.1075;
    twrEtaValues[24] = 2.247;
    twrEtaValues[25] = 2.411;
    twrEtaValues[26] = 2.575;
    twrEtaValues[27] = 2.825;
    myLSB = 1.0;

    // Initialize tuning parameters
    //associatedRegionEt = -1;
    //associatedJetPt = -1;
    //puLevel = -1;

    // Setup the reco::Candidate (physics) 4-vector
    math::PtEtaPhiMLorentzVector myP4(
        this->ptValue(), this->etaValue(), this->phiValue(), 0);
    this->setP4(myP4);
  }

  // Extra values for tuning UCT parameters - just public members, to
  // eventually be removed
  double associatedJetPt() const { return associatedJetPt_; }
  unsigned int puLevel() const { return puLevel_; }
  //double puLevelUIC() const { return puLevelUIC_; }
  unsigned int puLevelUIC() const { return puLevelUIC_; }
  double associatedRegionEt() const { return associatedRegionEt_; }
  bool ellIsolation() const { return ellIsolation_; };
  bool tauVeto() const { return tauVeto_; }
  bool mipBit() const { return mipBit_; }

  double associatedJetPt_;
  unsigned int puLevel_;
  //double puLevelUIC_;
  unsigned int  puLevelUIC_;
  double associatedRegionEt_;
  bool ellIsolation_;

  // For the EG objects, don't require this to build the object, just embed it.
  bool tauVeto_;
  bool mipBit_;

private:

  unsigned int myEt;
  unsigned int myEta;
  unsigned int myPhi;

  string myName;

  bool myTwrGranularity;

  double myLSB;
  double rgnEtaValues[11];
  double rgnPhiValues[18];
  double twrEtaValues[28];
  double twrPhiValues[72];

};

#endif
