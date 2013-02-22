#ifndef UCT2015L1ExtraProducer_h
#define UCT2015L1ExtraProducer_h
// -*- C++ -*-
//
// Package:     L1ExtraFromDigis
// Class  :     L1ExtraParticlesProd
// 
/**\class L1ExtraParticlesProd \file L1ExtraParticlesProd.h L1Trigger/L1ExtraFromDigis/interface/L1ExtraParticlesProd.h \author Werner Sun

 Description: producer of L1Extra particle objects from Level-1 hardware objects.

*/
//
// Original Author:  
//         Created:  Tue Oct 17 00:13:51 EDT 2006
// $Id: UCT2015L1ExtraProducer.h,v 1.2 2013/01/15 17:50:58 jbrooke Exp $
//

// system include files

// user include files
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

// for muon isolation
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "L1Trigger/UCT2015/src/L1GObject.h"


class UCT2015L1ExtraProducer : public edm::EDProducer {
   public:
      explicit UCT2015L1ExtraProducer(const edm::ParameterSet&);
      ~UCT2015L1ExtraProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      int deltaPhi18(int phi1,
		     int phi2);
      int deltaGctPhi(const L1CaloRegion& r1,
		      const L1CaloRegion& r2);
      L1CaloRegionCollection::const_iterator matchObjectToRegion(const L1CaloGeometry* geom,
								 const edm::Handle<L1CaloRegionCollection> &newRegions,
								 const float &eta, const float &phi);
      L1GObject buildJetAtIndex(const edm::Handle<L1CaloRegionCollection> &newRegions,
				const L1CaloRegionCollection::const_iterator &newRegion,
				const float &regionLSB_,
				float threshold);

      // ----------member data ---------------------------
      bool produceMuonParticles_ ;
      edm::InputTag muonSource_ ;

      bool produceCaloParticles_ ;
      edm::InputTag egSource_ ;
      edm::InputTag tauSource_ ;
      edm::InputTag jetSource_ ;
      edm::InputTag etMissSource_ ;
      edm::InputTag etTotSource_ ;
      edm::InputTag htMissSource_ ;
      edm::InputTag htTotSource_ ;

      static double muonMassGeV_ ;

      bool centralBxOnly_ ;

      // EG & tau isolation
      double egIso_;  // relative
      double tauIso_; // relative

      // muon isolation
      edm::InputTag regionSource_ ;
      double rgnThreshold_;
      double regionLSB_;
      double muFixedIso_;  // require 3x3 Et below this value
      double muRelIso_;    // placeholder for relative isolation

};

#endif
