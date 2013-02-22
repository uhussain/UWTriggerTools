// -*- C++ -*-
//
// Package:    UCT2015L1ExtraProducer
// Class:      UCT2015L1ExtraProducer
// 
/**\class UCT2015L1ExtraProducer \file UCT2015L1ExtraProducer.cc src/UCT2015L1ExtraProducer/src/UCT2015L1ExtraProducer.cc
*/
//
// Original Author:  Werner Sun
//         Created:  Mon Oct  2 22:45:32 EDT 2006
// $Id: UCT2015L1ExtraProducer.cc,v 1.3 2013/01/15 17:50:58 jbrooke Exp $
//
//


// system include files
#include <memory>

// user include files
#include "L1Trigger/UCT2015/plugins/UCT2015L1ExtraProducer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

// for muon isolation
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "L1Trigger/UCT2015/src/L1GObject.h"

// #include "FWCore/Utilities/interface/EDMException.h"

//
// class decleration
//


//
// constants, enums and typedefs
//

//
// static data member definitions
//

double UCT2015L1ExtraProducer::muonMassGeV_ = 0.105658369 ; // PDG06

//
// constructors and destructor
//
UCT2015L1ExtraProducer::UCT2015L1ExtraProducer(const edm::ParameterSet& iConfig)
   : produceMuonParticles_( iConfig.getParameter< bool >(
      "produceMuonParticles" ) ),
     muonSource_( iConfig.getParameter< edm::InputTag >(
	"muonSource" ) ),
     produceCaloParticles_( iConfig.getParameter< bool >(
	"produceCaloParticles" ) ),
     egSource_( iConfig.getParameter< edm::InputTag >("egLabel") ),
     tauSource_( iConfig.getParameter< edm::InputTag >("tauLabel") ),
     jetSource_( iConfig.getParameter< edm::InputTag >("jetLabel") ),
     etMissSource_( iConfig.getParameter< edm::InputTag >("etMissLabel") ),
     etTotSource_( iConfig.getParameter< edm::InputTag >("etTotLabel") ),
     htMissSource_( iConfig.getParameter< edm::InputTag >("htMissLabel") ),
     htTotSource_( iConfig.getParameter< edm::InputTag >("htTotLabel") ),
     centralBxOnly_( iConfig.getParameter< bool >( "centralBxOnly" ) ),
     egIso_( iConfig.getParameter< double >( "egIso" ) ),
     tauIso_( iConfig.getParameter< double >( "tauIso" ) ),
     regionSource_( iConfig.getParameter< edm::InputTag >("regionLabel") ),
     rgnThreshold_( iConfig.getParameter<double>("rgnThreshold") ),
     regionLSB_( iConfig.getParameter<double>("regionLSB") ),
     muFixedIso_( iConfig.getParameter<double>("muFixedIso") ),
     muRelIso_( iConfig.getParameter<double>("muRelIso") )
{
   using namespace l1extra ;

   //register your products
   produces< L1EmParticleCollection >( "Isolated" ) ;
   produces< L1EmParticleCollection >( "Relaxed" ) ;
   produces< L1JetParticleCollection >( "Jets" ) ;
   produces< L1JetParticleCollection >( "FwdJets" ) ;
   produces< L1JetParticleCollection >( "IsolatedTau" ) ;
   produces< L1JetParticleCollection >( "RelaxedTau" ) ;
   produces< L1MuonParticleCollection >( "" ) ;
   //   produces< L1MuonParticleCollection >( "Isolated" ) ;
   produces< L1EtMissParticleCollection >( "MET" ) ;
   produces< L1EtMissParticleCollection >( "MHT" ) ;

   //now do what ever other initialization is needed
}


UCT2015L1ExtraProducer::~UCT2015L1ExtraProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
UCT2015L1ExtraProducer::produce( edm::Event& iEvent,
			       const edm::EventSetup& iSetup)
{
   using namespace edm ;
   using namespace l1extra ;
   using namespace std ;

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~ Muons ~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   auto_ptr< L1MuonParticleCollection > muColl(
     new L1MuonParticleCollection );

   //   auto_ptr< L1MuonParticleCollection > isoMuColl(
   //     new L1MuonParticleCollection );

   if( produceMuonParticles_ ) {
      ESHandle< L1MuTriggerScales > muScales ;
      iSetup.get< L1MuTriggerScalesRcd >().get( muScales ) ;

      ESHandle< L1MuTriggerPtScale > muPtScale ;
      iSetup.get< L1MuTriggerPtScaleRcd >().get( muPtScale ) ;

      Handle< L1MuGMTReadoutCollection > hwMuCollection ;
      iEvent.getByLabel( muonSource_, hwMuCollection ) ;

      vector< L1MuGMTExtendedCand > hwMuCands ;

      // for isolation
      edm::Handle<L1CaloRegionCollection> regions;
      iEvent.getByLabel(regionSource_, regions);
      
      ESHandle< L1CaloGeometry > caloGeomESH ;
      iSetup.get< L1CaloGeometryRecord >().get( caloGeomESH ) ;
      const L1CaloGeometry* caloGeom = &( *caloGeomESH ) ;

      if( !hwMuCollection.isValid() )
	{
	  LogDebug("UCT2015L1ExtraProducer")
	    << "\nWarning: L1MuGMTReadoutCollection with " << muonSource_
	    << "\nrequested in configuration, but not found in the event."
	    << std::endl;
	}
      else
	{
	  if( centralBxOnly_ )
	    {
	      // Get GMT candidates from central bunch crossing only
	      hwMuCands = hwMuCollection->getRecord().getGMTCands() ;
	    }
	  else
	    {
	      // Get GMT candidates from all bunch crossings
	      vector< L1MuGMTReadoutRecord > records = hwMuCollection->getRecords();
	      vector< L1MuGMTReadoutRecord >::const_iterator rItr = records.begin();
	      vector< L1MuGMTReadoutRecord >::const_iterator rEnd = records.end();

	      for( ; rItr != rEnd ; ++rItr )
		{
		  vector< L1MuGMTExtendedCand > tmpCands = rItr->getGMTCands() ;

		  hwMuCands.insert( hwMuCands.end(),
				    tmpCands.begin(),
				    tmpCands.end() ) ;
		}
	    }

//       cout << "HW muons" << endl ;
	  vector< L1MuGMTExtendedCand >::const_iterator muItr = hwMuCands.begin() ;
	  vector< L1MuGMTExtendedCand >::const_iterator muEnd = hwMuCands.end() ;
	  for( int i = 0 ; muItr != muEnd ; ++muItr, ++i )
	    {
// 	 cout << "#" << i
// 	      << " name " << muItr->name()
// 	      << " empty " << muItr->empty()
// 	      << " pt " << muItr->ptIndex()
// 	      << " eta " << muItr->etaIndex()
// 	      << " phi " << muItr->phiIndex()
// 	      << " iso " << muItr->isol()
// 	      << " mip " << muItr->mip()
// 	      << " bx " << muItr->bx()
// 	      << endl ;

	      if( !muItr->empty() )
		{
		  // keep x and y components non-zero and protect against roundoff.
		  double pt =
		    muPtScale->getPtScale()->getLowEdge( muItr->ptIndex() ) + 1.e-6 ;

// 	    cout << "L1Extra pt " << pt << endl ;

		  double eta =
		    muScales->getGMTEtaScale()->getCenter( muItr->etaIndex() ) ;

		  double phi =
		    muScales->getPhiScale()->getLowEdge( muItr->phiIndex() ) ;

		  math::PtEtaPhiMLorentzVector p4( pt,
						   eta,
						   phi,
						   muonMassGeV_ ) ;

		  // calculate isolation
		  bool isolated = false;

		  L1CaloRegionCollection::const_iterator matchedRegion = 
		    matchObjectToRegion(caloGeom,
					regions,
					eta,
					phi);
		  L1GObject seededJet = buildJetAtIndex(regions,
							matchedRegion,
							regionLSB_,
							rgnThreshold_);

		  isolated = (muFixedIso_ > 0.) && (seededJet.pt() < muFixedIso_);

		  L1MuonParticle muonParticle( muItr->charge(),
					       p4,
					       *muItr,
					       muItr->bx() );

		  muonParticle.setIsolated(isolated);
		    
		  // push to collection
		  //if (isolated) isoMuColl->push_back( muonParticle );
		  muColl->push_back( muonParticle );

		}
	    }
	}
   }
   
   OrphanHandle< L1MuonParticleCollection > muHandle =
     iEvent.put( muColl );

   //   OrphanHandle< L1MuonParticleCollection > isoMuHandle =
   //     iEvent.put( isoMuColl, "Isolated" );


   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~ Calorimeter ~~~~~~~~~~~~~~~~~~~~
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   auto_ptr< L1EmParticleCollection > isoEGColl(
      new L1EmParticleCollection );

   auto_ptr< L1EmParticleCollection > relEGColl(
      new L1EmParticleCollection );

   auto_ptr< L1JetParticleCollection > isoTauColl(
      new L1JetParticleCollection );

   auto_ptr< L1JetParticleCollection > relTauColl(
      new L1JetParticleCollection );

   auto_ptr< L1JetParticleCollection > jetColl(
      new L1JetParticleCollection );

   auto_ptr< L1JetParticleCollection > fwdJetColl(
      new L1JetParticleCollection );

   auto_ptr< L1EtMissParticleCollection > etMissColl(
      new L1EtMissParticleCollection );

   auto_ptr< L1EtMissParticleCollection > htMissColl(
      new L1EtMissParticleCollection );

   if( produceCaloParticles_ )
   {

      // EG
      Handle< std::vector<L1GObject> > egObjs ;
      iEvent.getByLabel( egSource_, egObjs ) ;

      if( !egObjs.isValid() ) {
	LogDebug("UCT2015L1ExtraProducer")
	  << "\nWarning: std::vector<L1GObject> with " << egSource_
	  << "\nrequested in configuration, but not found in the event."
	  << std::endl;
      }
      else {
	
	std::vector<L1GObject>::const_iterator itr = egObjs->begin() ;
	std::vector<L1GObject>::const_iterator end = egObjs->end() ;
	for( int i = 0 ; itr != end ; ++itr, ++i ) {
	  double pt = itr->ptValue();
	  double jetPt = itr->associatedJetPt();
	  double phi = itr->phiValue();
	  double eta = itr->etaValue();
	  math::PtEtaPhiMLorentzVector p4( pt,
					   eta,
					   phi,
					   0. );
	  Ref< L1GctEmCandCollection > ref();

	  // check relaxed
	  if ( (!itr->tauVeto() && !itr->mipBit()) || pt > 62) {

	    relEGColl->push_back( L1EmParticle( p4,
						Ref< L1GctEmCandCollection >(), 
						0 ) );

	  // check isolated
	    if ( (jetPt - pt)/pt < egIso_ ) {
	      isoEGColl->push_back( L1EmParticle( p4,
						  Ref< L1GctEmCandCollection >(), 
						  0 ) );
	    }

	  }

	}
      }
      

      // tau
      Handle< std::vector<L1GObject> > tauObjs ;
      iEvent.getByLabel( tauSource_, tauObjs ) ;

      if( !tauObjs.isValid() ) {
	LogDebug("UCT2015L1ExtraProducer")
	  << "\nWarning: std::vector<L1GObject> with " << tauSource_
	  << "\nrequested in configuration, but not found in the event."
	  << std::endl;
      }
      else {
	
	std::vector<L1GObject>::const_iterator itr = tauObjs->begin() ;
	std::vector<L1GObject>::const_iterator end = tauObjs->end() ;
	for( int i = 0 ; itr != end ; ++itr, ++i ) {
	  double pt = max(itr->pt(), itr->associatedRegionEt());
	  double jetPt = itr->associatedJetPt();
	  double phi = itr->phiValue();
	  double eta = itr->etaValue();
	  math::PtEtaPhiMLorentzVector p4( pt,
					   eta,
					   phi,
					   0. );

	  // check relaxed
	  if ( true ) {  // take any tau as relaxed !

	    relTauColl->push_back( L1JetParticle( p4, 
						  Ref< L1GctJetCandCollection >(), 
						  0 ) );

	  // check isolated
	    if ( (jetPt - pt)/pt < tauIso_ ) {
	      isoTauColl->push_back( L1JetParticle( p4, 
						    Ref< L1GctJetCandCollection >(), 
						    0 ) );
	    }

	  }

	}
      }

      // jet
      Handle< std::vector<L1GObject> > jetObjs ;
      iEvent.getByLabel( jetSource_, jetObjs ) ;
      
      if( !jetObjs.isValid() ) {
	LogDebug("UCT2015L1ExtraProducer")
	  << "\nWarning: std::vector<L1GObject> with " << jetSource_
	  << "\nrequested in configuration, but not found in the event."
	  << std::endl;
      }
      else {
	
	std::vector<L1GObject>::const_iterator itr = jetObjs->begin() ;
	std::vector<L1GObject>::const_iterator end = jetObjs->end() ;
	for( int i = 0 ; itr != end ; ++itr, ++i ) {
	  double pt = itr->pt();
	  double phi = itr->phiValue();
	  double eta = itr->etaValue();
	  math::PtEtaPhiMLorentzVector p4( pt,
					   eta,
					   phi,
					   0. );
	  
	  jetColl->push_back( L1JetParticle( p4, 
					     Ref< L1GctJetCandCollection >(), 
					     0 ) );
	  
	}
      }
      
      /// sums
      Handle< std::vector<L1GObject> > metObjs ;
      iEvent.getByLabel( etMissSource_, metObjs ) ;

      Handle< std::vector<L1GObject> > setObjs ;
      iEvent.getByLabel( etTotSource_, setObjs ) ;

      if( !metObjs.isValid() ) {
        LogDebug("UCT2015L1ExtraProducer")
          << "\nWarning: std::vector<L1GObject> with " << etMissSource_
          << "\nrequested in configuration, but not found in the event."
          << std::endl;
      }
      else if ( !setObjs.isValid() ) {
        LogDebug("UCT2015L1ExtraProducer")
          << "\nWarning: std::vector<L1GObject> with " << etTotSource_
          << "\nrequested in configuration, but not found in the event."
          << std::endl;
      }
      else {

	if (metObjs->size()>0 && setObjs->size()>0) {

	  L1GObject met = metObjs->at(0);
	  L1GObject set = setObjs->at(0);

	  math::PtEtaPhiMLorentzVector p4( met.ptValue(),
                                           0.,
                                           met.phiValue(),
                                           0. );

	  double etTot = set.ptValue();

	  etMissColl->push_back( L1EtMissParticle( p4,
						L1EtMissParticle::kMET,
						etTot,
						Ref< L1GctEtMissCollection >(),
						Ref< L1GctEtTotalCollection >(),
						Ref< L1GctHtMissCollection >(),
						Ref< L1GctEtHadCollection >(),
						0 ) );

        }
      }

      Handle< std::vector<L1GObject> > mhtObjs ;
      iEvent.getByLabel( htMissSource_, mhtObjs ) ;
      
      Handle< std::vector<L1GObject> > shtObjs ;
      iEvent.getByLabel( htTotSource_, shtObjs ) ;

      if( !mhtObjs.isValid() ) {
        LogDebug("UCT2015L1ExtraProducer")
          << "\nWarning: std::vector<L1GObject> with " << htMissSource_
          << "\nrequested in configuration, but not found in the event."
          << std::endl;
      }
      else if ( !shtObjs.isValid() ) {
        LogDebug("UCT2015L1ExtraProducer")
          << "\nWarning: std::vector<L1GObject> with " << htTotSource_
          << "\nrequested in configuration, but not found in the event."
          << std::endl;
      }
      else {

	if (mhtObjs->size()>0 && shtObjs->size()>0) {

	  L1GObject mht = mhtObjs->at(0);
	  L1GObject sht = shtObjs->at(0);

	  math::PtEtaPhiMLorentzVector p4( mht.ptValue(),
                                           0.,
                                           mht.phiValue(),
                                           0. );

	  double htTot = sht.ptValue();

	  htMissColl->push_back( L1EtMissParticle( p4,
						L1EtMissParticle::kMET,
						htTot,
						Ref< L1GctEtMissCollection >(),
						Ref< L1GctEtTotalCollection >(),
						Ref< L1GctHtMissCollection >(),
						Ref< L1GctEtHadCollection >(),
						0 ) );

        }
      }

   }

   OrphanHandle< L1EmParticleCollection > isoEGHandle =
     iEvent.put( isoEGColl, "Isolated" ) ;

   OrphanHandle< L1EmParticleCollection > relEGHandle =
     iEvent.put( relEGColl, "Relaxed" ) ;

   OrphanHandle< L1JetParticleCollection > isoTauHandle =
     iEvent.put( isoTauColl, "IsolatedTau" ) ;

   OrphanHandle< L1JetParticleCollection > relTauHandle =
     iEvent.put( relTauColl, "RelaxedTau" ) ;

   OrphanHandle< L1JetParticleCollection > jetHandle =
     iEvent.put( jetColl, "Jets" ) ;

   OrphanHandle< L1JetParticleCollection > fwdJetHandle =
     iEvent.put( fwdJetColl, "FwdJets" ) ;

   OrphanHandle< L1EtMissParticleCollection > etMissCollHandle =
     iEvent.put( etMissColl, "MET" ) ;

   OrphanHandle< L1EtMissParticleCollection > htMissCollHandle =
     iEvent.put( htMissColl, "MHT" ) ;

}

// ------------ method called once each job just before starting event loop  ------------
void 
UCT2015L1ExtraProducer::beginJob() {

}

// ------------ method called once each job just after ending the event loop  ------------
void 
UCT2015L1ExtraProducer::endJob() {

}

int UCT2015L1ExtraProducer::deltaPhi18(int phi1, int phi2)
{   
    // Compute the difference in phi between two towers, wrapping at phi = 18
    int difference = phi1 - phi2;
    if (std::abs(phi1 - phi2) == 17) {
        difference = -difference/std::abs(difference);
    }
    return difference;
}   


int UCT2015L1ExtraProducer::deltaGctPhi(const L1CaloRegion& r1, const L1CaloRegion& r2)
{
    return deltaPhi18(r1.gctPhi(), r2.gctPhi());
}   


L1CaloRegionCollection::const_iterator UCT2015L1ExtraProducer::matchObjectToRegion(const L1CaloGeometry* geom,
										   const edm::Handle<L1CaloRegionCollection> &regions,
										   const float &eta,
										   const float &phi)
{

  // find closest region to eta value...
  L1CaloRegionCollection::const_iterator matchedRegion = regions->begin();
  float dRMin = 999.9;
  for(L1CaloRegionCollection::const_iterator region = regions->begin();
      region != regions->end(); region++)
    {
      
      double regionEta = geom->etaBinCenter( region->id() ) ;
      double regionPhi = geom->emJetPhiBinCenter( region->id() ) ;
      float dEta = (eta - regionEta);
      float dPhi = (phi - regionPhi);
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      if (dR > dRMin) continue;
      matchedRegion = region;
      dRMin = dR;
    }
  
  return matchedRegion;
  
}


L1GObject UCT2015L1ExtraProducer::buildJetAtIndex(const edm::Handle<L1CaloRegionCollection> &newRegions,
						  const L1CaloRegionCollection::const_iterator &newRegion, 
						  const float &regionLSB_, 
						  float threshold) {

  double regionET = newRegion->et() * regionLSB_;
  if (regionET < threshold) regionET = 0.0;
  
  double neighborN_et = 0;
  double neighborS_et = 0;
  double neighborE_et = 0;
  double neighborW_et = 0;
  double neighborNE_et = 0;
  double neighborSW_et = 0;
  double neighborNW_et = 0;
  double neighborSE_et = 0;
  unsigned int nNeighbors = 0;
  
  //std::cout << "Looking for seed @ " << newRegion->gctPhi() << " " << newRegion->gctEta() << std::endl;
  
  for(L1CaloRegionCollection::const_iterator neighbor = newRegions->begin();
      neighbor != newRegions->end(); neighbor++)
    {
      
      double neighborET = neighbor->et() * regionLSB_;
      if (neighborET < threshold) neighborET = 0.0;
      
      if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	 (newRegion->gctEta()    ) == neighbor->gctEta()) {
	neighborN_et = neighborET;
	nNeighbors++;
	//debug("N", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta()    ) == neighbor->gctEta()) {
	neighborS_et = neighborET;
	nNeighbors++;
	//debug("S", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborE_et = neighborET;
	nNeighbors++;
	//debug("E", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborW_et = neighborET;
	nNeighbors++;
	//debug("W", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborNE_et = neighborET;
	nNeighbors++;
	//debug("NE", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborSW_et = neighborET;
	nNeighbors++;
	//debug("SW", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborNW_et = neighborET;
	nNeighbors++;
	//debug("NW", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborSE_et = neighborET;
	nNeighbors++;
	//debug("SE", *newRegion, *neighbor);
	continue;
      }
    }
  
  //   ---- this should probably be removed for constructing "isolation" jets
  
  //    if(regionET > neighborN_et &&
  //            regionET > neighborNW_et &&
  //            regionET > neighborW_et &&
  //            regionET > neighborSW_et &&
  //            regionET >= neighborNE_et &&
  //            regionET >= neighborE_et &&
  //            regionET >= neighborSE_et &&
  //            regionET >= neighborS_et) 
  //    {

  unsigned int jetET = regionET +
    neighborN_et + neighborS_et + neighborE_et + neighborW_et +
    neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;
  
  // Temporarily use the region granularity -- we will try to improve as above when code is debugged
  int jetPhi = newRegion->gctPhi();
  int jetEta = newRegion->gctEta();
  
  bool neighborCheck = (nNeighbors == 8);
  // On the eta edge we only expect 5 neighbors
  if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
    neighborCheck = true;
  if (!neighborCheck) {
    std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
    assert(false);
  }
  
  L1GObject jet(jetET, jetEta, jetPhi, "MuonSeededJet");
  jet.associatedRegionEt_ = regionET;
  return jet;
  //    }
  
  std::cout << "it's all gone a little bit wrong" << std::endl;
  assert(false);
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(UCT2015L1ExtraProducer);
