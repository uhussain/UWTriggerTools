#include "L1Trigger/UCT2015/plugins/UCT2015GctCandsProducer.h"

// system includes
#include <memory>
#include <vector>

// EDM includes
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Trigger configuration includes
#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/L1TObjects/interface/L1GctJetFinderParams.h"
#include "CondFormats/L1TObjects/interface/L1GctChannelMask.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1HtMissScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1HfRingEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1GctJetFinderParamsRcd.h"
#include "CondFormats/DataRecord/interface/L1GctChannelMaskRcd.h"

// GCT include files
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetEtCalibrationLut.h"

// RCT data includes
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "L1Trigger/UCT2015/interface/UCTCandidate.h"

// GCT data includes
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"


// Scales
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1HtMissScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1HfRingEtScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1GctJetFinderParams.h"
#include "CondFormats/DataRecord/interface/L1GctJetFinderParamsRcd.h"

using std::vector;

UCT2015GctCandsProducer::UCT2015GctCandsProducer(const edm::ParameterSet& ps) :
  egSourceRlx_(ps.getParameter<edm::InputTag>("egRelaxed")),
  egSourceIso_(ps.getParameter<edm::InputTag>("egIsolated")),
  tauSourceRlx_(ps.getParameter<edm::InputTag>("tauRelaxed")),
  tauSourceIso_(ps.getParameter<edm::InputTag>("tauIsolated")),
  jetSource_(ps.getParameter<edm::InputTag>("jetSource")),
  setSource_(ps.getParameter<edm::InputTag>("setSource")),
  shtSource_(ps.getParameter<edm::InputTag>("shtSource")),
  metSource_(ps.getParameter<edm::InputTag>("metSource")),
  mhtSource_(ps.getParameter<edm::InputTag>("mhtSource")),
  maxEGs_(ps.getUntrackedParameter<int>("maxEGs",4)),
  maxIsoEGs_(ps.getUntrackedParameter<int>("maxIsoEGs",4)),
  maxTaus_(ps.getUntrackedParameter<int>("maxTaus",4)),
  maxIsoTaus_(ps.getUntrackedParameter<int>("maxIsoTaus",4)),
  maxJets_(ps.getUntrackedParameter<int>("maxJets",4))
 {

  // list of products
  produces<L1GctEmCandCollection>("isoEm");
  produces<L1GctEmCandCollection>("rlxEm");
  produces<L1GctJetCandCollection>("isoTau");
//  produces<L1GctJetCandCollection>("rlxTau");

  produces<L1GctJetCandCollection>("cenJets");
  produces<L1GctJetCandCollection>("forJets");
 
  produces<L1GctEtTotalCollection>();
  produces<L1GctEtHadCollection>();

  produces<L1GctEtMissCollection>();
  produces<L1GctHtMissCollection>();

  produces<L1GctHFBitCountsCollection>();
  produces<L1GctHFRingEtSumsCollection>();
}

UCT2015GctCandsProducer::~UCT2015GctCandsProducer() {
}


void UCT2015GctCandsProducer::beginJob()
{
}

void UCT2015GctCandsProducer::endJob()
{
}

void UCT2015GctCandsProducer::produce(edm::Event& e, const edm::EventSetup& c) {

  // The emulator will always produce output collections, which get filled as long as
  // the setup and input data are present. Start by making empty output collections.

  // create the em and tau collections
  std::auto_ptr<L1GctEmCandCollection> isoEmResult   (new L1GctEmCandCollection( ) );
  std::auto_ptr<L1GctEmCandCollection> rlxEmResult(new L1GctEmCandCollection( ) );
  std::auto_ptr<L1GctJetCandCollection> isoTauResult   (new L1GctJetCandCollection( ) );
  std::auto_ptr<L1GctJetCandCollection> rlxTauResult(new L1GctJetCandCollection( ) );

  // create the jet collections
  std::auto_ptr<L1GctJetCandCollection> cenJetResult(new L1GctJetCandCollection( ) );
  std::auto_ptr<L1GctJetCandCollection> forJetResult(new L1GctJetCandCollection( ) );

  // create the energy sum digis
  std::auto_ptr<L1GctEtTotalCollection> etTotResult (new L1GctEtTotalCollection( ) );
  std::auto_ptr<L1GctEtHadCollection>   etHadResult (new L1GctEtHadCollection  ( ) );
  std::auto_ptr<L1GctEtMissCollection>  etMissResult(new L1GctEtMissCollection ( ) );
  std::auto_ptr<L1GctHtMissCollection>  htMissResult(new L1GctHtMissCollection ( ) );

  // create the Hf sums digis
  std::auto_ptr<L1GctHFBitCountsCollection>  hfBitCountResult (new L1GctHFBitCountsCollection ( ) );
  std::auto_ptr<L1GctHFRingEtSumsCollection> hfRingEtSumResult(new L1GctHFRingEtSumsCollection( ) );

 
  // Getting the scales 

   edm::ESHandle< L1CaloEtScale > emScale ;
   c.get< L1EmEtScaleRcd >().get( emScale ) ;

   edm::ESHandle< L1CaloEtScale > jetScale ;
   c.get< L1JetEtScaleRcd >().get( jetScale ) ;

   edm::ESHandle< L1CaloEtScale > hwForJetScale ;
   c.get< L1JetEtScaleRcd >().get( hwForJetScale ) ;

   edm::ESHandle< L1GctJetFinderParams > jetFinderParams ;
   c.get< L1GctJetFinderParamsRcd >().get( jetFinderParams ) ;

   edm::ESHandle< L1CaloEtScale > htMissScale ;
   std::vector< bool > htMissMatched ;
   c.get< L1HtMissScaleRcd >().get( htMissScale ) ;

   double etSumLSB = jetScale->linearLsb() ;
   double htSumLSB = jetFinderParams->getHtLsbGeV();
   
  // And here we go! fro UCT objects to something that the GT will understand

     // EG
      edm::Handle< UCTCandidateCollection > egObjs ;
      e.getByLabel( egSourceRlx_, egObjs ) ;

      if( !egObjs.isValid() ) {
                LogDebug("")<<"EG Collection not found - check name";
      }
      else {
                for( unsigned int i = 0 ; i<egObjs->size() && i<maxEGs_; i++){
                        UCTCandidate itr=egObjs->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");
                        unsigned gctEta=((rctEta & 0x7) | (iEta<11 ? 0x8 : 0x0));
                        unsigned rank = emScale->rank( itr.pt() ) ;
                        L1GctEmCand gctEmCand=L1GctEmCand(rank,iPhi,gctEta,0);        
                        rlxEmResult->push_back( gctEmCand  );
                }
                if(rlxEmResult->size()<maxEGs_)
                        for ( unsigned int j = 0 ; j<(maxEGs_-egObjs->size()); j++){
                                L1GctEmCand gctEmCand=L1GctEmCand( 0,(unsigned)0,(unsigned)0,0);
                                isoEmResult->push_back( gctEmCand  );
                        }

      }


     // isoEG
      edm::Handle< UCTCandidateCollection > egObjsIso ;
      e.getByLabel( egSourceIso_, egObjsIso ) ;

      if( !egObjsIso.isValid() ) {
                LogDebug("")<<"isoEG Collection not found - check name";
      }
      else {
                for( unsigned int i = 0 ; i<egObjsIso->size() && i<maxIsoEGs_; i++){
                        UCTCandidate itr=egObjsIso->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");

                        unsigned gctEta=((rctEta & 0x7) | (iEta<11 ? 0x8 : 0x0));

                        unsigned rank = emScale->rank( itr.pt() ) ;

                        L1GctEmCand gctEmCand=L1GctEmCand(rank,iPhi,gctEta,1);
                        isoEmResult->push_back( gctEmCand  );
                }
                if(isoEmResult->size()<maxIsoEGs_)
                        for ( unsigned int j = 0 ; j<(maxIsoEGs_-egObjsIso->size()); j++){
                                L1GctEmCand gctEmCand=L1GctEmCand( 0,(unsigned)0,(unsigned)0,1);
                                isoEmResult->push_back( gctEmCand  );
                        }
      }


     // TAU 
     // They need to be treated as Jets! 
/*   
      edm::Handle< UCTCandidateCollection > tauObjs ;
      e.getByLabel( tauSourceRlx_, tauObjs ) ;

      if( !tauObjs.isValid() ) {
                LogDebug("")<<"TAU Collection not found - check name";
      }
      else {
                for( unsigned int i = 0 ; i<tauObjs->size() && i<maxTaus_; i++){
                        UCTCandidate itr=tauObjs->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");
                        unsigned hwEta=(((rctEta % 7) & 0x7) | (iEta<11 ? 0x8 : 0));
                        unsigned hwPhi= iPhi& 0x1f;
                        const int16_t bx=0; 
                        double pt=itr.getFloat("associatedRegionEt");
                        unsigned rank = jetScale->rank(pt);
                        bool isFor=false;
                        bool isTau=true;
                        L1GctJetCand gctJetCand=L1GctJetCand(rank, hwPhi, hwEta, isTau , isFor,(uint16_t) 0, (uint16_t) 0, bx);
                        rlxTauResult->push_back( gctJetCand  );
                }
                if(rlxTauResult->size()<maxTaus_)
                        for ( unsigned int j = 0 ; j<(maxTaus_-tauObjs->size()); j++){
                                L1GctJetCand gctJetCand=L1GctJetCand( 0,(unsigned)0,(unsigned)0,1,0,(uint16_t) 0, (uint16_t) 0, 0);
                                isoTauResult->push_back( gctJetCand  );
                        }
      }
*/

     // isoTAU
     // We'll only use this one    
      edm::Handle< UCTCandidateCollection > tauObjsIso ;
      e.getByLabel( tauSourceIso_, tauObjsIso ) ;

      if( !tauObjsIso.isValid() ) {
                LogDebug("")<<"isoTAU Collection not found - check name";
      }
      else {
                for( unsigned int i = 0 ; i<tauObjsIso->size() && i<maxIsoTaus_; i++){
                        UCTCandidate itr=tauObjsIso->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");
                        unsigned hwEta=(((rctEta % 7) & 0x7) | (iEta<11 ? 0x8 : 0));
                        unsigned hwPhi= iPhi& 0x1f;
                        const int16_t bx=0; 
                        double pt=itr.getFloat("associatedRegionEt");
                        unsigned rank = jetScale->rank(pt);
                        bool isFor=false;
                        bool isTau=true;
                        L1GctJetCand gctJetCand=L1GctJetCand(rank, hwPhi, hwEta, isTau , isFor,(uint16_t) 0, (uint16_t) 0, bx);
                        isoTauResult->push_back( gctJetCand  );
                }
                if(isoTauResult->size()<maxIsoTaus_)
                        for ( unsigned int j = 0 ; j<(maxIsoTaus_-tauObjsIso->size()); j++){
                                L1GctJetCand gctJetCand=L1GctJetCand( 0,(unsigned)0,(unsigned)0,1,0,(uint16_t) 0, (uint16_t) 0, 0);
                                isoTauResult->push_back( gctJetCand  );
                        }

      }


    // Jets
      edm::Handle< UCTCandidateCollection > jetObjs ;
      e.getByLabel( jetSource_, jetObjs) ;

      if( !jetObjs.isValid() ) {
                LogDebug("")<<"JET Collection not found - check name";
      }
      else {
                for( unsigned int i = 0 ; i<jetObjs->size() &&  i<maxJets_; i++){
                        UCTCandidate itr=jetObjs->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");
                        unsigned hwEta=(((rctEta % 7) & 0x7) | (iEta<11 ? 0x8 : 0));
                        unsigned hwPhi= iPhi& 0x1f;
                        bool isTau=false;
                        const int16_t bx=0; 
                        unsigned rank = jetScale->rank(itr.pt()) ;

                        bool isFor=(rctEta>=7);
                        L1GctJetCand gctJetCand=L1GctJetCand(rank, hwPhi, hwEta, isTau , isFor,(uint16_t) 0, (uint16_t) 0, bx);
                        if (!isFor) cenJetResult->push_back( gctJetCand  );
                }
                if(cenJetResult->size()<maxJets_)
                        for ( unsigned int j = 0 ; j<(maxJets_-cenJetResult->size()); j++){
                                L1GctJetCand gctJetCand=L1GctJetCand(0,0,0,0,0,(uint16_t) 0, (uint16_t) 0,0);
                                cenJetResult->push_back( gctJetCand  );
                        }
                for( unsigned int i = 0 ; i<maxJets_ && i<jetObjs->size(); i++){
                        UCTCandidate itr=jetObjs->at(i);
                        unsigned iEta=itr.getInt("rgnEta");
                        unsigned iPhi=itr.getInt("rgnPhi");
                        unsigned rctEta=itr.getInt("rctEta");
                        unsigned hwEta=(((rctEta % 7) & 0x7) | (iEta<11 ? 0x8 : 0));
                        unsigned hwPhi= iPhi& 0x1f;
                        bool isTau=false;
                        const int16_t bx=0; 
                        unsigned rank = jetScale->rank(itr.pt()) ;

                        bool isFor=(rctEta>=7);
                        L1GctJetCand gctJetCand=L1GctJetCand(rank, hwPhi, hwEta, isTau , isFor,(uint16_t) 0, (uint16_t) 0, bx);
                        if (isFor) forJetResult->push_back( gctJetCand  );
                }
                if(forJetResult->size()<maxJets_)
                        for ( unsigned int j = 0 ; j<(maxJets_-forJetResult->size()); j++){
                                L1GctJetCand gctJetCand=L1GctJetCand(0,0,0,0,1,(uint16_t) 0, (uint16_t) 0,0);
                                forJetResult->push_back( gctJetCand  );
                        }

      }


    // Sums
    // SumET:

      edm::Handle< UCTCandidateCollection > setObjs ;
      e.getByLabel( setSource_, setObjs) ;

      if( !setObjs.isValid() ) {
                LogDebug("")<<"SET Collection not found - check name";
      }
      else {
                if(setObjs->size()>0){ // This is just for safety        
                        UCTCandidate itr=setObjs->at(0);
                        const int16_t bx=0; // ???
                        double convert=itr.pt()/etSumLSB;
                        unsigned rank=(unsigned)convert;
                        L1GctEtTotal gctSumEt=L1GctEtTotal(rank, 0, bx);
                        etTotResult->push_back(gctSumEt);        

                        }

      }

    // SumHT:

      edm::Handle< UCTCandidateCollection > shtObjs ;
      e.getByLabel( shtSource_, shtObjs) ;

      if( !shtObjs.isValid() ) {
                LogDebug("")<<"SHT Collection not found - check name";
      }
      else {
                if(shtObjs->size()>0){ // This is just for safety
                        UCTCandidate itr=shtObjs->at(0);
                        double convert=itr.pt()/htSumLSB;
                        unsigned rank=(unsigned)convert;
                        const int16_t bx=0; // ???
                        L1GctEtHad gctSumHt=L1GctEtHad(rank, 0, bx);
                        etHadResult->push_back(gctSumHt);

                        }

      }


    // MET:

      edm::Handle< UCTCandidateCollection > metObjs ;
      e.getByLabel( metSource_, metObjs) ;

      if( !metObjs.isValid() ) {
                LogDebug("")<<"MET Collection not found - check name";
      }
      else {
                if(metObjs->size()>0){ // This is just for safety
                        UCTCandidate itr=metObjs->at(0);
                        double phiMod=36.*itr.phi()/M_PI;
                        if(phiMod<0) phiMod  += 72;
                        unsigned iPhi = (unsigned)phiMod;
                        double convert=itr.pt()/etSumLSB;
                        unsigned rank=(unsigned)convert;
                        L1GctEtMiss gctMET=L1GctEtMiss(rank, iPhi, 0);  
                        etMissResult->push_back(gctMET);
                        }
  
      }


    // MHT:

      edm::Handle< UCTCandidateCollection > mhtObjs ;
      e.getByLabel( mhtSource_, mhtObjs) ;

      if( !mhtObjs.isValid() ) {
                LogDebug("")<<"MHT Collection not found - check name";
      }
      else {
                if(mhtObjs->size()>0){ // This is just for safety
                        UCTCandidate itr=mhtObjs->at(0);
                        double phiMod=9.*itr.phi()/M_PI;
                        if(phiMod<0) phiMod  += 18.0;
                        unsigned iPhi = (unsigned)phiMod;
                        unsigned rank=htMissScale->rank(itr.pt());
                        L1GctHtMiss gctMHT=L1GctHtMiss(rank, iPhi, 0);  
                        htMissResult->push_back(gctMHT);
                        }

      }




  // put the collections into the event
  e.put(rlxEmResult,"rlxEm");
  e.put(isoEmResult,"isoEm");

  e.put(isoTauResult,"isoTau");
//  e.put(rlxTauResult,"rlxTau");

  e.put(cenJetResult,"cenJets");
  e.put(forJetResult,"forJets");

  e.put(etTotResult);
  e.put(etHadResult);
  e.put(etMissResult);
  e.put(htMissResult);

  e.put(hfBitCountResult); // Empty for now
  e.put(hfRingEtSumResult);


}

DEFINE_FWK_MODULE(UCT2015GctCandsProducer);
