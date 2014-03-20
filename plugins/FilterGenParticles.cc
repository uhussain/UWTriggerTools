#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

class FilterGenParticles : public edm::EDProducer {

public:
  FilterGenParticles (const edm::ParameterSet &);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
private:
  edm::InputTag GenParticleTag_;
  double minPtThreshold_;
  double maxEtaThreshold_;
  int genLevelSelect_;
  double maxIsolation_;
  double isolationCone_;

  double nall;
  double nsel;
};
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <map>
#include <memory>

using namespace edm;
using namespace std;
using namespace reco;


FilterGenParticles::FilterGenParticles( const ParameterSet & cfg ) :
      GenParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles"))),
      minPtThreshold_(cfg.getUntrackedParameter<double> ("MinPtThreshold",5)),
      maxEtaThreshold_(cfg.getUntrackedParameter<double> ("MaxEtaThreshold",5)),
      genLevelSelect_(cfg.getUntrackedParameter<int> ("GenLevelSelect",11)),
      maxIsolation_(cfg.getUntrackedParameter<double> ("MaxIsolation",5)),
      isolationCone_(cfg.getUntrackedParameter<double> ("IsolationCone",0.4))
{
          produces<GenParticleCollection>();
}

void FilterGenParticles::beginJob() {
      nall=0;
      nsel=0;	

}

void FilterGenParticles::endJob() {
     cout<<"********************************************************************"<<endl;
     cout<<"GEN LEVEL FILTERING"<<endl<<endl;
     cout<<"Total Analyzed =   "<<nall<<endl;
     cout<<"GEN Selection  =   "<<nsel<<endl;
     cout<<"********************************************************************"<<endl;
}

void FilterGenParticles::produce (Event & ev, const EventSetup &) {
  nall++;

  bool found=false;

  std::auto_ptr<GenParticleCollection> cleanPart(new GenParticleCollection);


  edm::Handle< vector<reco::GenParticle> >pGenPart;
  if(ev.getByLabel(GenParticleTag_, pGenPart)){

  for( size_t i = 0; i < pGenPart->size(); ++ i ) {
        const reco::GenParticle& genpart = (*pGenPart)[i];
                //cout<<genpart.status()<<" "<<genpart.pt()<<"   "<<fabs(genpart.eta())<<"   "<<genpart.pdgId()<<endl;
        if ( genpart.status()!=1) continue;
        if ( genpart.pt()<minPtThreshold_) continue;
        if ( fabs(genpart.eta())>maxEtaThreshold_) continue;
        if ( fabs(genpart.pdgId())!=genLevelSelect_) continue;

        double sumPT=0;
        for( size_t j = 0; j < pGenPart->size(); ++ j ) {
                const reco::GenParticle& genpart2 = (*pGenPart)[j];
                if (i==j) continue;
                if ( genpart2.status()!=1) continue;
                double DeltaR=deltaR( genpart.eta(), genpart.phi(), genpart2.eta(), genpart2.phi());
                if(DeltaR<isolationCone_) sumPT+=genpart2.pt();
        }
        if(sumPT/genpart.pt()>maxIsolation_) continue;          
        cleanPart->push_back(genpart);
        found=true;
  }
  }

  ev.put(cleanPart); 
  if (found) nsel++;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterGenParticles);
