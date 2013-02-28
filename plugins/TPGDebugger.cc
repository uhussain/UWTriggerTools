#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "L1Trigger/UCT2015/interface/helpers.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

class TPGDebugger : public edm::EDFilter {
  public:
    TPGDebugger(const edm::ParameterSet& pset);
    virtual ~TPGDebugger(){}
    bool filter(edm::Event& evt, const edm::EventSetup& es);
  private:
    typedef std::vector<edm::ParameterSet> VPSet;
    edm::InputTag ecalSrc_;
    edm::InputTag hcalSrc_;
    VPSet toPrint_;
    bool filter_;
};

TPGDebugger::TPGDebugger(const edm::ParameterSet& pset) {
  filter_ = pset.exists("filter") ? pset.getParameter<bool>("filter") : false;
  toPrint_ = pset.getParameter<VPSet>("toPrint");
  ecalSrc_ = pset.getParameter<edm::InputTag>("ecalSrc");
  hcalSrc_ = pset.getParameter<edm::InputTag>("hcalSrc");
}

double getEta(int ieta) {
  unsigned int iEta = 0;
  // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
  // And we have mapped [-28, -1], [1, 28] onto [0, 55]
  if(ieta < 0)
    iEta = ieta + 28;
  else if(ieta > 0)
    iEta = ieta + 27;
  return convertTPGEta(iEta);
}

double getPhi(int iPhi) {
  // TPG iPhi starts at 1 and goes to 72.  Let's index starting at zero.
  return convertTPGPhi(iPhi-1);
}

bool TPGDebugger::filter(edm::Event& evt, const edm::EventSetup& es) {

  edm::Handle<EcalTrigPrimDigiCollection> ecal;
  edm::Handle<HcalTrigPrimDigiCollection> hcal;

  const edm::ParameterSet* matchedPSet = NULL;

  for (size_t i = 0; i < toPrint_.size(); ++i) {
    const edm::ParameterSet& pset = toPrint_[i];
    // if run=0 take anything
    if (pset.getParameter<unsigned int>("run") != 0) {
      if (pset.getParameter<unsigned int>("run") != evt.run())
        break;
      if (pset.getParameter<unsigned int>("ls") != evt.id().luminosityBlock())
        break;
      if (pset.getParameter<unsigned int>("event") != evt.id().event())
        break;
    }
    matchedPSet = &pset;
    break;
  }

  // check if we don't care about this event.
  if (!matchedPSet) {
    if (filter_)
      return false;
    return true;
  }

  std::cout << std::endl <<  "[TPG DEBUGGER] run: " << evt.run()
    << " ls: " << evt.id().luminosityBlock() << " evt: "
    << " evt: " << evt.id().event() << std::endl;


  // check if we want to look anywhere
  int minIEta = matchedPSet->getParameter<int>("minIEta");
  int maxIEta = matchedPSet->getParameter<int>("maxIEta");
  int minIPhi = matchedPSet->getParameter<int>("minIPhi");
  int maxIPhi = matchedPSet->getParameter<int>("maxIPhi");

  evt.getByLabel(ecalSrc_, ecal);
  evt.getByLabel(hcalSrc_, hcal);

  std::cout << "ECAL TPGS" << std::endl;
  for (size_t i = 0; i < ecal->size(); ++i) {
    int ieta = (*ecal)[i].id().ieta();
    int iphi = (*ecal)[i].id().iphi();

    if (ieta >= minIEta && ieta <= maxIEta &&
        iphi >= minIPhi && ieta <= maxIPhi) {
      if ((*ecal)[i].compressedEt() > 5) {
        std::cout << "ecal eta/phi=" << ieta << "/" << iphi
          << " = (" << getEta(ieta) << "/" << getPhi(iphi) << ") "
          << " et="<< (*ecal)[i].compressedEt() << " fg=" << (*ecal)[i].fineGrain()
          << std::endl;
      }
    }
  }

  edm::ESHandle<L1CaloHcalScale> hcalScale;
  es.get<L1CaloHcalScaleRcd>().get(hcalScale);

  std::cout << "HCAL TPGS" << std::endl;
  for (size_t i = 0; i < hcal->size(); ++i) {
    int ieta = (*hcal)[i].id().ieta();
    int iphi = (*hcal)[i].id().iphi();
    short absieta = std::abs((*hcal)[i].id().ieta());
    short zside = (*hcal)[i].id().zside();

    if (ieta >= minIEta && ieta <= maxIEta &&
        iphi >= minIPhi && ieta <= maxIPhi) {
      double energy = hcalScale->et(
          (*hcal)[i].SOI_compressedEt(), absieta, zside);
      if (energy > 5) {
        std::cout << "hcal eta/phi=" << ieta << "/" << iphi
          << " = (" << getEta(ieta) << "/" << getPhi(iphi) << ") "
          << " et=" << (*hcal)[i].SOI_compressedEt()
          << " energy=" << energy
          << " fg=" << (*hcal)[i].SOI_fineGrain() << std::endl;
      }
    }
  }
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TPGDebugger);
