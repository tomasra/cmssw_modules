// -*- C++ -*-
//
// Package:    CustomModules/PremixDataFetcher
// Class:      PremixDataFetcher
//
/**\class PremixDataFetcher PremixDataFetcher.cc CustomModules/PremixDataFetcher/plugins/PremixDataFetcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tomas Raila
//         Created:  Fri, 27 Jun 2025 08:02:50 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class PremixDataFetcher : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PremixDataFetcher(const edm::ParameterSet&);
  ~PremixDataFetcher() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  TTree* tree_;
  int event_, detId_, row_, col_, adc_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PremixDataFetcher::PremixDataFetcher(const edm::ParameterSet& iConfig) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  digiToken_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("src"));
}

PremixDataFetcher::~PremixDataFetcher() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void PremixDataFetcher::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  event_ = iEvent.id().event();

  edm::Handle<edm::DetSetVector<PixelDigi>> digis;
  iEvent.getByToken(digiToken_, digis);

  for (const auto& detSet : *digis) {
    detId_ = detSet.detId();
    for (const auto& digi : detSet) {
      row_ = digi.row();
      col_ = digi.column();
      adc_ = digi.adc();
      tree_->Fill();
    }
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void PremixDataFetcher::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("digiTree", "SiPixelDigis");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("detId", &detId_, "detId/I");
  tree_->Branch("row", &row_, "row/I");
  tree_->Branch("col", &col_, "col/I");
  tree_->Branch("adc", &adc_, "adc/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void PremixDataFetcher::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PremixDataFetcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PremixDataFetcher);
