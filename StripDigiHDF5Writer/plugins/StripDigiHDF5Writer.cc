// -*- C++ -*-
//
// Package:    CustomModules/StripDigiHDF5Writer
// Class:      StripDigiHDF5Writer
//
/**\class StripDigiHDF5Writer StripDigiHDF5Writer.cc CustomModules/StripDigiHDF5Writer/plugins/StripDigiHDF5Writer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tomas Raila
//         Created:  Mon, 26 May 2025 13:27:05 GMT
//
//

// system include files
#include <cstdint>
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class StripDigiHDF5Writer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit StripDigiHDF5Writer(const edm::ParameterSet&);
  ~StripDigiHDF5Writer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> digiToken_;
  std::ofstream outFile_;
  std::string outputFilename_;

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
StripDigiHDF5Writer::StripDigiHDF5Writer(const edm::ParameterSet& iConfig) {
  digiToken_ = consumes<edm::DetSetVector<SiStripDigi>>(
      edm::InputTag("simSiStripDigis", "ZeroSuppressed", "DIGI"));
  tkGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

  outputFilename_ = iConfig.getParameter<std::string>("outputFile");
  outFile_.open(outputFilename_);
  outFile_ << "event,detid,strip,adc,x,y,z\n";
}

StripDigiHDF5Writer::~StripDigiHDF5Writer() {
  outFile_.close();
}

//
// member functions
//

// ------------ method called for each event  ------------
void StripDigiHDF5Writer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get SiStripDigis
  edm::Handle<edm::DetSetVector<SiStripDigi>> digis;
  iEvent.getByToken(digiToken_, digis);

  // Get TrackerGeometry from EventSetup
  const auto& tkGeom = iSetup.getData(tkGeomToken_);

  for (const auto& detSet : *digis) {
    uint32_t detid = detSet.id;

    // Get the detector unit
    const GeomDet* geomDet = tkGeom.idToDet(detid);
    if (!geomDet) {
      edm::LogWarning("SiStripDigiToCSVAnalyzer") << "No GeomDet found for detid " << detid;
      continue;
    }

    const BoundPlane& surface = geomDet->surface();

    for (const auto& digi : detSet.data) {
      int strip = digi.strip();
      int adc = digi.adc();

      // Calculate local strip coordinate (center of the strip)
      LocalPoint local(strip * 0.01, 0.0, 0.0);  // 0.01 cm pitch approximation

      // Transform to global
      GlobalPoint global = surface.toGlobal(local);

      // Write to CSV
      outFile_ << iEvent.id().event() << ","
               << detid << ","
               << strip << ","
               << adc << ","
               << global.x() << ","
               << global.y() << ","
               << global.z() << "\n";

      /*
      The CMS global coordinate system is defined as:
      z axis → along the beam line
      x axis → pointing towards the center of the LHC ring (from CMS)
      y axis → upwards (perpendicular to Earth's surface)

      Units are cm.
      
      */
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void StripDigiHDF5Writer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void StripDigiHDF5Writer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void StripDigiHDF5Writer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(StripDigiHDF5Writer);
