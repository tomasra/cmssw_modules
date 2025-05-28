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


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

extern "C" {
#include "hdf5.h"
}

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

  // For HDF5 tracking
  hsize_t current_size_;
  hsize_t chunk_size_;

  hid_t h5file_;
  hid_t dataset_;
  hid_t dataspace_;
  hid_t memtype_;

  struct DigiEntry {
      int32_t adc;
      float x, y, z;
  };
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
StripDigiHDF5Writer::StripDigiHDF5Writer(const edm::ParameterSet& iConfig)
    : current_size_(0), chunk_size_(1024),
      h5file_(-1), dataset_(-1), dataspace_(-1), memtype_(-1) {
  digiToken_ = consumes<edm::DetSetVector<SiStripDigi>>(
      edm::InputTag("simSiStripDigis", "ZeroSuppressed", "DIGI"));
  tkGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

  outputFilename_ = iConfig.getParameter<std::string>("outputFile");
}

StripDigiHDF5Writer::~StripDigiHDF5Writer() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void StripDigiHDF5Writer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<edm::DetSetVector<SiStripDigi>> digis;
    iEvent.getByToken(digiToken_, digis);

    const auto& tkGeom = iSetup.getData(tkGeomToken_);

    std::vector<DigiEntry> entries;

    for (const auto& detSet : *digis) {
        uint32_t detid = detSet.id;
        const GeomDet* geomDet = tkGeom.idToDet(detid);
        if (!geomDet) {
            edm::LogWarning("StripDigiHDF5Writer") << "No GeomDet found for detid " << detid;
            continue;
        }

        const BoundPlane& surface = geomDet->surface();

        for (const auto& digi : detSet.data) {
            int strip = digi.strip();
            int adc = digi.adc();

            LocalPoint local(strip * 0.01, 0.0, 0.0);
            GlobalPoint global = surface.toGlobal(local);

            DigiEntry entry;
            entry.adc = adc;
            entry.x = global.x();
            entry.y = global.y();
            entry.z = global.z();

            entries.push_back(entry);
        }
    }

    if (!entries.empty()) {
        std::string groupName = "/event_" + std::to_string(iEvent.id().event());
        hid_t group_id = H5Gcreate2(h5file_, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t dims[1] = {entries.size()};
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

        hid_t dset_id = H5Dcreate2(group_id, "digis", memtype_, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset_id, memtype_, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries.data());

        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
    }
}

// ------------ method called once each job just before starting event loop  ------------
void StripDigiHDF5Writer::beginJob() {
    h5file_ = H5Fcreate(outputFilename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    memtype_ = H5Tcreate(H5T_COMPOUND, sizeof(DigiEntry));
    H5Tinsert(memtype_, "adc", HOFFSET(DigiEntry, adc), H5T_NATIVE_INT32);
    H5Tinsert(memtype_, "x", HOFFSET(DigiEntry, x), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype_, "y", HOFFSET(DigiEntry, y), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype_, "z", HOFFSET(DigiEntry, z), H5T_NATIVE_FLOAT);
}

// ------------ method called once each job just after ending the event loop  ------------
void StripDigiHDF5Writer::endJob() {
    H5Tclose(memtype_);
    H5Fclose(h5file_);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void StripDigiHDF5Writer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("outputFile", "strip_digis.h5");
  descriptions.add("stripDigiHDF5Writer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StripDigiHDF5Writer);
