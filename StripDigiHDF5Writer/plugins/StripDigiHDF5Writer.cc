// -*- C++ -*-
//
// Package:    CustomModules/StripDigiHDF5Writer
// Class:      StripDigiHDF5Writer
//
/**\class StripDigiHDF5Writer StripDigiHDF5Writer.cc CustomModules/StripDigiHDF5Writer/plugins/StripDigiHDF5Writer.cc

 Description: Writes SiStrip and SiPixel digi data to HDF5.

 Implementation:
     Writes both strip and pixel digis into per-event groups.
*/
//
// Original Author:  Tomas Raila
//         Updated:   ChatGPT, 26 May 2025
//

// system include files
#include <fstream>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

extern "C" {
#include "hdf5.h"
}

class StripDigiHDF5Writer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit StripDigiHDF5Writer(const edm::ParameterSet&);
    ~StripDigiHDF5Writer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    // Tokens
    edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
    edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> digiToken_;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigiToken_;

    // HDF5
    std::string outputFilename_;
    hid_t h5file_;
    hid_t strip_memtype_;
    hid_t pixel_memtype_;

    struct DigiEntry {
        int32_t adc;
        float x, y, z;
    };

    struct PixelEntry {
        int32_t adc;
        int32_t row;
        int32_t col;
        float x, y, z;
    };
};

StripDigiHDF5Writer::StripDigiHDF5Writer(const edm::ParameterSet& iConfig)
    : h5file_(-1), strip_memtype_(-1), pixel_memtype_(-1) {
    digiToken_ = consumes<edm::DetSetVector<SiStripDigi>>(edm::InputTag("simSiStripDigis", "ZeroSuppressed", "DIGI"));
    pixelDigiToken_ = consumes<edm::DetSetVector<PixelDigi>>(edm::InputTag("simSiPixelDigis"));
    tkGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

    outputFilename_ = iConfig.getParameter<std::string>("outputFile");
}

StripDigiHDF5Writer::~StripDigiHDF5Writer() {
    // Nothing here; handled in endJob
}

void StripDigiHDF5Writer::beginJob() {
    h5file_ = H5Fcreate(outputFilename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    strip_memtype_ = H5Tcreate(H5T_COMPOUND, sizeof(DigiEntry));
    H5Tinsert(strip_memtype_, "adc", HOFFSET(DigiEntry, adc), H5T_NATIVE_INT32);
    H5Tinsert(strip_memtype_, "x", HOFFSET(DigiEntry, x), H5T_NATIVE_FLOAT);
    H5Tinsert(strip_memtype_, "y", HOFFSET(DigiEntry, y), H5T_NATIVE_FLOAT);
    H5Tinsert(strip_memtype_, "z", HOFFSET(DigiEntry, z), H5T_NATIVE_FLOAT);

    pixel_memtype_ = H5Tcreate(H5T_COMPOUND, sizeof(PixelEntry));
    H5Tinsert(pixel_memtype_, "adc", HOFFSET(PixelEntry, adc), H5T_NATIVE_INT32);
    H5Tinsert(pixel_memtype_, "row", HOFFSET(PixelEntry, row), H5T_NATIVE_INT32);
    H5Tinsert(pixel_memtype_, "col", HOFFSET(PixelEntry, col), H5T_NATIVE_INT32);
    H5Tinsert(pixel_memtype_, "x", HOFFSET(PixelEntry, x), H5T_NATIVE_FLOAT);
    H5Tinsert(pixel_memtype_, "y", HOFFSET(PixelEntry, y), H5T_NATIVE_FLOAT);
    H5Tinsert(pixel_memtype_, "z", HOFFSET(PixelEntry, z), H5T_NATIVE_FLOAT);
}

void StripDigiHDF5Writer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const auto& tkGeom = iSetup.getData(tkGeomToken_);

    // STRIP DIGIS
    edm::Handle<edm::DetSetVector<SiStripDigi>> digis;
    iEvent.getByToken(digiToken_, digis);

    std::vector<DigiEntry> strip_entries;
    for (const auto& detSet : *digis) {
        uint32_t detid = detSet.id;
        const GeomDet* geomDet = tkGeom.idToDet(detid);
        if (!geomDet) continue;

        const BoundPlane& surface = geomDet->surface();
        for (const auto& digi : detSet.data) {
            int strip = digi.strip();
            int adc = digi.adc();

            LocalPoint local(strip * 0.01, 0.0, 0.0);
            GlobalPoint global = surface.toGlobal(local);

            strip_entries.push_back({adc, global.x(), global.y(), global.z()});
        }
    }

    // PIXEL DIGIS
    edm::Handle<edm::DetSetVector<PixelDigi>> pixelDigis;
    iEvent.getByToken(pixelDigiToken_, pixelDigis);

    std::vector<PixelEntry> pixel_entries;
    for (const auto& detSet : *pixelDigis) {
        uint32_t detid = detSet.id;
        const GeomDet* geomDet = tkGeom.idToDet(detid);
        if (!geomDet) continue;

        const BoundPlane& surface = geomDet->surface();
        for (const auto& digi : detSet.data) {
            int row = digi.row();
            int col = digi.column();
            int adc = digi.adc();

            LocalPoint local(col * 0.01, row * 0.01, 0.0);
            GlobalPoint global = surface.toGlobal(local);

            pixel_entries.push_back({adc, row, col, global.x(), global.y(), global.z()});
        }
    }

    // CREATE HDF5 GROUP FOR EVENT
    std::string groupName = "/event_" + std::to_string(iEvent.id().event());
    hid_t group_id = H5Gcreate2(h5file_, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (!strip_entries.empty()) {
        hsize_t dims[1] = {strip_entries.size()};
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dset_id = H5Dcreate2(group_id, "strip_digis", strip_memtype_, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset_id, strip_memtype_, H5S_ALL, H5S_ALL, H5P_DEFAULT, strip_entries.data());

        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
    }

    if (!pixel_entries.empty()) {
        hsize_t dims[1] = {pixel_entries.size()};
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dset_id = H5Dcreate2(group_id, "pixel_digis", pixel_memtype_, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset_id, pixel_memtype_, H5S_ALL, H5S_ALL, H5P_DEFAULT, pixel_entries.data());

        H5Dclose(dset_id);
        H5Sclose(dataspace_id);
    }

    H5Gclose(group_id);
}

void StripDigiHDF5Writer::endJob() {
    H5Tclose(strip_memtype_);
    H5Tclose(pixel_memtype_);
    H5Fclose(h5file_);
}

void StripDigiHDF5Writer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("outputFile", "strip_pixel_digis.h5");
    descriptions.add("stripDigiHDF5Writer", desc);
}

DEFINE_FWK_MODULE(StripDigiHDF5Writer);
