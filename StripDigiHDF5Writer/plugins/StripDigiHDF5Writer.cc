// -*- C++ -*-
//
// Package:    CustomModules/StripDigiHDF5Writer
// Class:      StripDigiHDF5Writer
//

#include <cstdint>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

extern "C" {
#include "hdf5.h"
}

#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

class StripDigiHDF5Writer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit StripDigiHDF5Writer(const edm::ParameterSet&);
  ~StripDigiHDF5Writer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> stripDigiToken_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigiToken_;

  std::string outputFilename_;
  std::vector<std::string> exportTypes_;
  bool exportStrip_;
  bool exportPixel_;

  hid_t h5file_;
  hid_t memtypeStrip_;
  hid_t memtypePixel_;

  struct StripDigiEntry {
    int32_t adc;
    float x, y, z;
  };

  struct PixelDigiEntry {
    // uint32_t detid;
    int32_t adc;
    // int32_t row, col;
    float x, y, z;
  };
};

StripDigiHDF5Writer::StripDigiHDF5Writer(const edm::ParameterSet& iConfig)
    : h5file_(-1) {

  stripDigiToken_ = consumes<edm::DetSetVector<SiStripDigi>>(
      edm::InputTag("simSiStripDigis", "ZeroSuppressed", "DIGI"));
  pixelDigiToken_ = consumes<edm::DetSetVector<PixelDigi>>(
      edm::InputTag("simSiPixelDigis"));

  tkGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  outputFilename_ = iConfig.getParameter<std::string>("outputFile");

  exportTypes_ = iConfig.getParameter<std::vector<std::string>>("exportTypes");
  exportStrip_ = std::find(exportTypes_.begin(), exportTypes_.end(), "strip") != exportTypes_.end();
  exportPixel_ = std::find(exportTypes_.begin(), exportTypes_.end(), "pixel") != exportTypes_.end();
}

StripDigiHDF5Writer::~StripDigiHDF5Writer() {}

void StripDigiHDF5Writer::beginJob() {
  h5file_ = H5Fcreate(outputFilename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (exportStrip_) {
    memtypeStrip_ = H5Tcreate(H5T_COMPOUND, sizeof(StripDigiEntry));
    H5Tinsert(memtypeStrip_, "adc", HOFFSET(StripDigiEntry, adc), H5T_NATIVE_INT32);
    H5Tinsert(memtypeStrip_, "x", HOFFSET(StripDigiEntry, x), H5T_NATIVE_FLOAT);
    H5Tinsert(memtypeStrip_, "y", HOFFSET(StripDigiEntry, y), H5T_NATIVE_FLOAT);
    H5Tinsert(memtypeStrip_, "z", HOFFSET(StripDigiEntry, z), H5T_NATIVE_FLOAT);
  }

  if (exportPixel_) {
    memtypePixel_ = H5Tcreate(H5T_COMPOUND, sizeof(PixelDigiEntry));
    // H5Tinsert(memtypePixel_, "detid", HOFFSET(PixelDigiEntry, detid), H5T_NATIVE_UINT32);
    H5Tinsert(memtypePixel_, "adc", HOFFSET(PixelDigiEntry, adc), H5T_NATIVE_INT32);
    // H5Tinsert(memtypePixel_, "row", HOFFSET(PixelDigiEntry, row), H5T_NATIVE_INT32);
    // H5Tinsert(memtypePixel_, "col", HOFFSET(PixelDigiEntry, col), H5T_NATIVE_INT32);
    H5Tinsert(memtypePixel_, "x", HOFFSET(PixelDigiEntry, x), H5T_NATIVE_FLOAT);
    H5Tinsert(memtypePixel_, "y", HOFFSET(PixelDigiEntry, y), H5T_NATIVE_FLOAT);
    H5Tinsert(memtypePixel_, "z", HOFFSET(PixelDigiEntry, z), H5T_NATIVE_FLOAT);
  }
}

void StripDigiHDF5Writer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& tkGeom = iSetup.getData(tkGeomToken_);

  std::string groupName = "/event_" + std::to_string(iEvent.id().event());
  hid_t group_id = H5Gcreate2(h5file_, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (exportStrip_) {
    edm::Handle<edm::DetSetVector<SiStripDigi>> stripDigis;
    iEvent.getByToken(stripDigiToken_, stripDigis);

    std::vector<StripDigiEntry> entries;
    for (const auto& detSet : *stripDigis) {
      uint32_t detid = detSet.id;
      const GeomDet* geomDet = tkGeom.idToDet(detid);
      if (!geomDet) continue;

      const BoundPlane& surface = geomDet->surface();
      for (const auto& digi : detSet.data) {
        LocalPoint local(digi.strip() * 0.01, 0.0, 0.0);
        GlobalPoint global = surface.toGlobal(local);

        entries.push_back({digi.adc(), global.x(), global.y(), global.z()});
      }
    }

    if (!entries.empty()) {
      hsize_t dims[1] = {entries.size()};
      hid_t dataspace = H5Screate_simple(1, dims, NULL);
      hid_t dset = H5Dcreate2(group_id, "strip_digis", memtypeStrip_, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dset, memtypeStrip_, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries.data());
      H5Dclose(dset);
      H5Sclose(dataspace);
    }
  }

  if (exportPixel_) {
    edm::Handle<edm::DetSetVector<PixelDigi>> pixelDigis;
    iEvent.getByToken(pixelDigiToken_, pixelDigis);

    std::vector<PixelDigiEntry> entries;
    for (const auto& detSet : *pixelDigis) {
      uint32_t detid = detSet.id;
      const GeomDet* geomDet = tkGeom.idToDet(detid);
      if (!geomDet) continue;

      const BoundPlane& surface = geomDet->surface();
      for (const auto& digi : detSet.data) {
        LocalPoint local(digi.row() * 0.01, digi.column() * 0.01, 0.0);
        GlobalPoint global = surface.toGlobal(local);

        entries.push_back({
          // detid, 
          digi.adc(), 
          // digi.row(), 
          // digi.column(), 
          global.x(), 
          global.y(), 
          global.z()
        });
      }
    }

    if (!entries.empty()) {
      hsize_t dims[1] = {entries.size()};
      hid_t dataspace = H5Screate_simple(1, dims, NULL);
      hid_t dset = H5Dcreate2(group_id, "pixel_digis", memtypePixel_, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dset, memtypePixel_, H5S_ALL, H5S_ALL, H5P_DEFAULT, entries.data());
      H5Dclose(dset);
      H5Sclose(dataspace);
    }
  }

  H5Gclose(group_id);
}

void StripDigiHDF5Writer::endJob() {
  if (exportStrip_) H5Tclose(memtypeStrip_);
  if (exportPixel_) H5Tclose(memtypePixel_);
  H5Fclose(h5file_);
}

void StripDigiHDF5Writer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("outputFile", "strip_pixel_digis.h5");
  desc.add<std::vector<std::string>>("exportTypes", {"strip", "pixel"});
  descriptions.add("stripDigiHDF5Writer", desc);
}

DEFINE_FWK_MODULE(StripDigiHDF5Writer);
