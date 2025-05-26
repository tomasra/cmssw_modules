#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <fstream>

class PixelDigiCSVWriter : public edm::one::EDAnalyzer<> {
public:
  explicit PixelDigiCSVWriter(const edm::ParameterSet&);
  ~PixelDigiCSVWriter() override;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigiToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;

  std::ofstream csvfile_;
};

PixelDigiCSVWriter::PixelDigiCSVWriter(const edm::ParameterSet& iConfig)
  : pixelDigiToken_(consumes<edm::DetSetVector<PixelDigi>>(edm::InputTag("simSiPixelDigis"))),
    trackerGeomToken_(esConsumes())
{
  csvfile_.open("pixel_digis.csv");
  csvfile_ << "event,detId,row,col,adc,global_x,global_y,global_z\n";
}

PixelDigiCSVWriter::~PixelDigiCSVWriter() {}

void PixelDigiCSVWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const auto& pixelDigis = iEvent.get(pixelDigiToken_);
  const auto& trackerGeom = iSetup.getData(trackerGeomToken_);
  const auto eventId = iEvent.id().event();

  for (const auto& detSet : pixelDigis) {
    uint32_t detId = detSet.detId();
    const auto* geomDet = dynamic_cast<const PixelGeomDetUnit*>(trackerGeom.idToDetUnit(DetId(detId)));
    if (!geomDet)
      continue;

    const PixelTopology* topo = &(geomDet->specificTopology());

    for (const auto& digi : detSet) {
      MeasurementPoint mp((float(digi.column()) + 0.5f), (float(digi.row()) + 0.5f));
      LocalPoint lp = topo->localPosition(mp);
      GlobalPoint gp = geomDet->surface().toGlobal(lp);

      csvfile_ << eventId << "," << detId << ","
               << digi.row() << "," << digi.column() << ","
               << digi.adc() << ","
               << gp.x() << "," << gp.y() << "," << gp.z() << "\n";
    }
  }
}

void PixelDigiCSVWriter::endJob()
{
  csvfile_.close();
}

DEFINE_FWK_MODULE(PixelDigiCSVWriter);

