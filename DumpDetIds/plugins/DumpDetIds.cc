// -*- C++ -*-
//
// Package:    CustomModules/DumpDetIds
// Class:      DumpDetIds
//
/**\class DumpDetIds DumpDetIds.cc CustomModules/DumpDetIds/plugins/DumpDetIds.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tomas Raila
//         Created:  Thu, 19 Jun 2025 10:00:24 GMT
//
//

// system include files
// #include <cstddef>
#include <memory>
#include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <unistd.h>

// user include files
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class DumpDetIds : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DumpDetIds(const edm::ParameterSet&);
  ~DumpDetIds() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void writeLocalCoordMap(const edm::EventSetup& iSetup, const std::string& filename);
  void writePixelDetJsonFragment(const PixelGeomDetUnit* pixelDet, 
                                 edm::ESHandle<TrackerTopology> trackerTopo,
                                 std::ofstream& out);
  // ----------member data ---------------------------
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  
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
DumpDetIds::DumpDetIds(const edm::ParameterSet& iConfig) {
    // : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  tkGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  tTopoToken_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
}

DumpDetIds::~DumpDetIds() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

GlobalPoint rowcolToGlobal_builtin(const PixelGeomDetUnit* pixelDet,
                                         const PixelTopology* pixelTopo,
                                         int row,
                                         int col) {
  MeasurementPoint mp(row, col);
  LocalPoint lp = pixelTopo->localPosition(mp);
  GlobalPoint gp = pixelDet->surface().toGlobal(lp);
  return gp;
}


void DumpDetIds::writeLocalCoordMap(const edm::EventSetup& iSetup, const std::string& filename) {
  edm::ESHandle<TrackerGeometry> tracker = iSetup.getHandle(tkGeomToken_);

  // Use the first BPIX module to get the topology
  const PixelGeomDetUnit* det = dynamic_cast<const PixelGeomDetUnit*>(tracker->detsPXB().front());
  if (!det) {
    edm::LogError("DumpDetIds") << "Failed to cast first det to PixelGeomDetUnit";
    return;
  }

  const PixelTopology* topo = &(det->specificTopology());

  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    edm::LogError("DumpDetIds") << "Failed to open file: " << filename;
    return;
  }

  // Write CSV header
  outFile << "row,col,local_x,local_y\n";

  for (int row = 0; row < topo->nrows(); row++) {
    for (int col = 0; col < topo->ncolumns(); col++) {
      MeasurementPoint mp(row, col);
      LocalPoint lp = topo->localPosition(mp);
      outFile << row << "," << col << "," << lp.x() << "," << lp.y() << "\n";
    }
  }

  outFile.close();
  edm::LogInfo("DumpDetIds") << "Local coordinate map written to: " << filename;
}

void DumpDetIds::writePixelDetJsonFragment(const PixelGeomDetUnit* pixelDet, 
                                           edm::ESHandle<TrackerTopology> trackerTopo,
                                           std::ofstream& out) {
  const DetId detId = pixelDet->geographicalId();
  const GlobalPoint position = pixelDet->position();
  const Surface& surface = pixelDet->surface();
  const Surface::RotationType& rot = surface.rotation();

  out << "  {\n";
  out << "    \"det_id\": " << detId.rawId() << ",\n";

  // Based on: https://github.com/CMSTrackerDPG/SiPixelTools-PixelTrees/blob/master/plugins/DetectorInformation.cc
  if (detId.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {
    PixelBarrelName barrelName(detId);
    out << "    \"subdet\": \"BPIX\",\n";
    out << "    \"layer\": "  << trackerTopo->pxbLayer(detId)  << ",\n";
    out << "    \"ladder\": " << trackerTopo->pxbLadder(detId) << ",\n";
    out << "    \"module\": " << trackerTopo->pxbModule(detId) << ",\n";
  } else if (detId.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
    PixelEndcapName endcapName(detId);
    out << "    \"subdet\": \"FPIX\",\n";
    out << "    \"disk\": "   << trackerTopo->pxfDisk(detId)  << ",\n";
    out << "    \"blade\": "  << trackerTopo->pxfBlade(detId)  << ",\n";
    out << "    \"panel\": "  << trackerTopo->pxfPanel(detId)  << ",\n";
    out << "    \"module\": " << trackerTopo->pxfModule(detId) << ",\n";
    out << "    \"side\": "   << trackerTopo->pxfSide(detId) << ",\n";
  }

  out << "    \"position\": [" << position.x() << ", " << position.y() << ", " << position.z() << "],\n";
  out << "    \"rotation\": [\n";
  out << "      [" << rot.xx() << ", " << rot.xy() << ", " << rot.xz() << "],\n";
  out << "      [" << rot.yx() << ", " << rot.yy() << ", " << rot.yz() << "],\n";
  out << "      [" << rot.zx() << ", " << rot.zy() << ", " << rot.zz() << "]\n";
  out << "    ]\n";
  out << "  }";
}


// ------------ method called for each event  ------------
void DumpDetIds::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // std::cout << "DumpDetIds: analyze called for event " << iEvent.id() << std::endl;
  using namespace edm;

  edm::ESHandle<TrackerGeometry> tracker = iSetup.getHandle(tkGeomToken_);
  edm::ESHandle<TrackerTopology> tTopo = iSetup.getHandle(tTopoToken_);
  writeLocalCoordMap(iSetup, "rowcol_to_local.csv");

  std::ofstream outBpix("detids_bpix.json");
  outBpix << std::fixed << std::setprecision(6);
  outBpix << "[\n";

  // BPIX
  const auto& detsBpix = tracker->detsPXB();
  for (size_t i = 0; i < detsBpix.size(); ++i) {
    const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(detsBpix[i]);
    if (!pixelDet) continue;
    
    writePixelDetJsonFragment(pixelDet, tTopo, outBpix);
    if (i != detsBpix.size() - 1) {
      outBpix << ",\n";
    } else {
      outBpix << "\n";
    }

    /*
    // Row/col -> global coordinates
    for (int row = 0; row <= 2; row++) {
      for (int col = 0; col <= 2; col++) {
        // Built-in method
        MeasurementPoint mp(row, col);
        LocalPoint lp = pixelTopo->localPosition(mp);
        Plane surface = pixelDet->surface();
        GlobalPoint gp1 = rowcolToGlobal_builtin(pixelDet, pixelTopo, row, col);

        std::cout << "Row: " << row << ", Col: " << col 
                  << " Local Position: (" << lp.x() << ", " << lp.y() << ", " << ")"
                  << " Global Position (builtin): (" << gp1.x() << ", " << gp1.y() << ", " << gp1.z() << ")" 
                  << std::endl;
      }
    }
    */
  }
  outBpix << "]\n";
  outBpix.close();

  // FPIX
  std::ofstream outFpix("detids_fpix.json");
  outFpix << std::fixed << std::setprecision(6);
  outFpix << "[\n";

  const auto& detsFpix = tracker->detsPXF();
  for (size_t i = 0; i < detsFpix.size(); ++i) {
    const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(detsFpix[i]);
    if (!pixelDet) continue;

    writePixelDetJsonFragment(pixelDet, tTopo, outFpix);
    if (i != detsFpix.size() - 1) {
      outFpix << ",\n";
    } else {
      outFpix << "\n";
    }
  }
  outFpix << "]\n";
  outFpix.close();


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void DumpDetIds::beginJob() {
  // please remove this method if not needed
  // std::cout << "DumpDetIds: beginJob called" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void DumpDetIds::endJob() {
  // please remove this method if not needed
  // std::cout << "DumpDetIds: endJob called" << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DumpDetIds::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(DumpDetIds);
