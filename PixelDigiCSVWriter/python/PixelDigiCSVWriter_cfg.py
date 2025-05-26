import FWCore.ParameterSet.Config as cms

process = cms.Process("DUMPDIGIS")

# Load services
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Geometry
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "140X_mcRun3_2024_realistic_v26"

# Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/mc/RunIIISummer24PrePremix/Neutrino_E-10_gun/PREMIX/Premixlib2024_140X_mcRun3_2024_realistic_v26-v1/140007/6b22f430-90d3-4e66-99f1-7fccdae45317.root'
    ),
    inputCommands = cms.untracked.vstring('keep *')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Analyzer
process.trackerDigiCSVWriter = cms.EDAnalyzer("PixelDigiCSVWriter")

# Path
process.p = cms.Path(process.trackerDigiCSVWriter)
