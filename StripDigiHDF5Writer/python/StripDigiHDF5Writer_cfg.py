import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Setup command line arguments safely
options = VarParsing.VarParsing('analysis')

options.setDefault('inputFiles', [])
options.setDefault('outputFile', 'out.h5')
options.setDefault('maxEvents', 2)

options.parseArguments()

# Define the process
process = cms.Process("STRIP2HDF5")

# Load geometry from the conditions DB
# process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
# process.load("Configuration.Geometry.GeometryDD4hep_cff")
# process.load("Configuration.Geometry.GeometryDD4hepReco_cff")

# Optional tracker parameters, if needed
# process.load("Geometry.TrackerGeometryBuilder.trackerParameters_cfi")

# Magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2023_realistic', '')

# Set maxEvents from the command line
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles)
)

# Configure your analyzer with output file name
process.analyzer = cms.EDAnalyzer("StripDigiHDF5Writer",
    outputFile = cms.string(options.outputFile)
)

process.p = cms.Path(process.analyzer)
