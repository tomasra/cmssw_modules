import FWCore.ParameterSet.Config as cms
from FWCore.MessageService.MessageLogger_cfi import *

process = cms.Process("DUMPDETIDS")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
# process.GlobalTag.globaltag = autoCond['mc']
process.GlobalTag.globaltag = '140X_mcRun3_2024_realistic_v26'

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# # Enable debug logging
# process.MessageLogger = cms.Service("MessageLogger",
#     destinations = cms.untracked.vstring('cout'),
#     debugModules = cms.untracked.vstring('*'),  # Enable debug for all modules
#     cout = cms.untracked.PSet(
#         threshold = cms.untracked.string('DEBUG'),  # Minimum level to print
#         default = cms.untracked.PSet(limit = cms.untracked.int32(-1))  # No limit
#     )
# )

process.dumpDetIds = cms.EDAnalyzer("DumpDetIds")
process.p = cms.Path(process.dumpDetIds)
