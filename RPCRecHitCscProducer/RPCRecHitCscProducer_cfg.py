import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

#process = cms.Process("testRPCDigiMerger")
process = cms.Process("REClUSTERIZATION")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"

"""#######################################################
### RPC RawToDigi
### RPC RawToDigi - from Legacy
process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

### RPC RawToDigi - from TwinMux
process.load("EventFilter.RPCRawToDigi.RPCTwinMuxRawToDigi_cff")

### RPC RawToDigi - from CPPF
process.load("EventFilter.RPCRawToDigi.RPCCPPFRawToDigi_cff")
# process.load("EventFilter.RPCRawToDigi.RPCCPPFRawToDigi_sqlite_cff") #to load CPPF link maps from the local DB

### RPC RawToDigi - from OMTF
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.omtfStage2Digis = cms.EDProducer("OmtfUnpacker",
  inputLabel = cms.InputTag('rawDataCollector'),
)

process.load("EventFilter.RPCRawToDigi.RPCDigiMerger_cff")
"""
########## RPC RecHit #############
#process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")
###################################

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Source
process.source = cms.Source("PoolSource"
                            , fileNames = cms.untracked.vstring(
#                                    'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/data/FEAF370D-6243-E811-9145-A0369FD0B122.root'
#                                    'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/data/FE9BC2AB-D63B-E811-A038-0CC47A4DEEF8.root'
#                                    '/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/FE54D6E5-0A42-E811-8FD4-48D539F3863E.root'
#                                    '/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/FC400A2C-8C43-E811-AC8D-48D539D3335B.root'
                                    '/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/F8DEC9CC-3A3B-E811-A81A-0CC47A4DED92.root'
                            )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

"""
process.p = cms.Path( 
                       (
                        process.rpcUnpackingModule 
                        + process.rpcTwinMuxRawToDigi 
                        + process.rpcCPPFRawToDigi 
                        + process.omtfStage2Digis
                        )
                       * process.rpcDigiMerger 
)
"""

process.p = cms.Path( process.rpcRecHits )

# Output
process.out = cms.OutputModule("PoolOutputModule"
                               , outputCommands = cms.untracked.vstring("drop *"
#                                                                        , "keep *_*_*_testRPCDigiMerger")
                                                                        , "keep *_simMuonRPCDigis__*"
                                                                        , "keep *_simCscTriggerPrimitiveDigis_MPCSORTED_*")
                               # , fileName = cms.untracked.string(options.outputFile)
                               , fileName = cms.untracked.string("testRPCoutput.root")
                               #, fileName = cms.untracked.string("testRPCDigiMerger.root")
                               , SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)


process.e = cms.EndPath(process.out)
