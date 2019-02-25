import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

process = cms.Process("REClUSTERIZATION")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"

#RPCRecHit
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits
process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Source
process.source = cms.Source("PoolSource"
                            , fileNames = cms.untracked.vstring(
#                                    'root://xrootd-cms.infn.it///store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/FEAF370D-6243-E811-9145-A0369FD0B122.root'
                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/F8DEC9CC-3A3B-E811-A81A-0CC47A4DED92.root'
                            )
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.p = cms.Path( process.rpcRecHits )

# Output
process.out = cms.OutputModule("PoolOutputModule"
                               , outputCommands = cms.untracked.vstring("drop *"
                                                                        , "keep *_simMuonRPCDigis__*"
                                                                        , "keep *_rpcRecHits__*"
                                                                        , "keep *_simCscTriggerPrimitiveDigis_MPCSORTED_*")
                               , fileName = cms.untracked.string("RPCoutput_withRecHitCSC.root")
                               , SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)


process.e = cms.EndPath(process.out)
