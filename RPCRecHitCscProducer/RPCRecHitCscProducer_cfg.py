import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

process = cms.Process("REClUSTERIZATION")

process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#for 93X sample
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D38Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D38_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#RPCRecHit
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits
process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Source
process.source = cms.Source("PoolSource"
                            , fileNames = cms.untracked.vstring(
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/F8DEC9CC-3A3B-E811-A81A-0CC47A4DED92.root'
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/30000/FEA28F45-58A6-E811-A000-0242AC130002.root'
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/30000/FA0DA254-10A7-E811-A7A2-1CB72C0A3DC5.root'
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/30000/F883D2F3-84A8-E811-848A-AC1F6B8DD1F8.root'
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/30000/F4B778FD-23A7-E811-92B4-509A4C74913F.root'
#                                    'root://cms-xrd-global.cern.ch//store/mc/PhaseIIFall17D/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_14TeV-pythia8/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/30000/F01D7F0E-A6A6-E811-A4E5-B49691091EA2.root'
#                                    'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/RPCRecHitCscProducer/02D5AE85-0DA7-E811-9C73-44A84225C911.root'
                                    'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_10_5_0_pre1/src/RPC-CSCTrigger/RPCRecHitCscProducer/FEAF370D-6243-E811-9145-A0369FD0B122.root'
                            )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.p = cms.Path( process.rpcRecHits )

# Output
process.out = cms.OutputModule("PoolOutputModule"
                               , outputCommands = cms.untracked.vstring("drop *"
                                                                     #   , "keep *_simMuonRPCDigis__*"
                                                                        , "keep *_g4SimHits_MuonCSCHits_*"
                                                                        , "keep *_g4SimHits_MuonRPCHits_*"
                                                                        , "keep *_rpcRecHits__*"
                                                                        , "keep *_simCscTriggerPrimitiveDigis_MPCSORTED_*")
                               , fileName = cms.untracked.string("RPCoutput_withRecHitLCT_GEN.root")
                               , SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)


process.e = cms.EndPath(process.out)
