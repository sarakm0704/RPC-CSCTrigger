import FWCore.ParameterSet.Config as cms

process = cms.Process("rpcNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"
#process.load("Geometry.CSCGeometry.cscGeometry_cfi")

# import of standard configurations
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

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#       'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/CSCExtrapoltoRPC/test/rpcoutput-testone.root'
#        'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/data/loggs/RPCoutput_withRecHitCSC_101.root',
#       'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/CSCExtrapoltoRPC/test/rpcJPSIoutput.root',
       'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/RPCRecHitCscProducer/RPCoutput_withRecHitLCT_GEN.root'
    )
)

process.rpcntupler = cms.EDAnalyzer("CSCExtrapoltoRPC",
#   rpcRecHits = cms.InputTag('rpcRecHits'),
#   simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
   simHitLl = cms.untracked.InputTag('g4SimHits','MuonCSCHits'),
   simMuonRPCDigis = cms.InputTag('rpcRecHits'),
   simCSCTriggerpreDigis = cms.InputTag('simCscTriggerPrimitiveDigis') 
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('rpcNtuple_sample.root')
 )

process.p = cms.Path(process.rpcntupler)
