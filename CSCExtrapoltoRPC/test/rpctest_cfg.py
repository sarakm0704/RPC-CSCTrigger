import FWCore.ParameterSet.Config as cms

process = cms.Process("rpcNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"
#process.load("Geometry.CSCGeometry.cscGeometry_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/CSCExtrapoltoRPC/test/rpcoutput-testone.root'
#        'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/data/loggs/RPCoutput_withRecHitCSC_101.root',
    )
)

process.rpcntupler = cms.EDAnalyzer("CSCExtrapoltoRPC",
#   rpcRecHits = cms.InputTag('rpcRecHits'),
#   simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
   simMuonRPCDigis = cms.InputTag('rpcRecHits'),
   simCSCTriggerpreDigis = cms.InputTag('simCscTriggerPrimitiveDigis') 
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('rpcNtuple1.root')
 )

process.p = cms.Path(process.rpcntupler)
