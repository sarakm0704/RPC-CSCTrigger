import FWCore.ParameterSet.Config as cms

process = cms.Process("rpcNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/RPC-CSCTrigger/CSCExtrapoltoRPC/test/rpcoutput-testone.root'
    )
)

process.rpcntupler = cms.EDAnalyzer("CSCExtrapoltoRPC",
   cscSegments = cms.InputTag('cscSegments'),
   tracks = cms.untracked.InputTag('tracks'),
#   rpcRecHits = cms.InputTag('rpcRecHits'),
#   simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
   simMuonRPCDigis = cms.InputTag('rpcRecHits'),
   simCSCTriggerpreDigis = cms.InputTag('simCscTriggerPrimitiveDigis') 
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('rpcNtuple.root')
 )

process.p = cms.Path(process.rpcntupler)
