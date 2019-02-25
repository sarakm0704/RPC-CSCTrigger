import FWCore.ParameterSet.Config as cms

process = cms.Process("RCNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/j/jichoi/public/RPC+CSCTrigger/skimmed_output/RPCoutput_withRecHitCSC.root'
    )
)

process.rcntupler = cms.EDAnalyzer("HitAnalyer",
   cscSegments = cms.InputTag('cscSegments'),
   tracks = cms.InputTag('tracks'),
   rpcRecHits = cms.InputTag('rpcRecHits'),
   simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
   simCSCTriggerpreDigis = cms.InputTag('simCscTriggerPrimitiveDigis') 
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('rpcNtuple.root')
 )

process.p = cms.Path(process.rcntupler)
