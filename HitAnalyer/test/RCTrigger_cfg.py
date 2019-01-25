import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/u/user/jichoi/WORK/RPC+CSCTrigger/CMSSW_9_3_13/src/data/RPCoutput.root'
    )
)

process.demo = cms.EDAnalyzer('HitAnalyer'
#  MuonDigiRPC = cms.InputTag('simMuonRPCDigis'),
#  MuonDigiCSC = cms.InputTag('simCscTriggerPrimitiveDigis') 
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('RCNtuple.root')
 )

process.p = cms.Path(process.demo)
