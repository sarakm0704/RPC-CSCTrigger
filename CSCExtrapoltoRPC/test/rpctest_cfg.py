import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Labels', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "1: 1DIGI, 4:4DIGI")
options.parseArguments()

process = cms.Process("rpcNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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
       'file:/afs/cern.ch/user/j/jichoi/public/RPC+CSCTrigger/CMSSW_10_5_0_pre1/src/RPC-CSCTrigger/CSCExtrapoltoRPC/test/output/RPCoutput_withRecHitLCT_GEN.root'
    )
)

process.rpcntupler = cms.EDAnalyzer("CSCExtrapoltoRPC",
#   rpcRecHits = cms.InputTag('rpcRecHits'),
#   simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
   simHitLl = cms.untracked.InputTag('g4SimHits','MuonCSCHits'),
   simMuonRPCDigis = cms.InputTag('rpcRecHits'),
   simCSCTriggerpreDigis = cms.InputTag('simCscTriggerPrimitiveDigis'),
   label = cms.untracked.int32(options.Labels)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('rpcNtuple_sample.root')
 )

process.p = cms.Path(process.rpcntupler)
