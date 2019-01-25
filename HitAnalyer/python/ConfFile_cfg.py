import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:FEAF370D-6243-E811-9145-A0369FD0B122.root'
    )
)

process.demo = cms.EDAnalyzer('HitAnalyer'
)


process.p = cms.Path(process.demo)
