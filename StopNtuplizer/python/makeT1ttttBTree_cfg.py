import FWCore.ParameterSet.Config as cms
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices

goodOfflinePrimaryVerticesDQM = goodOfflinePrimaryVertices.clone(
        src=cms.InputTag("offlineSlimmedPrimaryVertices"))

process = cms.Process("Stops")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/b/bravo/data/susy/T1ttttB/miniAOD-prod_PAT_1.root'
    )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("stopT1ttttBNtuple.root") )

process.demo = cms.EDAnalyzer('StopNtuplizer',
        convLabel = cms.InputTag("reducedEgamma","reducedConversions")
)


process.p = cms.Path(process.demo)
