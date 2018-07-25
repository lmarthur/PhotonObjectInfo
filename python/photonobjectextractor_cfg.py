import FWCore.ParameterSet.Config as cms

process = cms.Process("photonext")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("poolsource",
     # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                    'file:myfile.root'
                                    )
                            )

process.demo = cms.EDAnalyzer('PhotonObjectInfoExtractor',
                              InputCollection = cms.InputTag("photons")
                              )

process.p = cms.Path(process.photonextractor)
