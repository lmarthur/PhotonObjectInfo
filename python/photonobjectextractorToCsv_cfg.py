import FWCore.ParameterSet.Config as cms

process = cms.Process("photonext")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
     # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                    'file:myfile.root'
                                    )
                            )

process.photonextractorToCsv = cms.EDAnalyzer('PhotonInfoExtractorToCsv',
                                                InputCollection = cms.InputTag("photons"),
                                                maxNumberPhotons = cms.untracked.int32(10)
                                                )

process.p = cms.Path(process.photonextractorToCsv)
