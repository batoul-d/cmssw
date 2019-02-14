import FWCore.ParameterSet.Config as cms


pfCandComposites = cms.EDProducer('PFCandCompositeProducer',
                                  pfCandTag    = cms.InputTag('particleFlow'),
                                  compositeTag = cms.InputTag('onia2MuMuPatGlbGlb'),
                                  jpsiTrigFilter = cms.string("hltL3fL1DoubleMuOpenL3FilteredPsi"),
                                  isHI = cms.bool(True),
                                  replaceJMM = cms.bool(True),
                                  replaceDKPi = cms.bool(False)
                                  )
