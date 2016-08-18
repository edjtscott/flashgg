import FWCore.ParameterSet.Config as cms

flashggParamTemplateMaker = cms.EDAnalyzer('FlashggParamTemplateMaker',
					    DiPhotonTag     = cms.InputTag('flashggPreselectedDiPhotons'),
					    MVAResultTag    = cms.InputTag('flashggDiPhotonMVA'),
					    )

