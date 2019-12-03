import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import UnpackedJetCollectionVInputTag

HTXSInputTags = cms.PSet( ClassificationObj = cms.InputTag("rivetProducerHTXS","HiggsClassification") )

flashggStageOneCombinedTag = cms.EDProducer("FlashggStageOneCombinedTagProducer",
                               DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                               SystLabel=cms.string(""),
                               MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                               VBFMVAResultTag=cms.InputTag('flashggVBFMVA'),
                               GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ),
                               GenJetTag = cms.InputTag("slimmedGenJets"),
                               HTXSTags     = HTXSInputTags,
                               inputTagJets = UnpackedJetCollectionVInputTag,
                               rawDiphoBounds = cms.vdouble(),
                               rawDijetBounds = cms.vdouble()
                               )
