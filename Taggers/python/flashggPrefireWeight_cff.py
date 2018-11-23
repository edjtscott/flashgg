import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

# prefire weight calculator
flashggPrefireWeight = cms.EDProducer('FlashggPrefireWeightProducer',
                               DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                               inputTagJets= UnpackedJetCollectionVInputTag,
                               photonFilename = cms.FileInPath("flashgg/Taggers/data/L1prefiring_photon_2017BtoF.root"),
                               jetFilename = cms.FileInPath("flashgg/Taggers/data/L1prefiring_jet_2017BtoF.root")
)
