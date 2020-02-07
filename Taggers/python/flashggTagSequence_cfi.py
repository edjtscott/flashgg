import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggUnpackedJets
from flashgg.Taggers.flashggDiPhotonMVA_cfi import flashggDiPhotonMVA
from flashgg.Taggers.flashggVBFMVA_cff import flashggVBFMVA,flashggVBFDiPhoDiJetMVA
from flashgg.Taggers.flashggVHhadMVA_cff import flashggVHhadMVA
from flashgg.Taggers.flashggPrefireDiPhotons_cff import flashggPrefireDiPhotons
from flashgg.Taggers.flashggTags_cff import *
from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
from flashgg.Taggers.flashggTagSorter_cfi import flashggTagSorter
from flashgg.Taggers.flashggDifferentialPhoIdInputsCorrection_cfi import flashggDifferentialPhoIdInputsCorrection, setup_flashggDifferentialPhoIdInputsCorrection

def flashggPrepareTagSequence(process, options):
    setup_flashggDifferentialPhoIdInputsCorrection(process, options)
    flashggPreselectedDiPhotons.src = "flashggPrefireDiPhotons"

    if "flashggDiPhotonMVA" in options:
        flashggDiPhotonMVA.diphotonMVAweightfile = cms.FileInPath(str(options["flashggDiPhotonMVA"]["weightFile"]))
        flashggDiPhotonMVA.Version = cms.string(str(options["flashggDiPhotonMVA"]["version"]))
    if "flashggVBFMVA" in options:
        flashggVBFMVA.vbfMVAweightfile = cms.FileInPath(str(options["flashggVBFMVA"]["weightFile"]))
    if "flashggVHhadMVA" in options:
        flashggVHhadMVA.vbfMVAweightfile = cms.FileInPath(str(options["flashggVHhadMVA"]["weightFile"]))
    flashggTHQLeptonicTag.thqCatweightfile_ForNonPeakingBkg = cms.FileInPath(str(options['THQLeptonicTag']['MVAweights_VsAllBkg']))
    flashggTHQLeptonicTag.MVAThreshold_ForNonPeakingBkg = cms.double(options['THQLeptonicTag']['MVAThreshold_VsAllBkg'])
    flashggTHQLeptonicTag.thqleptonicMVAweightfile = cms.FileInPath(str(options['THQLeptonicTag']['MVAweights_VsttH']))
    flashggTHQLeptonicTag.MVAThreshold_thq = cms.double(options['THQLeptonicTag']['MVAThreshold_VsttH'])

    flashggTagSequence = cms.Sequence(flashggDifferentialPhoIdInputsCorrection
                                      * flashggPrefireDiPhotons
                                      * flashggPreselectedDiPhotons
                                      * flashggDiPhotonMVA
                                      * flashggUnpackedJets
                                      * flashggVBFMVA
                                      * flashggVHhadMVA
                                      * flashggVBFDiPhoDiJetMVA
                                      * ( flashggUntagged
                                      #                                  *( flashggSigmaMoMpToMTag
                                          + flashggVBFTag
                                          + flashggTTHDiLeptonTag
                                          + flashggTTHLeptonicTag
                      + flashggTHQLeptonicTag
#                                     + flashggTTHHadronicTTag                                      
#                                     + flashggTTHHadronicLTag                                      
                                          + flashggTTHHadronicTag
                                      #############old VH tags##############
                                      #                  + flashggVHEtTag
                                      #                  + flashggVHLooseTag
                                      #                  + flashggVHTightTag
                                      ###########updated VH tags############
                                          + flashggVHMetTag
                                          + flashggWHLeptonicTag
                                          + flashggZHLeptonicTag
                                          + flashggVHLeptonicLooseTag
                                          + flashggVHHadronicTag
                                      )
                                      * flashggTagSorter
                                  )

    return flashggTagSequence
