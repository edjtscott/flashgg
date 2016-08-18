# ================================================
#   Both jet tree and VBF tree at the same time
#   Y. Haddad 01/2015, S Zenz 04/2016
# ================================================

from sys import argv
print argv
myfilenum = str(argv[2])
outdir  = "/vols/cms/es811/Templates18Aug16/"
outname = "Template_file"

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("ParamTemplateMaker")

# +++++ the number of processed events
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1000 ) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
process.MessageLogger.cerr.threshold = 'ERROR'

# +++++ the source file
process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring("root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_145437/0000/myMicroAODOutputFile_"+myfilenum+".root"))
                            fileNames=cms.untracked.vstring("/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_145437/0000/myMicroAODOutputFile_"+myfilenum+".root"))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outdir+outname+myfilenum+".root")) 

process.load("flashgg.Taggers.flashggTagSequence_cfi")
process.load("flashgg.Taggers.flashggParamTemplateMaker_cfi")

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.p = cms.Path( 
    process.flashggUpdatedIdMVADiPhotons *
    process.flashggTagSequence *
    process.flashggParamTemplateMaker
    )

