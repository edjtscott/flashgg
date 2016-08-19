# ================================================
#   Both jet tree and VBF tree at the same time
#   Y. Haddad 01/2015, S Zenz 04/2016
# ================================================

#from sys import argv
#print argv
#myfilenum = str(argv[2])
outdir  = "/vols/cms/es811/Templates18Aug16/"
outname = "Template_file"
myrandseed = 235

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

from flashgg.MetaData.JobConfig import customize
customize.options.register('QCDParam',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'QCDParam')
customize.options.register('FileNum',
				  '8',
				  VarParsing.VarParsing.multiplicity.singleton,
				  VarParsing.VarParsing.varType.string,
				  "File number")
customize.parse()
print "customize.processId:",customize.processId

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
                            fileNames=cms.untracked.vstring("root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_145437/0000/myMicroAODOutputFile_"+customize.FileNum+".root"))
                            #fileNames=cms.untracked.vstring("/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_145437/0000/myMicroAODOutputFile_"+myfilenum+".root"))

#process.TFileService = cms.Service("TFileService", fileName = cms.string(outdir+outname+myfilenum+".root")) 
process.TFileService = cms.Service("TFileService", fileName = cms.string(outdir+outname+customize.FileNum+".root")) 

process.load("flashgg.Taggers.flashggParamTemplateMaker_cfi")

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

if customize.QCDParam:
    process.load("flashgg.Taggers.flashggPromptFakeTagSequence_cfi")
    process.load("flashgg.MicroAOD.flashggParameterisedFakePhotons_cfi")
    process.load("flashgg.MicroAOD.flashggParameterisedPromptFakeDiPhotons_cfi")
    process.RandomNumberGeneratorService   =  cms.Service ("RandomNumberGeneratorService",
                                                           flashggParameterisedFakePhotons = cms.PSet(
                                                               initialSeed = cms.untracked.uint32( myrandseed ),
                                                               engineName  = cms.untracked.string('TRandom3')
                                                           ),
                                                           flashggParameterisedDiPhotonMVA = cms.PSet(
                                                               initialSeed = cms.untracked.uint32( 45*myrandseed + 1 ),
                                                               engineName  = cms.untracked.string('TRandom3')
                                                           ),
                                                           flashggJetValidationTreeMakerPFCHS = cms.PSet(
                                                               initialSeed = cms.untracked.uint32( 5912*myrandseed + 179 ),
                                                               engineName  = cms.untracked.string('TRandom3')
                                                           ),)
    massSearchReplaceAnyInputTag( process.flashggPromptFakeTagSequence, cms.InputTag("flashggUpdatedIdMVADiPhotons"), cms.InputTag("flashggParameterisedPromptFakeDiPhotons") )
    process.flashggParamTemplateMaker.MVAResultTag = cms.InputTag("flashggParameterisedDiPhotonMVA")
    process.flashggParamTemplateMaker.doParam = cms.untracked.bool(True)
    process.flashggParameterisedFakePhotons.doPtReweighting = cms.untracked.bool(False)
    process.p = cms.Path(
        process.flashggParameterisedFakePhotons *
        process.flashggParameterisedPromptFakeDiPhotons *
        process.flashggPromptFakeTagSequence *
        process.flashggParamTemplateMaker
        )
else:
    process.load("flashgg.Taggers.flashggTagSequence_cfi")
    process.p = cms.Path( 
        process.flashggUpdatedIdMVADiPhotons *
        process.flashggTagSequence *
        process.flashggParamTemplateMaker
        )

#from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

#process.p = cms.Path( 
#    process.flashggUpdatedIdMVADiPhotons *
#    process.flashggTagSequence *
#    process.flashggParamTemplateMaker
#    )

customize(process)
