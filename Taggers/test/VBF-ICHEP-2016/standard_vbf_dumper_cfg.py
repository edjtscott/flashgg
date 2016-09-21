#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables

# ========================================================================
# SYSTEMATICS SECTION
process = cms.Process("FLASHggSyst")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
if   os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.maxEvents   = cms.untracked.PSet( input  = cms.untracked.int32( 50000 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 50000 )

from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)

# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
# Register forwardJetRMSCut to be used from customize

customize.options.register('forwardJetRMSCut',
                           0.03,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.float,
                           'forwardJetRMSCut')

customize.options.register('runOnZee',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'runOnZee')

customize.options.register('dumpJetSysTrees',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'dumpJetSysTrees')

customize.options.register('QCDParam',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'QCDParam')

customize.options.register('Debug',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'Debug')

customize.options.register('DoPtReweighting',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'DoPtReweighting')

customize.parse()
print "customize.processId:",customize.processId

process.load("flashgg.Taggers.flashggTagSequence_cfi")
process.load("flashgg.Taggers.flashggTagTester_cfi")

print "-----------------------------------------"
from datetime import datetime
myrandseed = int(datetime.now().strftime('%H%M%S'))
print " ramdom seed   == ", myrandseed
print "-----------------------------------------"

if customize.QCDParam :
    process.load("flashgg.Taggers.flashggPromptFakeTagSequence_cfi")
    process.load("flashgg.MicroAOD.flashggParameterisedFakePhotons_cfi")
    process.load("flashgg.MicroAOD.flashggParameterisedPromptFakeDiPhotons_cfi")
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

    # include line below to omit pt reweighting. Only needs to be done when making templates
    if customize.DoPtReweighting==False:
	process.flashggParameterisedFakePhotons.doPtReweighting = cms.untracked.bool(False)

    if customize.Debug:
	process.flashggParameterisedFakePhotons.debug = cms.untracked.bool(True)
	process.flashggParameterisedDiPhotonMVA.debug = cms.untracked.bool(True)
    

    process.flashggDiPhotonSystematics.src =  cms.InputTag("flashggParameterisedPromptFakeDiPhotons")
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
    
    

# Keep an extra category as 'would go elsewhere instead', ignore preselection

process.flashggVBFTag.Boundaries             = cms.vdouble(-999.,0.62, 0.94)
process.flashggVBFTag.SetArbitraryNonGoldMC  = cms.bool(False)
process.flashggVBFTag.DropNonGoldData        = cms.bool(False)
process.flashggVBFTag.RequireVBFPreselection = cms.bool(False)

process.flashggVBFTagMerger = cms.EDProducer("VBFTagMerger",src=cms.VInputTag("flashggVBFTag"))
modifyTagSequenceForSystematics(process,jetSystematicsInputTags,customize.dumpJetSysTrees)

systlabels    = [""]
phosystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels  = []


useEGMTools(process)

#== Only run systematics for signal events
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

# load the correctors
process.load("JetMETCorrections.Configuration.JetCorrectors_cff")

if customize.processId == "Data":
    print "Data, so turn of all shifts and systematics, with some exceptions"
    variablesToUse = minimalNonSignalVariables

    #process.flashggElectronSystematics.src = cms.InputTag("flashggElectrons")
    process.flashggElectronSystematics.src = cms.InputTag("flashggSelectedElectrons")
    customizeSystematicsForData(process)
else:
    print "Background MC, so store mgg and central only"
    variablesToUse = minimalNonSignalVariables

    if customize.dumpJetSysTrees:
        print "Running jet systematics and putting them in ntuples because doJetSystTrees is set"
        for direction in ["Up","Down"]:
            jetsystlabels.append("JEC%s01sigma" % direction)
            jetsystlabels.append("JER%s01sigma" % direction)
            jetsystlabels.append("RMSShift%s01sigma" % direction)
        systlabels += jetsystlabels
    else:
        print "Background MC, so store mgg and central only"
        customizeSystematicsForBackground(process)
        
print "--- Systematics  with independent collections ---"
print systlabels
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"

cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,jetsystlabels,jetSystematicsInputTags,customize.dumpJetSysTrees)

# ========================================================================
    
# Dumper section
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
                                #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_0_0-25ns/2_0_0/VBFHToGG_M-120_13TeV_powheg_pythia8/RunIISpring16DR80X-2_0_0-25ns-2_0_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/160524_093617/0000/myMicroAODOutputFile_1.root"
                                 #"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_0_0-25ns/2_0_0/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISpring16DR80X-2_0_0-25ns-2_0_0-v0-RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/160524_084452/0000/myMicroAODOutputFile_110.root"
                                 #"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160127_023814/0000/myMicroAODOutputFile_1.root"
                                 "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_145437/0000/myMicroAODOutputFile_23.root"
                             )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"))

import flashgg.Taggers.dumperConfigTools as cfgTools
from   flashgg.Taggers.tagsDumpers_cfi   import createTagDumper

process.vbfTagDumper = createTagDumper("VBFTag")
process.vbfTagDumper.dumpTrees     = True
process.vbfTagDumper.dumpHistos    = True
process.vbfTagDumper.dumpWorkspace = False

# Use JetID
process.flashggVBFMVA.UseJetID      = cms.bool(True)
process.flashggVBFMVA.JetIDLevel    = cms.string("Tight")
#process.flashggVBFMVA.JetIDLevel    = cms.string("Loose")

process.flashggVBFTag.Boundaries  = cms.vdouble(-2)
process.systematicsTagSequences   = cms.Sequence()

process.flashggVBFMVA.rmsforwardCut = cms.double(customize.forwardJetRMSCut)
print '------------------------------------------------------------'
print ' formward RMS Cut value ::' , customize.forwardJetRMSCut
print '------------------------------------------------------------'
print ' running on Zee         ::' , customize.runOnZee
print '------------------------------------------------------------'

# run on Drell-Yan 
if customize.runOnZee:
    process.flashggPreselectedDiPhotons.variables =  cms.vstring('pfPhoIso03', 
                                                                 'trkSumPtHollowConeDR03', 
                                                                 'full5x5_sigmaIetaIeta', 
                                                                 'full5x5_r9', 
                                                                 '1-passElectronVeto')

for syst in jetsystlabels:
    getattr(process,"flashggPreselectedDiPhotons%s"%syst).variables = process.flashggPreselectedDiPhotons.variables
    
# get the variable list
import flashgg.Taggers.VBFTagVariables as var
new_variables = [
    "n_jets           := VBFMVA.n_rec_jets",
    "dijet_jet1_RMS   := leading_rms",
    "dijet_jet2_RMS   := subLeading_rms",
    "dijet_jet1_QGL   := leading_QGL",
    "dijet_jet2_QGL   := subLeading_QGL",
    "dipho_fakeWeight := fakeWeight",
    "dipho_leadPt     := diPhoton.leadingPhoton.pt()",
    "dipho_subPt      := diPhoton.subLeadingPhoton.pt()",
    "dipho_leadEnergy := diPhoton.leadingPhoton.energy()",
    "dipho_subEnergy  := diPhoton.subLeadingPhoton.energy()",
    "dipho_leadIDMVA  := diPhoton.leadPhotonId()",
    "dipho_subIDMVA   := diPhoton.subLeadPhotonId()",
    "dipho_vtxprob    := diPhotonMVA.vtxprob",
    "dipho_CosPhi     := diPhotonMVA.CosPhi",
    "dipho_sigmarv    := diPhotonMVA.sigmarv",
    "dipho_sigmawv    := diPhotonMVA.sigmawv",
    "dipho_leadEta    := diPhotonMVA.leadeta",
    "dipho_subEta     := diPhotonMVA.subleadeta"
]

matching_photon = [
    #"dijet_jet1_match := leadingJet_match",
    #"dijet_jet2_match := subLeadingJet_match",
    #"qcd_pt1_weight   := diPhoton().leadingPhoton().weight(\"fakePtReweight\")",
    #"qcd_pt2_weight   := diPhoton().subLeadingPhoton().weight(\"fakePtReweight\")"#,
    #"qcd_pt2_weight   := diPhoton().subLeadingPhoton().centralWeight()"#,
    "qcd_pt1_weight  := qcd_leadPho_pt_weight()",
    "qcd_pt2_weight  := qcd_sublPho_pt_weight()",
    "prompt_pho_1    := diPhoton.leadingPhoton.genMatchType()",
    "prompt_pho_2    := diPhoton.subLeadingPhoton.genMatchType()",
    "fakeGenJetEnergy := genJetNearestFake_ptr().energy()",
    "fakeGenJetEta    := genJetNearestFake_ptr().eta()",
    "fakeGenJetPhi    := genJetNearestFake_ptr().phi()"
]

all_variables = var.dipho_variables + var.dijet_variables +  new_variables

if customize.processId != "Data":
    all_variables += var.truth_variables + matching_photon
    
cats = [
    ("VBFDiJet","1",0)#,
]
if customize.dumpJetSysTrees :
    for syst in jetsystlabels:
        systcutstring = "hasSyst(\"%s\") "%syst
        cats += [
            ("VBFDiJet_%s"%syst,"%s"%systcutstring,0)]#,
        

cfgTools.addCategories(process.vbfTagDumper,
                       cats,
                       variables  = all_variables,
                       histograms = []
)
process.vbfTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
if customize.runOnZee:
    if customize.processId == "Data":
        process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*") )
        #else:
        #process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WP75_Gsf_v*") )
else:
    process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"))
    
process.options      = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits")

process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter

# Split WH and ZH
process.genFilter = cms.Sequence()
if (customize.processId.count("wh") or customize.processId.count("zh")) and not customize.processId.count("wzh"):
    process.load("flashgg/Systematics/VHFilter_cfi")
    process.genFilter += process.VHFilter
    process.VHFilter.chooseW = bool(customize.processId.count("wh"))
    process.VHFilter.chooseZ = bool(customize.processId.count("zh"))


print '----------- check source ------------------'
for p in dir(process):
    if 'flashgg' not in p : continue
    dd = getattr(process, '%s' % p)
    if 'src' in dir(dd):
        print '--- process :', p , ' :src: ', dd.src
print '------------------------------------------'

if customize.QCDParam :
    print '----------- MVAresultTag ------------------'
    for p in dir(process):
        if 'flashgg' not in p : continue
        dd = getattr(process, '%s' % p)
        if 'MVAResultTag' in dir(dd):
            dd.MVAResultTag = cms.InputTag('flashggParameterisedDiPhotonMVA')
            
    print '------------------------------------------'
    process.NewTagSeg  =  cms.Sequence(
        process.flashggPreselectedDiPhotons
        * process.flashggParameterisedDiPhotonMVA
        * process.flashggUnpackedJets
        * process.flashggVBFMVA
        * process.flashggVBFDiPhoDiJetMVA
        * (process.flashggUntagged
           + process.flashggVBFTag
           #+ process.flashggTTHLeptonicTag
           #+ process.flashggTTHHadronicTag
           #+ process.flashggVHEtTag
           #+ process.flashggVHLooseTag
           #+ process.flashggVHTightTag
           #+ process.flashggVHHadronicTag
        )
        * process.flashggTagSorter
    )

    process.p = cms.Path( process.dataRequirements
                          * process.genFilter
                          * process.flashggParameterisedFakePhotons
                          * process.flashggParameterisedPromptFakeDiPhotons
                          * process.flashggDiPhotonSystematics
                          * process.flashggMuonSystematics
                          * process.flashggElectronSystematics
                          * (  process.flashggUnpackedJets
                               * process.ak4PFCHSL1FastL2L3CorrectorChain
                               * process.jetSystematicsSequence)
                          * (
                              process.NewTagSeg
                              + process.systematicsTagSequences
                          ) 
                          #* process.flashggVBFTagMerger
                          * process.vbfTagDumper
    )

else:
    process.p = cms.Path(process.dataRequirements
                         * process.genFilter
                         * process.flashggUpdatedIdMVADiPhotons
                         * process.flashggDiPhotonSystematics
                         * process.flashggMuonSystematics
                         * process.flashggElectronSystematics
                         * (process.flashggUnpackedJets
                            * process.ak4PFCHSL1FastL2L3CorrectorChain
                            * process.jetSystematicsSequence)
                         * (process.flashggTagSequence
                            + process.systematicsTagSequences)
                         * process.flashggVBFTagMerger
                         * process.vbfTagDumper
    )
    
print "--- Dumping modules that take diphotons as input: ---"
mns = process.p.moduleNames()
for mn in mns:
    module = getattr(process,mn)
    if hasattr(module,"src") and type(module.src) == type(cms.InputTag("")) and module.src.value().count("DiPhoton"):
        print str(module),module.src
    elif hasattr(module,"DiPhotonTag"):
        print str(module),module.DiPhotonTag
print
printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents"  ,50000   )
customize.setDefault("targetLumi" ,1.00e+3)

# call the customization
customize(process)
