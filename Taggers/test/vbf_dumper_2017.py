#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
from flashgg.Systematics.SystematicDumperDefaultVariables import jetStudyVariables,minimalVariablesStage1

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
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0_20"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0_25"):
        process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0_26"):
        process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
elif os.environ["CMSSW_VERSION"].count("CMSSW_9_4"):
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v12'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.maxEvents   = cms.untracked.PSet( input  = cms.untracked.int32( 100 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)

process.load("flashgg.Taggers.flashggTagSequence_cfi")
process.load("flashgg.Taggers.flashggTagTester_cfi")

mva_wp = {
    "none"  : [
        [],
        [],
        []
    ],
    "tight" : [
        [0.69, -0.35, -0.26, -0.21],
        [0.86, -0.1 , -0.05, -0.01],
        [0.95,  0.28,  0.31,  0.28]
    ],
    "medium": [
        [0.18, -0.55, -0.42, -0.36],
        [0.61, -0.35, -0.23, -0.17],
        [0.87,  0.03,  0.13,  0.12]
    ],
    "loose" :[
        [-0.97, -0.68, -0.53, -0.47],
        [-0.89, -0.52, -0.38, -0.3 ],
        [-0.56, -0.17, -0.04, -0.01],
    ],
    "forward_tight" : [
        [-1, -0.35, -0.26, -0.21],
        [-1, -0.1 , -0.05, -0.01],
        [-1,  0.28,  0.31,  0.28]
    ],
    "forward_medium": [
        [-1, -0.55, -0.42, -0.36],
        [-1, -0.35, -0.23, -0.17],
        [-1,  0.03,  0.13,  0.12]
    ],
    "forward_loose" :[
        [-1, -0.68, -0.53, -0.47],
        [-1, -0.52, -0.38, -0.3 ],
        [-1, -0.17, -0.04, -0.01],
    ]
}


# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
# Register forwardJetRMSCut to be used from customize
customize.options.register('forwardJetRMSCut',
                           3.00,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.float,
                           'forwardJetRMSCut')

customize.options.register('pujidWP',
                           'tight', 
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.string,
                           'pujidWP')

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
customize.parse()
print "customize.processId:",customize.processId

modifyTagSequenceForSystematics(process,jetSystematicsInputTags,2)

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
    customizeSystematicsForData(process)
else:
    print "Background MC, so store mgg and central only"
    variablesToUse = minimalNonSignalVariables

    if customize.dumpJetSysTrees:
        print "Running jet systematics and putting them in ntuples because doJetSystTrees is set"
        for direction in ["Up","Down"]:
            jetsystlabels.append("JEC%s01sigma" % direction)
            jetsystlabels.append("JER%s01sigma" % direction)
            jetsystlabels.append("PUJIDShift%s01sigma" % direction)
        systlabels += jetsystlabels
        for direction in ["Up","Down"]:
            phosystlabels += ["MvaShift%s01sigma" % direction,
                           "SigmaEOverEShift%s01sigma" % direction
                           ]
        systlabels += phosystlabels
    else:
        print "Background MC, so store mgg and central only"
        customizeSystematicsForBackground(process)

#        KEY: TTreeDY_13TeV_VBFDiJet;68DY_13TeV_VBFDiJet
#        KEY: TTreeDY_13TeV_VBFDiJet_JECUp01sigma;68DY_13TeV_VBFDiJet_JECUp01sigma
#        KEY: TTreeDY_13TeV_VBFDiJet_JERUp01sigma;68DY_13TeV_VBFDiJet_JERUp01sigma
#        KEY: TTreeDY_13TeV_VBFDiJet_MvaShiftUp01sigma;67DY_13TeV_VBFDiJet_MvaShiftUp01sigma
#        KEY: TTreeDY_13TeV_VBFDiJet_PUJIDShiftUp01sigma;68DY_13TeV_VBFDiJet_PUJIDShiftUp01sigma
#        KEY: TTreeDY_13TeV_VBFDiJet_SigmaEOverEShiftUp01sigma;68DY_13TeV_VBFDiJet_SigmaEOverEShiftUp01sigma

      
print "--- Turning on central value for UnmatchedPUweight---"
for i in range(len(jetSystematicsInputTags)):
    prodname = 'flashggJetSystematics%i'%i
    vpset = getattr(process,prodname).SystMethods
    for pset in vpset:
        syst = pset.Label.value()
        if syst.count("UnmatchedPUWeight"):
            pset.ApplyCentralValue = False # default to false
            pset.Debug = False

#process.flashggVBFTag.GetQCDWeights = True
#process.flashggVBFTag.GetJetVetoWeights = True
        
print "--- Systematics  with independent collections ---"
print systlabels
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"

# ========================================================================
# Dumper section
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180605_202403/0000/myMicroAODOutputFile_24.root"
"root://eoscms.cern.ch//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_1/3_1_1/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIIFall17-3_1_1-3_1_1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180723_162219/0000/myMicroAODOutputFile_2.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180706_100515/0000/myMicroAODOutputFile_607.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DoubleEG/RunIIFall17-3_1_0-3_1_0-v0-Run2017B-31Mar2018-v1/180606_155530/0000/myMicroAODOutputFile_39.root"
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
process.vbfTagDumper.src = "flashggSystTagMerger"

# Use JetID
process.flashggVBFMVA.UseJetID      = cms.bool(True) #fixme
process.flashggVBFMVA.JetIDLevel    = cms.string("Tight2017") #cms.string("Loose")
process.flashggVBFMVA.DrJetPhoton   = cms.double(0.4) # this is the right number

process.flashggVBFTag.Boundaries             = cms.vdouble(-2.0,2.0)
process.flashggVBFTag.SetArbitraryNonGoldMC  = cms.bool(False)
process.flashggVBFTag.DropNonGoldData        = cms.bool(False)
process.flashggVBFTag.RequireVBFPreselection = cms.bool(False)


process.flashggVBFMVA.rmsforwardCut = cms.double(customize.forwardJetRMSCut)
process.flashggVBFMVA.pujidWpPtBin1 = cms.vdouble(mva_wp[customize.pujidWP][0])
process.flashggVBFMVA.pujidWpPtBin2 = cms.vdouble(mva_wp[customize.pujidWP][1])
process.flashggVBFMVA.pujidWpPtBin3 = cms.vdouble(mva_wp[customize.pujidWP][2])

print '------------------------------------------------------------'
print ' PUJID Working point    ::' , customize.pujidWP
print '------------------------------------------------------------'
print ' PUJID wp cuts          ::' , mva_wp[customize.pujidWP]
print '------------------------------------------------------------'
print ' formward RMS Cut value ::' , customize.forwardJetRMSCut
print '------------------------------------------------------------'
print ' running on Zee         ::' , customize.runOnZee
print '------------------------------------------------------------'

# use custum TMVA weights
# process.flashggVBFMVA.vbfMVAweightfile = cms.FileInPath("flashgg/Taggers/data/TMVAClassification_dijetMVA_76x_17_02_15_BDTG.weights.xml")
#process.flashggVBFMVA.vbfMVAweightfile = cms.FileInPath("flashgg/Taggers/data/sklearn_training_moriond17_v7.xml")
#process.flashggVBFMVA.MVAMethod        = cms.string("BDTG")
#process.flashggDiPhotonMVA.diphotonMVAweightfile = cms.FileInPath("flashgg/Taggers/data/Flashgg_DiPhoton_Moriond17.weights.xml")
# QCD Recovery
# process.flashggVBFMVA.merge3rdJet   = cms.untracked.bool(False)
# process.flashggVBFMVA.thirdJetDRCut = cms.untracked.double(1.5)

# run on Drell-Yan
if customize.runOnZee:
    process.flashggPreselectedDiPhotons.variables =  cms.vstring('pfPhoIso03',
                                                                 'trkSumPtHollowConeDR03',
                                                                 'full5x5_sigmaIetaIeta',
                                                                 'full5x5_r9',
                                                                 '1-passElectronVeto')

#for syst in jetsystlabels:
#    getattr(process,"flashggPreselectedDiPhotons%s"%syst).variables = process.flashggPreselectedDiPhotons.variables
# get the variable list
import flashgg.Taggers.VBFTagVariables as var
new_variables = [
    "n_jet_30             := VBFMVA.n_rec_jets",
    "dijet_jet1_RMS       := leading_rms",
    "dijet_jet2_RMS       := subLeading_rms",
    "dijet_jet1_QGL       := leading_QGL",
    "dijet_jet2_QGL       := subLeading_QGL",
    "dijet_jet1_pujid_mva := leading_pujidMVA()",
    "dijet_jet2_pujid_mva := subleading_pujidMVA()",
    "dipho_pt             := diPhoton.pt",
    "dijet_pt             := VBFMVA.dijet_pt"#,
#    "n_gen_jet_30         := nGenJet30()"
]

matching_photon = [
    "dijet_jet1_match := leadingJet_match",
    "dijet_jet2_match := subLeadingJet_match",
    "prompt_pho_1 := diPhoton.leadingPhoton.genMatchType()",
    "prompt_pho_2 := diPhoton.subLeadingPhoton.genMatchType()"
]

cloneTagSequenceForEachSystematic(process,
                                  systlabels=systlabels,
                                  phosystlabels=phosystlabels,
                                  jetsystlabels=jetsystlabels,
                                  jetSystematicsInputTags=jetSystematicsInputTags,
                                  ZPlusJetMode=2)

all_variables = jetStudyVariables
all_variables.append("prefireProbability := weight(\"prefireProbability\")")

if customize.processId != "Data":
    all_variables += minimalVariablesStage1
    all_variables += [
                     #STSX 1.1: gen-level quantities to define STXS 1.1 bin
                    "gen_pTH        := tagTruth().HTXSpTH()",
                    "n_gen_jets     := tagTruth().HTXSnjets()",
                    "gen_dijet_Mjj  := VBFMVA().gen_dijet_Mjj",
                    "gen_ptHjj      := VBFMVA().gen_ptHjj",
                    "gen_njets_vbfmva := VBFMVA().gen_njets_vbfmva",
                    "gen_dijet_LeadJPt := VBFMVA().gen_dijet_LeadJPt",
                    "gen_dijet_SubJPt := VBFMVA().gen_dijet_SubJPt",
                    "gen_dijet_abs_dEta := VBFMVA().gen_dijet_abs_dEta",
                    "gen_dijet_centrality_gg := VBFMVA().gen_dijet_centrality_gg",
                    "gen_dijet_dphi_trunc := VBFMVA().gen_dijet_dphi_trunc",
                    "gen_dijet_dphi := VBFMVA().gen_dijet_dphi",
                    "gen_dijet_minDRJetPho := VBFMVA().gen_dijet_minDRJetPho",
                    "gen_leadPho_PToM := VBFMVA().gen_leadPho_PToM",
                    "gen_sublPho_PToM := VBFMVA().gen_sublPho_PToM"
    ] 
else:
    all_variables += minimalNonSignalVariables

cats = []

if customize.dumpJetSysTrees and customize.processId != "Data" :
    for syst in (jetsystlabels+phosystlabels):
        systcutstring = "hasSyst(\"%s\") "%syst
        cats += [
            ("VBFDiJet_%s"%syst,"%s"%systcutstring,0)]#,
cats += [
    ("VBFDiJet","1",0)
]

cfgTools.addCategories(process.vbfTagDumper,
                       cats,
                       variables  = all_variables,
                       histograms = []
)
process.vbfTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
if customize.runOnZee:
    if customize.processId == "Data":
        #process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v**") )
        process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_WPTight_Gsf_v*") )
else:
    #process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"))
    process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"))

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

# Split out prompt-fake or fake-fake
process.finalFilter = cms.Sequence()
if (customize.processId.count("qcd") or customize.processId.count("gjet")) and customize.processId.count("fake"):
    process.load("flashgg/Systematics/PromptFakeFilter_cfi")
    process.finalFilter += process.PromptFakeFilter
    if (customize.processId.count("promptfake")):
        process.PromptFakeFilter.doPromptFake = cms.bool(True)
        process.PromptFakeFilter.doFakeFake =cms.bool(False)
    elif (customize.processId.count("fakefake")):
        process.PromptFakeFilter.doPromptFake =cms.bool(False)
        process.PromptFakeFilter.doFakeFake =cms.bool(True)
    else:
        raise Exception,"Mis-configuration of python for prompt-fake filter"

process.p = cms.Path(process.dataRequirements
                     * process.genFilter
                     * process.flashggUpdatedIdMVADiPhotons
                     * process.flashggDiPhotonSystematics
                     * process.flashggMetSystematics
                     * process.flashggMuonSystematics
                     * process.flashggElectronSystematics
                     * (process.flashggUnpackedJets
                        * process.ak4PFCHSL1FastL2L3CorrectorChain
                        * process.jetSystematicsSequence)
                     * (process.flashggTagSequence
                        + process.systematicsTagSequences)
                     * process.flashggSystTagMerger
                     * process.finalFilter
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
customize.setDefault("maxEvents"  , 100   )
customize.setDefault("targetLumi" ,1.00e+3)

# call the customization
customize(process)
