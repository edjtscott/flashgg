minimalVariables = ["CMS_hgg_mass[160,100,180]:=diPhoton().mass",
                    "dZ[40,-20.,20.]:=(tagTruth().genPV().z-diPhoton().vtx().z)", # store actual value
                                                                               #when doing systematics, variables need to have a binning
                                                                               #specified, otherwise the rooDataHist end up empty.
            								       #an assert in the code prevents you from doing this.
                    "centralObjectWeight[1,-999999.,999999.] := centralWeight"]

minimalHistograms = []

minimalNonSignalVariables = ["CMS_hgg_mass[160,100,180]:=diPhoton().mass"]#,"centralObjectWeight[1,-999999.,999999.] := centralWeight"]

minimalVariablesHTXS = minimalVariables+["stage0cat[72,9.5,81.5] := tagTruth().HTXSstage0cat"]
minimalVariablesStage1 = minimalVariables+["stage1cat[39,-8.5,30.5] := tagTruth().HTXSstage1orderedBin"]

defaultVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass", 
                                    "leadPt                   :=diPhoton().leadingPhoton.pt",
                                    "subleadPt                :=diPhoton().subLeadingPhoton.pt",
                                    "diphoMVA                 :=diPhotonMVA().result",    
                                    "maxEta                   :=max(abs(diPhoton().leadingPhoton.superCluster.eta),abs(diPhoton().leadingPhoton.superCluster.eta))",
                                    "genZ           :=tagTruth().genPV().z",
                                    "vtxZ           :=diPhoton().vtx().z",
                                    "dZ             :=(tagTruth().genPV().z-diPhoton().vtx().z)"]


defaultHistograms=["CMS_hgg_mass>>mass(160,100,180)",
                                     "subleadPt:leadPt>>ptLeadvsSub(180,20,200:180,20,200)",
                                     "diphoMVA>>diphoMVA(50,0,1)",
                                     "maxEta>>maxEta[0.,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.8,2.,2.2,2.3,2.5]"
                                     ]

systematicVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass"]#,"centralObjectWeight[1,-999999.,999999.] := centralWeight"]
systematicHistograms=["CMS_hgg_mass>>mass(160,100,180)"]

systematicVariablesHTXS = systematicVariables+["stage0cat[72,9.5,81.5] := tagTruth().HTXSstage0cat"]
systematicVariablesStage1 = systematicVariables+["stage1cat[39,-8.5,30.5] := tagTruth().HTXSstage1orderedBin"]

jetStudyVariables = [
                    #global variables
                    "n_rec_jets := VBFMVA().n_rec_jets",
                    "diphoptom  := VBFDiPhoDiJetMVA().dipho_PToM",
                    "diphomva   := diPhotonMVA().result",
                    "diphomvaxgb   := 1. / ( 1. + exp( 0.5*log( 2./(diPhotonMVA().xgbResult+1.) - 1 ) ) )",
                    "dipho_pt   := VBFDiPhoDiJetMVA().dipho_PToM*diPhoton().mass",
                    #diphoton BDT inputs
                    "leadmva     := diPhotonMVA().leadmva",
                    "subleadmva  := diPhotonMVA().subleadmva",
                    "leadptom    := diPhotonMVA().leadptom",
                    "subleadptom := diPhotonMVA().subleadptom",
                    "leadeta     := diPhotonMVA().leadeta",
                    "subleadeta  := diPhotonMVA().subleadeta",
                    "CosPhi      := diPhotonMVA().CosPhi",
                    "vtxprob     := diPhotonMVA().vtxprob",
                    "sigmarv     := diPhotonMVA().sigmarv",
                    "sigmawv     := diPhotonMVA().sigmawv",
                    #jet-related variables
                    "VBFMVAValue              := VBFMVA().VBFMVAValue()",
                    "ptHjj                    := VBFMVA().ptHjj()",
                    "dijet_Mjj                := VBFMVA().dijet_Mjj",
                    "dijet_leadEta            := VBFMVA().dijet_leadEta",
                    "dijet_subleadEta         := VBFMVA().dijet_subleadEta",
                    "dijet_subsubleadEta      := VBFMVA().dijet_subsubleadEta",
                    "dijet_LeadJPt            := VBFMVA().dijet_LeadJPt",
                    "dijet_SubJPt             := VBFMVA().dijet_SubJPt",
                    "dijet_SubsubJPt          := VBFMVA().dijet_SubsubJPt",
                    "dijet_leadPUMVA          := VBFMVA().dijet_leadPUMVA",
                    "dijet_subleadPUMVA       := VBFMVA().dijet_subleadPUMVA",
                    "dijet_subsubleadPUMVA    := VBFMVA().dijet_subsubleadPUMVA",
                    "dijet_leadDeltaPhi       := VBFMVA().dijet_leadDeltaPhi",
                    "dijet_subleadDeltaPhi    := VBFMVA().dijet_subleadDeltaPhi",
                    "dijet_subsubleadDeltaPhi := VBFMVA().dijet_subsubleadDeltaPhi",
                    "dijet_leadDeltaEta       := VBFMVA().dijet_leadDeltaEta",
                    "dijet_subleadDeltaEta    := VBFMVA().dijet_subleadDeltaEta",
                    "dijet_subsubleadDeltaEta := VBFMVA().dijet_subsubleadDeltaEta",
                    #STSX 1.1: gen-level quantities to define STXS 1.1 bin
                    "gen_pTH        := tagTruth().HTXSpTH()",
                    "n_gen_jets     := tagTruth().HTXSnjets()"#,
                    #"gen_dijet_Mjj  := VBFMVA().gen_dijet_Mjj",
                    #"gen_ptHjj      := VBFMVA().gen_ptHjj"
                    ]
