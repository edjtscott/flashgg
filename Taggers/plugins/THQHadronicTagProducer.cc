#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/THQHadronicTag.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

#include "flashgg/Taggers/interface/LeptonSelection2018.h"
#include "flashgg/DataFormats/interface/Met.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/RefToPtr.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "flashgg/Taggers/interface/TTH_DNN_Helper.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>

#include "TMath.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace edm;

namespace flashgg {

    class THQHadronicTagProducer : public EDProducer
    {

    public:
        THQHadronicTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int  chooseCategory( float );
        int  computeStage1Kinematics( const THQHadronicTag );
        int  chooseCategory_pt( float, float );

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<double> rhoTag_;
        EDGetTokenT<edm::TriggerResults> triggerRECO_;
        string systLabel_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
        unsigned int inputJetsCollSize_;
        bool applyMETfilters_;

        //---thresholds---
        //---photons
        double MVAThreshold_;
        double MVATTHHMVAThreshold_;
        double PhoMVAThreshold_;
        double leadPhoPtThreshold_;
        bool   leadPhoUseVariableTh_;
        double leadPhoOverMassThreshold_;
        double leadPhoOverMassTTHHMVAThreshold_;
        double subleadPhoPtThreshold_;
        bool   subleadPhoUseVariableTh_;
        double subleadPhoOverMassThreshold_;

        //---jets
        double jetPtThreshold_;
        double jetEtaThreshold_;
        double dRJetPhoLeadCut_;
        double dRJetPhoSubleadCut_;
        vector<double> bDiscriminator_;
        double jetsNumberThreshold_;
        double bjetsNumberThreshold_;
        double bjetsLooseNumberThreshold_;
        double jetsNumberTTHHMVAThreshold_;
        double bjetsNumberTTHHMVAThreshold_;
        double bjetsLooseNumberTTHHMVAThreshold_;
        double secondMaxBTagTTHHMVAThreshold_;
        string bTag_;

        //leptons
        double MuonEtaCut_;
        double MuonPtCut_;
        double MuonIsoCut_;
        double MuonPhotonDrCut_;

        double ElePtCut_;
        std::vector<double> EleEtaCuts_;
        double ElePhotonDrCut_;
        double ElePhotonZMassCut_;
        double DeltaRTrkEle_;
        bool debug_;

        unique_ptr<TMVA::Reader> THQvsBkgMva_;
        FileInPath thqVsBkgMVAweightfile_;
        string _MVAMethod;
        unique_ptr<TMVA::Reader> THQvsTTHMva_;
        FileInPath thqVsTTHMVAweightfile_;

        int jetcount_;
        float nJets_;
        int njets_btagloose_;
        int njets_btagmedium_;
        int njets_btagtight_;
        double idmva1_;
        double idmva2_;
        float leadJetPt_;
        float leadJetBTag_;
        float subLeadJetPt_;
        float sumJetPt_;
        float maxBTagVal_;
        float secondMaxBTagVal_;
        float thirdMaxBTagVal_;
        float fourthMaxBTagVal_;

        float mindRPhoLeadJet_;
        float maxdRPhoLeadJet_;

        float minPhoID_;
        float maxPhoID_;
        float pho1_ptoM_;
        float pho2_ptoM_;
        float pho1_sceta_;
        float pho2_sceta_;
        float pho1_sigmaEOverE_;
        float pho2_sigmaEOverE_;
        float pho1_scphi_;
        float pho2_scphi_;
        float pho1_hasPixelSeed_;
        float pho2_hasPixelSeed_;
        float diPhoY_;
        float diPhoPtoM_;
        float diPhoCosPhi_;
        float nbloose_;

        float btag_1_;
        float jetPt_1_;
        float jetEta_1_;
        float jetPhi_1_;
        float btag_2_;
        float jetPt_2_;
        float jetEta_2_;
        float jetPhi_2_;
        float btag_3_;
        float jetPt_3_;
        float jetEta_3_;
        float jetPhi_3_;
        float btag_4_;
        float jetPt_4_;
        float jetEta_4_;
        float jetPhi_4_;
      
        float MET_;
      
        float thqVsBkgMvaVal_;
        float thqVsTTHMvaVal_;
        
        float maxBTagVal_noBB_;
        float secondMaxBTagVal_noBB_;

        float pho1_eta_;
        float pho2_eta_;

        float diPhoDeltaR_;
        float btag_noBB_1_;
        float btag_noBB_2_;
        float btag_noBB_3_;
        float btag_noBB_4_;

        float ht_;
        float helicity_angle_;
        float top_tag_score_;

        double thqVsBkgMvaThreshold_;
        double thqVsTTHMvaThreshold_;

    };

    THQHadronicTagProducer::THQHadronicTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag>( "MVAResultTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
        rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
	    triggerRECO_( consumes<edm::TriggerResults>(iConfig.getParameter<InputTag>("RECOfilters") ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
        _MVAMethod( iConfig.getParameter<string> ( "MVAMethod" ) )
    {
        thqVsBkgMvaThreshold_ = iConfig.getParameter<double>( "thqVsBkgMvaThreshold");
        thqVsTTHMvaThreshold_ = iConfig.getParameter<double>( "thqVsTTHMvaThreshold");

        MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold");
        MVATTHHMVAThreshold_ = iConfig.getParameter<double>( "MVATTHHMVAThreshold");
        PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold");

        leadPhoPtThreshold_ = iConfig.getParameter<double>( "leadPhoPtThreshold");
        leadPhoUseVariableTh_ = iConfig.getParameter<bool>( "leadPhoUseVariableThreshold");
        leadPhoOverMassThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassThreshold");
        leadPhoOverMassTTHHMVAThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassTTHHMVAThreshold");
        subleadPhoPtThreshold_ = iConfig.getParameter<double>( "subleadPhoPtThreshold");
        subleadPhoUseVariableTh_ = iConfig.getParameter<bool>( "subleadPhoUseVariableThreshold");
        subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold");
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
        dRJetPhoLeadCut_ = iConfig.getParameter<double>( "dRJetPhoLeadCut");
        dRJetPhoSubleadCut_ = iConfig.getParameter<double>( "dRJetPhoSubleadCut");
        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");
        jetsNumberThreshold_ = iConfig.getParameter<int>( "jetsNumberThreshold");
        bjetsNumberThreshold_ = iConfig.getParameter<int>( "bjetsNumberThreshold");
        bTag_ = iConfig.getParameter<string> ( "bTag");
        MuonEtaCut_ = iConfig.getParameter<double>( "MuonEtaCut");
        MuonPtCut_ = iConfig.getParameter<double>( "MuonPtCut");
        MuonIsoCut_ = iConfig.getParameter<double>( "MuonIsoCut");
        MuonPhotonDrCut_ = iConfig.getParameter<double>( "MuonPhotonDrCut");
 
        EleEtaCuts_ = iConfig.getParameter<std::vector<double>>( "EleEtaCuts");
        ElePtCut_ = iConfig.getParameter<double>( "ElePtCut");
        ElePhotonDrCut_ = iConfig.getParameter<double>( "ElePhotonDrCut");
        ElePhotonZMassCut_ = iConfig.getParameter<double>( "ElePhotonZMassCut");
        DeltaRTrkEle_ = iConfig.getParameter<double>( "DeltaRTrkEle");

        debug_ = iConfig.getParameter<bool>( "debug" );

        applyMETfilters_   = iConfig.getParameter<bool>( "applyMETfilters");
        bjetsLooseNumberThreshold_ = iConfig.getParameter<int>( "bjetsLooseNumberThreshold");
        jetsNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "jetsNumberTTHHMVAThreshold");
        bjetsNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "bjetsNumberTTHHMVAThreshold");
        bjetsLooseNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "bjetsLooseNumberTTHHMVAThreshold");
        secondMaxBTagTTHHMVAThreshold_ = iConfig.getParameter<double>( "secondMaxBTagTTHHMVAThreshold");

        thqVsBkgMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "thqVsBkgMVAweightfile" ); 
        thqVsTTHMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "thqVsTTHMVAweightfile" );

        nJets_ = 0;
        leadJetPt_ = 0.;
        leadJetBTag_ = -1.;
        subLeadJetPt_ = 0.;
        sumJetPt_ = 0.;

        maxBTagVal_ = -3.;
        secondMaxBTagVal_ = -3.;
        thirdMaxBTagVal_ = -3.;
        fourthMaxBTagVal_ = -3.;

        mindRPhoLeadJet_ = -999;
        maxdRPhoLeadJet_= -999;

        minPhoID_= -999.;
        maxPhoID_= -999.;
        pho1_ptoM_= -999.;
        pho2_ptoM_= -999.;
        pho1_sceta_= -999.;
        pho2_sceta_= -999.;
        pho1_scphi_= -999.;
        pho2_scphi_= -999.;
        pho1_sigmaEOverE_= -999.;
        pho2_sigmaEOverE_= -999.;

        pho1_hasPixelSeed_=-999;
        pho2_hasPixelSeed_=-999;

        diPhoY_= -999.;
        diPhoPtoM_= -999.;
        diPhoCosPhi_= -999.;
        nbloose_=-999;

        btag_1_=-1;
        jetPt_1_=-1;
        jetEta_1_=-6;
        jetPhi_1_=-6;
        btag_2_=-1;
        jetPt_2_=-1;
        jetEta_2_=-6;
        jetPhi_2_=-6;
        btag_3_=-1;
        jetPt_3_=-1;
        jetEta_3_=-6;
        jetPhi_3_=-6;
        btag_4_=-1;
        jetPt_4_=-1;
        jetEta_4_=-6;
        jetPhi_4_=-6;
                
        MET_=-1;

        maxBTagVal_noBB_ = -3.;
        secondMaxBTagVal_noBB_ = -3.; 

        pho1_eta_ = -999.;
        pho2_eta_ = -999.;

        diPhoDeltaR_ = -999.;
        
        btag_noBB_1_ = -1;
        btag_noBB_2_ = -1;
        btag_noBB_3_ = -1;
        btag_noBB_4_ = -1;

        ht_ = 0.;
        helicity_angle_ = -999.;

        if (_MVAMethod != ""){
            THQvsBkgMva_.reset( new TMVA::Reader( "!Color:Silent" ) );

            THQvsBkgMva_->AddVariable( "nJets", &nJets_);
            THQvsBkgMva_->AddVariable( "sumJetPt", &sumJetPt_);
            THQvsBkgMva_->AddVariable( "maxBTagVal",&maxBTagVal_);
            THQvsBkgMva_->AddVariable( "secondMaxBTagVal", &secondMaxBTagVal_);
            THQvsBkgMva_->AddVariable( "pho1_ptoM", &pho1_ptoM_);
            THQvsBkgMva_->AddVariable( "pho2_ptoM", &pho2_ptoM_);
            THQvsBkgMva_->AddVariable( "pho1_sceta", &pho1_sceta_);
            THQvsBkgMva_->AddVariable( "pho2_sceta", &pho2_sceta_);
            THQvsBkgMva_->AddVariable( "pho1_scphi", &pho1_scphi_);
            THQvsBkgMva_->AddVariable( "pho2_scphi", &pho2_scphi_);
            THQvsBkgMva_->AddVariable( "diPhoY", &diPhoY_);

            THQvsBkgMva_->AddVariable( "minPhoID", &minPhoID_);
            THQvsBkgMva_->AddVariable( "maxPhoID", &maxPhoID_);
            THQvsBkgMva_->AddVariable( "diPhoPtoM", &diPhoPtoM_);  

            THQvsBkgMva_->AddVariable( "btag_1", &btag_1_);       
            THQvsBkgMva_->AddVariable( "btag_2", &btag_2_);     
            THQvsBkgMva_->AddVariable( "btag_3", &btag_3_);   
            THQvsBkgMva_->AddVariable( "btag_4", &btag_4_);    
            THQvsBkgMva_->AddVariable( "jetPt_1", &jetPt_1_);      
            THQvsBkgMva_->AddVariable( "jetPt_2", &jetPt_2_);     
            THQvsBkgMva_->AddVariable( "jetPt_3", &jetPt_3_);     
            THQvsBkgMva_->AddVariable( "jetPt_4", &jetPt_4_);   
            THQvsBkgMva_->AddVariable( "jetEta_1", &jetEta_1_);   
            THQvsBkgMva_->AddVariable( "jetEta_2", &jetEta_2_);   
            THQvsBkgMva_->AddVariable( "jetEta_3", &jetEta_3_);   

            THQvsBkgMva_->AddVariable( "jetEta_4", &jetEta_4_);       
            THQvsBkgMva_->AddVariable( "pho1_hasPixelSeed",&pho1_hasPixelSeed_);
            THQvsBkgMva_->AddVariable( "pho2_hasPixelSeed",&pho2_hasPixelSeed_);
            THQvsBkgMva_->AddVariable( "thirdMaxBTagVal", &thirdMaxBTagVal_);          
            THQvsBkgMva_->AddVariable( "MET",&MET_);

            THQvsBkgMva_->BookMVA( _MVAMethod.c_str() , thqVsBkgMVAweightfile_.fullPath() );

            THQvsTTHMva_.reset( new TMVA::Reader( "!Color:Silent" ) );
     
            //THQvsTTHMva_->AddVariable("maxIDMVA_", &maxPhoID_);
            //THQvsTTHMva_->AddVariable("minIDMVA_", &minPhoID_);
            //THQvsTTHMva_->AddVariable("max2_btag_", &secondMaxBTagVal_noBB_);
            //THQvsTTHMva_->AddVariable("max1_btag_", &maxBTagVal_noBB_);
            //THQvsTTHMva_->AddVariable("dipho_delta_R", &diPhoDeltaR_);
            //THQvsTTHMva_->AddVariable("njets_", &nJets_);
            //THQvsTTHMva_->AddVariable("ht_", &ht_);
            //THQvsTTHMva_->AddVariable("leadptoM_", &pho1_ptoM_);
            //THQvsTTHMva_->AddVariable("subleadptoM_", &pho2_ptoM_);
            //THQvsTTHMva_->AddVariable("lead_eta_", &pho1_eta_);
            //THQvsTTHMva_->AddVariable("sublead_eta_", &pho2_eta_);
     
            //THQvsTTHMva_->AddVariable("jet1_pt_", &jetPt_1_);
            //THQvsTTHMva_->AddVariable("jet1_eta_", &jetEta_1_);
            //THQvsTTHMva_->AddVariable("jet1_btag_", &btag_noBB_1_);
            //THQvsTTHMva_->AddVariable("jet2_pt_", &jetPt_2_);
            //THQvsTTHMva_->AddVariable("jet2_eta_", &jetEta_2_);
            //THQvsTTHMva_->AddVariable("jet2_btag_", &btag_noBB_2_);
            //THQvsTTHMva_->AddVariable("jet3_pt_", &jetPt_3_);
            //THQvsTTHMva_->AddVariable("jet3_eta_", &jetEta_3_);
            //THQvsTTHMva_->AddVariable("jet3_btag_", &btag_noBB_3_);
            //THQvsTTHMva_->AddVariable("jet4_pt_", &jetPt_4_);
            //THQvsTTHMva_->AddVariable("jet4_eta_", &jetEta_4_);
            //THQvsTTHMva_->AddVariable("jet4_btag_", &btag_noBB_4_);
     
            //THQvsTTHMva_->AddVariable("leadPSV_", &pho1_hasPixelSeed_);
            //THQvsTTHMva_->AddVariable("subleadPSV_", &pho2_hasPixelSeed_);
     
            //THQvsTTHMva_->AddVariable("dipho_cosphi_", &diPhoCosPhi_);
            //THQvsTTHMva_->AddVariable("dipho_rapidity_", &diPhoY_);
            //THQvsTTHMva_->AddVariable("met_", &MET_);
     
            //THQvsTTHMva_->AddVariable("dipho_pt_over_mass_", &diPhoPtoM_);
     
            //THQvsTTHMva_->AddVariable("helicity_angle_", &helicity_angle_);
            
            THQvsTTHMva_->AddVariable( "nJets", &nJets_);
            THQvsTTHMva_->AddVariable( "sumJetPt", &sumJetPt_);
            THQvsTTHMva_->AddVariable( "maxBTagVal",&maxBTagVal_);
            THQvsTTHMva_->AddVariable( "secondMaxBTagVal", &secondMaxBTagVal_);
            THQvsTTHMva_->AddVariable( "pho1_ptoM", &pho1_ptoM_);
            THQvsTTHMva_->AddVariable( "pho2_ptoM", &pho2_ptoM_);
            THQvsTTHMva_->AddVariable( "pho1_sceta", &pho1_sceta_);
            THQvsTTHMva_->AddVariable( "pho2_sceta", &pho2_sceta_);
            THQvsTTHMva_->AddVariable( "pho1_scphi", &pho1_scphi_);
            THQvsTTHMva_->AddVariable( "pho2_scphi", &pho2_scphi_);
            THQvsTTHMva_->AddVariable( "diPhoY", &diPhoY_);

            THQvsTTHMva_->AddVariable( "minPhoID", &minPhoID_);
            THQvsTTHMva_->AddVariable( "maxPhoID", &maxPhoID_);
            THQvsTTHMva_->AddVariable( "diPhoPtoM", &diPhoPtoM_);  

            THQvsTTHMva_->AddVariable( "btag_1", &btag_1_);       
            THQvsTTHMva_->AddVariable( "btag_2", &btag_2_);     
            THQvsTTHMva_->AddVariable( "btag_3", &btag_3_);   
            THQvsTTHMva_->AddVariable( "btag_4", &btag_4_);    
            THQvsTTHMva_->AddVariable( "jetPt_1", &jetPt_1_);      
            THQvsTTHMva_->AddVariable( "jetPt_2", &jetPt_2_);     
            THQvsTTHMva_->AddVariable( "jetPt_3", &jetPt_3_);     
            THQvsTTHMva_->AddVariable( "jetPt_4", &jetPt_4_);   
            THQvsTTHMva_->AddVariable( "jetEta_1", &jetEta_1_);   
            THQvsTTHMva_->AddVariable( "jetEta_2", &jetEta_2_);   
            THQvsTTHMva_->AddVariable( "jetEta_3", &jetEta_3_);   

            THQvsTTHMva_->AddVariable( "jetEta_4", &jetEta_4_);       
            THQvsTTHMva_->AddVariable( "pho1_hasPixelSeed",&pho1_hasPixelSeed_);
            THQvsTTHMva_->AddVariable( "pho2_hasPixelSeed",&pho2_hasPixelSeed_);
            THQvsTTHMva_->AddVariable( "thirdMaxBTagVal", &thirdMaxBTagVal_);          
            THQvsTTHMva_->AddVariable( "MET",&MET_);

            THQvsTTHMva_->BookMVA(_MVAMethod.c_str(), thqVsTTHMVAweightfile_.fullPath()); 
        
        }       

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<THQHadronicTag> >();
    }

    void THQHadronicTagProducer::produce( Event &evt, const EventSetup & )
    {

        JetCollectionVector Jets(inputJetsCollSize_);
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        edm::Handle<double>  rho;
        evt.getByToken(rhoTag_,rho);
        double rho_    = *rho;

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );

        Handle<View<flashgg::Met> > METs;
        evt.getByToken( METToken_, METs );


	    //Get trigger results relevant to MET filters
        bool passMETfilters = 1;
        edm::Handle<edm::TriggerResults> triggerBits;
        evt.getByToken( triggerRECO_, triggerBits );
        const edm::TriggerNames &triggerNames = evt.triggerNames( *triggerBits );
        std::vector<std::string> flagList = {"Flag_goodVertices", "Flag_globalTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoieIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter","Flag_eeBadScFilter", "Flag_ecalBadCalibFilter"};
        if( ! evt.isRealData() ) {
          flagList = {"Flag_goodVertices", "Flag_globalTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter"};
        }
        for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ ) {
	      if(!triggerBits->accept(i)) {
	        for(size_t j=0;j<flagList.size();j++) {
	          if(flagList[j]==triggerNames.triggerName(i)) {
	            passMETfilters=0;
	            break;
	          }
	        }
	      }
	    }

        std::unique_ptr<vector<THQHadronicTag> > tthhtags( new vector<THQHadronicTag> );

        assert( diPhotons->size() == mvaResults->size() );

        for(unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );

            if(!passMETfilters && applyMETfilters_) continue;

            std::vector<edm::Ptr<flashgg::Muon> >     Muons;
            std::vector<edm::Ptr<flashgg::Electron> > Electrons;

            if(theMuons->size()>0)
                Muons = selectMuons(theMuons->ptrs(), dipho, vertices->ptrs(), MuonPtCut_, MuonEtaCut_, MuonIsoCut_, MuonPhotonDrCut_, debug_);
            if(theElectrons->size()>0)
                Electrons = selectElectrons(theElectrons->ptrs(), dipho, ElePtCut_, EleEtaCuts_, ElePhotonDrCut_, ElePhotonZMassCut_, DeltaRTrkEle_, debug_);

            if( (Muons.size() + Electrons.size()) != 0) continue;

            jetcount_ = 0;

            nJets_ = 0;
            njets_btagloose_ = 0;
            njets_btagmedium_ = 0;
            njets_btagtight_ = 0;

            idmva1_ = -999.;
            idmva2_ = -999.;

            leadJetPt_ = 0.;
            leadJetBTag_ = -1.;
            subLeadJetPt_ = 0.;
            sumJetPt_ = 0.;

            maxBTagVal_ = -3.;
            secondMaxBTagVal_ = -3.;
            thirdMaxBTagVal_ = -3.;
            fourthMaxBTagVal_ = -3.;
        
            mindRPhoLeadJet_ = -999;
            maxdRPhoLeadJet_= -999;
            
            minPhoID_= -999.;
            maxPhoID_= -999.;
            pho1_ptoM_= -999.;
            pho2_ptoM_= -999.;
            pho1_sceta_= -999.;
            pho2_sceta_= -999.;
            pho1_scphi_= -999.;
            pho2_scphi_= -999.;
            pho1_hasPixelSeed_=-999;
            pho2_hasPixelSeed_=-999;

            pho1_sigmaEOverE_= -999.;
            pho2_sigmaEOverE_= -999.;
            diPhoY_= -999.;
            diPhoPtoM_= -999.;
            diPhoCosPhi_= -999.;
            nbloose_=-999;

            btag_1_=-999;
            jetPt_1_=-999;
            jetEta_1_=-999;
            jetPhi_1_=-999;
            btag_2_=-999;
            jetPt_2_=-999;
            jetEta_2_=-999;
            jetPhi_2_=-999;
            btag_3_=-999;
            jetPt_3_=-999;
            jetEta_3_=-999;
            jetPhi_3_=-999;
            btag_4_=-999;
            jetPt_4_=-999;
            jetEta_4_=-999;
            jetPhi_4_=-999;
                    
            MET_=-1;
        
            thqVsBkgMvaVal_ = -999.;
            thqVsTTHMvaVal_ = -999.;

            maxBTagVal_noBB_ = -3.;
            secondMaxBTagVal_noBB_ = -3.;

            diPhoDeltaR_ = -999.;

            ht_ = 0.;
            helicity_angle_ = -999.;

            unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();

            std::vector<edm::Ptr<flashgg::Jet> > JetVect;
            JetVect.clear();
            std::vector<edm::Ptr<flashgg::Jet> > BJetVect;
            BJetVect.clear();
            std::vector<edm::Ptr<flashgg::Jet> > BJetTTHHMVAVect;
            BJetTTHHMVAVect.clear();
            
            std::vector<float> JetBTagVal;
            JetBTagVal.clear();
        
            idmva1_ = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            idmva2_ = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );

            if( idmva1_ <= PhoMVAThreshold_ || idmva2_ <= PhoMVAThreshold_ ) { continue; }

            edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );
        
            double leadPhoPtCut = leadPhoPtThreshold_;
            double subleadPhoPtCut = subleadPhoPtThreshold_;
            if( leadPhoUseVariableTh_ )
            { 
                leadPhoPtCut = leadPhoOverMassTTHHMVAThreshold_ * dipho->mass();
            }
            if( subleadPhoUseVariableTh_ ) { subleadPhoPtCut = subleadPhoOverMassThreshold_ * dipho->mass(); }

            if( dipho->leadingPhoton()->pt() < leadPhoPtCut || dipho->subLeadingPhoton()->pt() < subleadPhoPtCut ) { continue; }

            for( unsigned int jetIndex = 0; jetIndex < Jets[jetCollectionIndex]->size() ; jetIndex++ ) {
                edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( jetIndex );
                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }
                if(!thejet->passesJetID  ( flashgg::Tight2017 ) ) { continue; }
                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                dipho->subLeadingPhoton()->superCluster()->phi() );

                if( dRPhoLeadJet < dRJetPhoLeadCut_ || dRPhoSubLeadJet < dRJetPhoSubleadCut_ ) { continue; }
                if( thejet->pt() < jetPtThreshold_ ) { continue; }

                jetcount_++;
                nJets_ = jetcount_;
                JetVect.push_back( thejet );

                ht_ += thejet->pt();
                
                float bDiscriminatorValue = -2.;
                if(bTag_ == "pfDeepCSV") bDiscriminatorValue = thejet->bDiscriminator("pfDeepCSVJetTags:probb")+thejet->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                else  bDiscriminatorValue = thejet->bDiscriminator( bTag_ );

                float bDiscriminatorValue_noBB = -2;
                if(bTag_ == "pfDeepCSV") bDiscriminatorValue_noBB = thejet->bDiscriminator("pfDeepCSVJetTags:probb");
                else  bDiscriminatorValue_noBB = thejet->bDiscriminator( bTag_ );
         
                float jetPt = thejet->pt();
                if(jetPt > leadJetPt_){
                    if(leadJetPt_ > subLeadJetPt_) { subLeadJetPt_ = leadJetPt_; }
                    leadJetPt_ = jetPt;
                    leadJetBTag_=  bDiscriminatorValue;
                } else if(jetPt > subLeadJetPt_){
                    subLeadJetPt_ = jetPt;
                }
                sumJetPt_ += jetPt;
                
                if(bDiscriminatorValue_noBB > maxBTagVal_noBB_){
                     if(maxBTagVal_noBB_ > secondMaxBTagVal_noBB_) { secondMaxBTagVal_noBB_ = maxBTagVal_noBB_; }
                     maxBTagVal_noBB_ = bDiscriminatorValue_noBB;
     
                } else if(bDiscriminatorValue_noBB > secondMaxBTagVal_noBB_){
                     secondMaxBTagVal_noBB_ = bDiscriminatorValue_noBB;
                }
                
                if(bDiscriminatorValue > maxBTagVal_){ 

                    BJetTTHHMVAVect.insert( BJetTTHHMVAVect.begin(), thejet );
                    if(BJetTTHHMVAVect.size() >= 3){ BJetTTHHMVAVect.pop_back(); }
                    
                    if(thirdMaxBTagVal_>fourthMaxBTagVal_) { fourthMaxBTagVal_= thirdMaxBTagVal_;}
                    if(secondMaxBTagVal_>thirdMaxBTagVal_){ thirdMaxBTagVal_= secondMaxBTagVal_; }
                    if(maxBTagVal_ > secondMaxBTagVal_) { secondMaxBTagVal_ = maxBTagVal_; }

                    maxBTagVal_ = bDiscriminatorValue;

                } else if(bDiscriminatorValue > secondMaxBTagVal_){

                    if(BJetTTHHMVAVect.size() >= 2){BJetTTHHMVAVect.pop_back();} 
                    BJetTTHHMVAVect.push_back( thejet );

                    if(thirdMaxBTagVal_>fourthMaxBTagVal_) { fourthMaxBTagVal_= thirdMaxBTagVal_;}
                    if(secondMaxBTagVal_>thirdMaxBTagVal_) { thirdMaxBTagVal_= secondMaxBTagVal_;}
                    secondMaxBTagVal_ = bDiscriminatorValue;

                }else if(bDiscriminatorValue > thirdMaxBTagVal_){

                    if(thirdMaxBTagVal_>fourthMaxBTagVal_) { fourthMaxBTagVal_= thirdMaxBTagVal_;}
                    thirdMaxBTagVal_ = bDiscriminatorValue;

                }else if(bDiscriminatorValue > fourthMaxBTagVal_){
                   fourthMaxBTagVal_ = bDiscriminatorValue;
                }

                JetBTagVal.push_back( bDiscriminatorValue );

                if( bDiscriminatorValue > bDiscriminator_[0] ) njets_btagloose_++;
                if( bDiscriminatorValue > bDiscriminator_[1] ){
                    
                    njets_btagmedium_++;
                    BJetVect.push_back( thejet );
                }
                if( bDiscriminatorValue > bDiscriminator_[2] ) njets_btagtight_++;
            }

            if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
            Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
            
            BJetVect.clear();
            BJetVect = BJetTTHHMVAVect;

            if(JetVect.size()>1) mindRPhoLeadJet_=TMath::Min(deltaR( dipho->leadingPhoton()->eta(),dipho->leadingPhoton()->phi(), JetVect[0]->eta(),JetVect[0]->phi()) ,
                                                            deltaR( dipho->leadingPhoton()->eta(),dipho->leadingPhoton()->phi(), JetVect[1]->eta(),JetVect[1]->phi()));
            if(JetVect.size()>1) maxdRPhoLeadJet_=TMath::Max(deltaR( dipho->leadingPhoton()->eta(),dipho->leadingPhoton()->phi(), JetVect[0]->eta(),JetVect[0]->phi()) ,
                                                            deltaR( dipho->leadingPhoton()->eta(),dipho->leadingPhoton()->phi(), JetVect[1]->eta(),JetVect[1]->phi()));

            minPhoID_=TMath::Min( idmva1_, idmva2_);
            maxPhoID_=TMath::Max( idmva1_, idmva2_);
            pho1_ptoM_= dipho->leadingPhoton()->pt()/dipho->mass();
            pho2_ptoM_= dipho->subLeadingPhoton()->pt()/dipho->mass();
            pho1_sceta_= dipho->leadingPhoton()->superCluster()->eta();
            pho2_sceta_= dipho->subLeadingPhoton()->superCluster()->eta();
            pho1_scphi_= dipho->leadingPhoton()->superCluster()->phi();
            pho2_scphi_= dipho->subLeadingPhoton()->superCluster()->phi();

            pho1_hasPixelSeed_= dipho->leadingPhoton()->hasPixelSeed();
            pho2_hasPixelSeed_= dipho->subLeadingPhoton()->hasPixelSeed();

            pho1_sigmaEOverE_= dipho->leadingPhoton()->sigEOverE();
            pho2_sigmaEOverE_= dipho->subLeadingPhoton()->sigEOverE();

            diPhoY_= dipho->rapidity();
            diPhoPtoM_= dipho->pt()/dipho->mass();
            diPhoCosPhi_=  abs(TMath::Cos( deltaPhi( dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->phi() ) ));
            nbloose_=float(njets_btagloose_);
            MET_ = theMET->getCorPt();

            pho1_eta_= dipho->leadingPhoton()->eta();
            pho2_eta_= dipho->subLeadingPhoton()->eta();

            diPhoDeltaR_ = deltaR( dipho->leadingPhoton()->eta(),dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->eta(),dipho->subLeadingPhoton()->phi());


            if(JetVect.size()>0){
                if(bTag_ == "pfDeepCSV") btag_1_=JetVect[0]->bDiscriminator("pfDeepCSVJetTags:probb")+JetVect[0]->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                else  btag_1_ = JetVect[0]->bDiscriminator( bTag_ );
                if(bTag_ == "pfDeepCSV") btag_noBB_1_=JetVect[0]->bDiscriminator("pfDeepCSVJetTags:probb");
                else  btag_noBB_1_ = JetVect[0]->bDiscriminator( bTag_ );
                jetPt_1_=JetVect[0]->pt();
                jetEta_1_=JetVect[0]->eta();
                jetPhi_1_=JetVect[0]->phi();
            }

            if(JetVect.size()>1){
                if(bTag_ == "pfDeepCSV") btag_2_=JetVect[1]->bDiscriminator("pfDeepCSVJetTags:probb")+JetVect[1]->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                else  btag_2_ = JetVect[1]->bDiscriminator( bTag_ );
                if(bTag_ == "pfDeepCSV") btag_noBB_2_=JetVect[1]->bDiscriminator("pfDeepCSVJetTags:probb");
                else  btag_noBB_2_ = JetVect[1]->bDiscriminator( bTag_ );
                jetPt_2_=JetVect[1]->pt();
                jetEta_2_=JetVect[1]->eta();
                jetPhi_2_=JetVect[1]->phi();
            }

            if(JetVect.size()>2){
                if(bTag_ == "pfDeepCSV") btag_3_=JetVect[2]->bDiscriminator("pfDeepCSVJetTags:probb")+JetVect[2]->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                else  btag_3_ = JetVect[2]->bDiscriminator( bTag_ );
                if(bTag_ == "pfDeepCSV") btag_noBB_3_=JetVect[2]->bDiscriminator("pfDeepCSVJetTags:probb");
                else  btag_noBB_3_ = JetVect[2]->bDiscriminator( bTag_ );
                jetPt_3_=JetVect[2]->pt();
                jetEta_3_=JetVect[2]->eta();
                jetPhi_3_=JetVect[2]->phi();
            }
            if(JetVect.size()>3){
                if(bTag_ == "pfDeepCSV") btag_4_=JetVect[3]->bDiscriminator("pfDeepCSVJetTags:probb")+JetVect[3]->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                else  btag_4_ = JetVect[3]->bDiscriminator( bTag_ );
                if(bTag_ == "pfDeepCSV") btag_noBB_4_=JetVect[3]->bDiscriminator("pfDeepCSVJetTags:probb");
                else  btag_noBB_4_ = JetVect[3]->bDiscriminator( bTag_ );	    
                jetPt_4_=JetVect[3]->pt();
                jetEta_4_=JetVect[3]->eta();
                jetPhi_4_=JetVect[3]->phi();
            }

            if(secondMaxBTagVal_ >= secondMaxBTagTTHHMVAThreshold_ && njets_btagloose_ >= bjetsLooseNumberTTHHMVAThreshold_ && njets_btagmedium_ >= bjetsNumberTTHHMVAThreshold_ && jetcount_ >= jetsNumberTTHHMVAThreshold_ && _MVAMethod != ""){
                
                if(debug_){
                    cout << "THQHadronicTag -- input MVA variables : " << endl;
                    cout << "----------------------------------------" << endl;

                    cout << "nJets_ = " << nJets_ <<" jetcount"<< jetcount_<< endl;
                    cout << "sumJetPt_ = " << sumJetPt_ << endl;                    
                    cout << "maxBTagVal_ = " << maxBTagVal_ << endl;
                    cout << "secondMaxBTagVal_ = " << secondMaxBTagVal_ << endl;
                    cout << "thirdMaxBTagVal_ = " << thirdMaxBTagVal_ << endl;

                    cout << "minPhoID_ = " << minPhoID_ << endl;
                    cout << "maxPhoID_ = " << maxPhoID_ << endl;
                    cout << "pho1_ptoM_ = " << pho1_ptoM_ << endl;
                    cout << "pho2_ptoM_ = " << pho2_ptoM_ << endl;
                    cout << "pho1_sceta_ = " << pho1_sceta_ << endl;
                    cout << "pho2_sceta_ = " << pho2_sceta_ << endl;
                    cout << "pho1_scphi_ = " << pho1_scphi_ << endl;
                    cout << "pho2_scphi_ = " << pho2_scphi_ << endl;

                    cout << "pho1_hasPixelSeed_ = " << pho1_hasPixelSeed_ << endl;
                    cout << "pho2_hasPixelSeed_ = " << pho2_hasPixelSeed_ << endl;

                    cout << "diPhoY_ = " << diPhoY_ << endl;
                    cout << "diPhoPtoM_ = " << diPhoPtoM_ << endl;

                    cout << "btag_1_ = " <<btag_1_ << endl;
                    cout << "jetPt_1_ = " <<jetPt_1_ << endl;
                    cout << "jetEta_1_ = " <<jetEta_1_ << endl;

                    cout << "btag_2_ = " <<btag_2_ << endl;
                    cout << "jetPt_2_ = " <<jetPt_2_ << endl;
                    cout << "jetEta_2_ = " <<jetEta_2_ << endl;

                    cout << "btag_3_ = " <<btag_3_ << endl;
                    cout << "jetPt_3_ = " <<jetPt_3_ << endl;
                    cout << "jetEta_3_ = " <<jetEta_3_ << endl;

                    cout << "btag_4_ = " <<btag_4_ << endl;
                    cout << "jetPt_4_ = " <<jetPt_4_ << endl;
                    cout << "jetEta_4_ = " <<jetEta_4_ << endl;

                    cout << "MET_ = " <<MET_ << endl;
                    cout << "---------------------------------------" << endl;
                }
                thqVsBkgMvaVal_ = THQvsBkgMva_->EvaluateMVA( _MVAMethod.c_str() );
                if(debug_)  cout << " THQHadronicTag -- output MVA value = " << thqVsBkgMvaVal_  << endl;
            }
            
            std::vector<double> global_features;
            global_features.resize(18);
            global_features[0] = dipho->leadingPhoton()->eta();
            global_features[1] = dipho->subLeadingPhoton()->eta();
            global_features[2] = dipho->leadingPhoton()->phi();
            global_features[3] = dipho->subLeadingPhoton()->phi();
            global_features[4] = pho1_ptoM_;
            global_features[5] = pho2_ptoM_;
            global_features[6] = maxPhoID_;
            global_features[7] = minPhoID_;
            global_features[8] = log((float)theMET->pt());
            global_features[9] = (float)theMET->phi();
            global_features[10] = pho1_hasPixelSeed_;
            global_features[11] = pho2_hasPixelSeed_;
            global_features[12] = diPhoY_;
            global_features[13] = diPhoPtoM_;
            global_features[14] = diPhoDeltaR_;
            global_features[15] = maxBTagVal_noBB_;
            global_features[16] = secondMaxBTagVal_noBB_;
            global_features[17] = nJets_;

            TLorentzVector pho1, pho2;
            pho1.SetPtEtaPhiE(dipho->leadingPhoton()->pt(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), dipho->leadingPhoton()->energy());
            pho2.SetPtEtaPhiE(dipho->subLeadingPhoton()->pt(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi(), dipho->subLeadingPhoton()->energy());
            helicity_angle_ = helicity(pho1, pho2);

            thqVsTTHMvaVal_ = convert_tmva_to_prob(THQvsTTHMva_->EvaluateMVA( _MVAMethod.c_str() ));
            if (debug_) {
                cout << "TTH Hadronic Tag -- input MVA variables for Run II MVA: " << endl;
                cout << "--------------------------------------------------------" << endl;
                cout << "maxIDMVA_: " << maxPhoID_ << endl;
                cout << "minIDMVA_: " << minPhoID_ << endl;
                cout << "max1_btag_: " << maxBTagVal_noBB_ << endl;
                cout << "max2_btag_: " << secondMaxBTagVal_noBB_ << endl;
                cout << "dipho_delta_R_: " << diPhoDeltaR_ << endl;

                cout << "njets_: " << nJets_ << endl;
                cout << "ht_: " << ht_ << endl;
                cout << "leadptoM_: " << pho1_ptoM_ << endl;
                cout << "subleadptoM_: " << pho2_ptoM_ << endl;
                cout << "lead_eta_: " << pho1_eta_ << endl;
                cout << "sublead_eta_: " << pho2_eta_ << endl;

                cout << "jet1_pt_: " << jetPt_1_ << endl;
                cout << "jet1_eta_: " << jetEta_1_ << endl;
                cout << "jet1_btag_: " << btag_noBB_1_ << endl;
                cout << "jet2_pt_: " << jetPt_2_ << endl;
                cout << "jet2_eta_: " << jetEta_2_ << endl;
                cout << "jet2_btag_: " << btag_noBB_2_ << endl;
                cout << "jet3_pt_: " << jetPt_3_ << endl;
                cout << "jet3_eta_: " << jetEta_3_ << endl;
                cout << "jet3_btag_: " << btag_noBB_3_ << endl;
                cout << "jet4_pt_: " << jetPt_4_ << endl;
                cout << "jet4_eta_: " << jetEta_4_ << endl;
                cout << "jet4_btag_: " << btag_noBB_4_ << endl;

                cout << "leadPSV_: " << pho1_hasPixelSeed_ << endl;
                cout << "subleadPSV_: " << pho2_hasPixelSeed_ << endl;

                cout << "dipho_cosphi_: " << diPhoCosPhi_ << endl;
                cout << "dipho_rapidity_: " << diPhoY_ << endl;
                cout << "met_: " << MET_ << endl;
                cout << "dipho_pt_over_mass_: " << diPhoPtoM_ << endl;
                cout << "helicity_angle_: " << helicity_angle_ << endl;
                cout << "top_tag_score_: " << top_tag_score_ << endl;

                cout << endl;
                cout << "BDT Score: " << thqVsTTHMvaVal_ << endl;
            }

            global_features.clear();

            bool isTHQHadronicTagged = false;
            if ( thqVsTTHMvaVal_ > thqVsTTHMvaThreshold_ && thqVsBkgMvaVal_ > thqVsBkgMvaThreshold_ && njets_btagloose_ >= bjetsLooseNumberTTHHMVAThreshold_ && njets_btagmedium_ >= bjetsNumberTTHHMVAThreshold_ && jetcount_ >= jetsNumberTTHHMVAThreshold_ ) {
                isTHQHadronicTagged = true;
            }
            
            if( isTHQHadronicTagged ) {

                THQHadronicTag tthhtags_obj( dipho, mvares, JetVect, BJetVect );
                tthhtags_obj.setNjet( jetcount_ );
                tthhtags_obj.setNBLoose( njets_btagloose_ );
                tthhtags_obj.setNBMedium( njets_btagmedium_ );
                tthhtags_obj.setNBTight( njets_btagtight_ );
                tthhtags_obj.setDiPhotonIndex( diphoIndex );
                tthhtags_obj.setLeadJetPt( leadJetPt_ );
                tthhtags_obj.setSubLeadJetPt( subLeadJetPt_ );
                tthhtags_obj.setSumJetPt( sumJetPt_ );
                tthhtags_obj.setMaxBTagVal( maxBTagVal_ );
                tthhtags_obj.setSecondMaxBTagVal( secondMaxBTagVal_ );
                tthhtags_obj.setThirdMaxBTagVal( thirdMaxBTagVal_ );
                tthhtags_obj.setFourthMaxBTagVal( fourthMaxBTagVal_ );
                std::string syst_label = systLabel_;
                tthhtags_obj.setSystLabel( syst_label ); 
                tthhtags_obj.setMVAres(thqVsBkgMvaVal_);
                tthhtags_obj.setMET( theMET );

                tthhtags_obj.setStage1recoTag(DiPhotonTagBase::stage1recoTag::RECO_THQ_HAD);

                float btagReshapeNorm = 1.;

                for( unsigned num = 0; num < JetVect.size(); num++ ) {
                    tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagReshapeWeight");
                    btagReshapeNorm *= JetVect[num]->weight("JetBTagReshapeWeightCentral");
                }                    

                tthhtags_obj.includeWeights( *dipho );

                tthhtags_obj.setWeight( "btagReshapeNorm_TTH_HAD", btagReshapeNorm );

                tthhtags->push_back( tthhtags_obj );
            }
        }
        evt.put( std::move( tthhtags ) );
    }

}
typedef flashgg::THQHadronicTagProducer FlashggTHQHadronicTagProducer;
DEFINE_FWK_MODULE( FlashggTHQHadronicTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
