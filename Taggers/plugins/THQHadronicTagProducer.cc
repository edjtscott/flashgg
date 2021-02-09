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
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "flashgg/Taggers/interface/TTH_DNN_Helper.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
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
        int  chooseCategory( float, float );

        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<double> rhoTag_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        string systLabel_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

        vector<double> bDiscriminator_;

        //Thresholds
        double leadPhoOverMassThreshold_      ; 
        double subleadPhoOverMassThreshold_   ; 
        double PhoMVAThreshold_               ; 
        double MVAThreshold_                  ; 
        double leptonPtThreshold_             ; 
        vector<double> electronEtaThresholds_ ; 
        bool useElectronMVARecipe_            ; 
        bool useElectronLooseID_              ; 
        double deltaRPhoElectronThreshold_    ; 
        double DeltaRTrkElec_                 ; 
        double deltaMassElectronZThreshold_   ; 
        double muonEtaThreshold_              ; 
        double muPFIsoSumRelThreshold_        ; 
        double deltaRMuonPhoThreshold_        ; 
        double jetsNumberThreshold_           ; 
        double jetPtThreshold_                ; 
        double jetEtaThreshold_               ; 
        double deltaRPhoLeadJet_              ; 
        double deltaRPhoSubLeadJet_           ; 
        double deltaRJetMuonThreshold_        ; 
        double METThreshold_                  ; 
        double tthMvaThreshold_               ; 
        double bkgMvaThreshold_               ; 

        //MVAs vs ttH and vs bkg
        FileInPath TTHMVAweightfile_;
        FileInPath BkgMVAweightfile_;
        unique_ptr<TMVA::Reader> TTHMva_;
        unique_ptr<TMVA::Reader> BkgMva_;

        // input vars for tHq vs ttH BDT
        float _lead_ptoM;                                 
        float _sublead_ptoM;                                 
        float _leadEta;                                 
        float _subleadEta;                                 
        float _leadPhi;                                 
        float _subleadPhi;                                 
        float _leadPixelSeed;                                 
        float _subleadPixelSeed;                                 
        float _bjet1_csv;                                 
        float _bjet2_csv;                                 
        float _n_jets;                                 
        float _forward_pt;                                 
        float _forward_eta;                                 
        float _dipho_ptoM;                                 
        float _deltaR_gg;                                 
        float _dipho_rapidity;                                 
        float _dipho_cosphi;                                 
        float _jet1_bdiscriminant;                                  
        float _jet2_bdiscriminant;                                 
        float _jet3_bdiscriminant;                                 
        float _jet1_pt;                                 
        float _jet2_pt;                                 
        float _jet3_pt;                                 
        float _jet1_eta;                                 
        float _jet2_eta;                                 
        float _jet3_eta;                                 
        float _min_id;                                 
        float _max_id;                                 
        float _nbjets;                                 
        float _ncentraljets;                                 
        float _bjet1_pt;                                 
        float _dr_leadpho_fwdjet;                                 
        float _dr_subleadpho_fwdjet;                                 
        float _dr_bjet_fwdjet;                                 

        // extra input vars for tHq vs non-H bkg BDT
        float _dipho_helicity;
        float _bjet2_pt;
        float _bjet3_pt;
        float _bjet1_eta;
        float _bjet2_eta;
        float _bjet3_eta;
        float _bjet3_csv;
    };

    THQHadronicTagProducer::THQHadronicTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag> ( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {

        //Photon
        leadPhoOverMassThreshold_    = iConfig.getParameter<double>( "leadPhoOverMassThreshold");
        subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold");
        PhoMVAThreshold_             = iConfig.getParameter<double>( "PhoMVAThreshold");
        MVAThreshold_                = iConfig.getParameter<double>( "MVAThreshold");

        //Electron
        leptonPtThreshold_           = iConfig.getParameter<double>( "leptonPtThreshold");
        electronEtaThresholds_       = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
        useElectronMVARecipe_        = iConfig.getParameter<bool>("useElectronMVARecipe");
        useElectronLooseID_          = iConfig.getParameter<bool>("useElectronLooseID");
        deltaRPhoElectronThreshold_  = iConfig.getParameter<double>( "deltaRPhoElectronThreshold");
        DeltaRTrkElec_               = iConfig.getParameter<double>( "DeltaRTrkElec");
        deltaMassElectronZThreshold_ = iConfig.getParameter<double>( "deltaMassElectronZThreshold");
        
        //Muon
        muonEtaThreshold_            = iConfig.getParameter<double>( "muonEtaThreshold");
        muPFIsoSumRelThreshold_      = iConfig.getParameter<double>( "muPFIsoSumRelThreshold");
        deltaRMuonPhoThreshold_      = iConfig.getParameter<double>( "deltaRMuonPhoThreshold");

        //Jet
        jetsNumberThreshold_         = iConfig.getParameter<double>( "jetsNumberThreshold");
        jetPtThreshold_              = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_             = iConfig.getParameter<double>( "jetEtaThreshold");
        deltaRPhoLeadJet_            = iConfig.getParameter<double>( "deltaRPhoLeadJet");
        deltaRPhoSubLeadJet_         = iConfig.getParameter<double>( "deltaRPhoSubLeadJet");
        deltaRJetMuonThreshold_      = iConfig.getParameter<double>( "deltaRJetMuonThreshold");

        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");

        //Met
        METThreshold_                = iConfig.getParameter<double>( "METThreshold");

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }
        
        // tHq vs ttH MVA
        TTHMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "TTHMVAweightfile" );
        tthMvaThreshold_  = iConfig.getParameter<double>( "tthMvaThreshold" );
        BkgMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "BkgMVAweightfile" );
        bkgMvaThreshold_  = iConfig.getParameter<double>( "bkgMvaThreshold" );

        _lead_ptoM = -999.;
        _sublead_ptoM = -999.;
        _leadEta = -999.;
        _subleadEta = -999.;
        _leadPhi = -999.;
        _subleadPhi = -999.;
        _leadPixelSeed = -999.;
        _subleadPixelSeed = -999.;
        _bjet1_csv = -999.;
        _bjet2_csv = -999.;
        _n_jets = -999.;
        _forward_pt = -999.;
        _forward_eta = -999.;
        _dipho_ptoM = -999.;
        _deltaR_gg = -999.;
        _dipho_rapidity = -999.;
        _dipho_cosphi = -999.;
        _jet1_bdiscriminant = -999.;
        _jet2_bdiscriminant = -999.;
        _jet3_bdiscriminant = -999.;
        _jet1_pt = -999.;
        _jet2_pt = -999.;
        _jet3_pt = -999.;
        _jet1_eta = -999.;
        _jet2_eta = -999.;
        _jet3_eta = -999.;
        _min_id = -999.;
        _max_id = -999.;
        _nbjets = -999.;
        _ncentraljets = -999.;
        _bjet1_pt = -999.;
        _dr_leadpho_fwdjet = -999.;
        _dr_subleadpho_fwdjet = -999.;
        _dr_bjet_fwdjet = -999.;

        _dipho_helicity = -999.;
        _bjet2_pt = -999.;
        _bjet3_pt = -999.;
        _bjet1_eta = -999.;
        _bjet2_eta = -999.;
        _bjet3_eta = -999.;
        _bjet3_csv = -999.;

        TTHMva_.reset( new TMVA::Reader( "!Color:!Silent" ) );
        TTHMva_->AddVariable( "lead_ptoM",  &_lead_ptoM);
        TTHMva_->AddVariable( "sublead_ptoM",  &_sublead_ptoM);
        TTHMva_->AddVariable( "leadEta",  &_leadEta);
        TTHMva_->AddVariable( "subleadEta",  &_subleadEta);
        TTHMva_->AddVariable( "leadPhi",  &_leadPhi);
        TTHMva_->AddVariable( "subleadPhi",  &_subleadPhi);
        TTHMva_->AddVariable( "leadPixelSeed",  &_leadPixelSeed);
        TTHMva_->AddVariable( "subleadPixelSeed",  &_subleadPixelSeed);
        TTHMva_->AddVariable( "bjet1_csv",  &_bjet1_csv);
        TTHMva_->AddVariable( "bjet2_csv",  &_bjet2_csv);
        TTHMva_->AddVariable( "n_jets",  &_n_jets);
        TTHMva_->AddVariable( "forward_pt",  &_forward_pt);
        TTHMva_->AddVariable( "forward_eta",  &_forward_eta);
        TTHMva_->AddVariable( "dipho_ptoM",  &_dipho_ptoM);
        TTHMva_->AddVariable( "deltaR_gg",  &_deltaR_gg);
        TTHMva_->AddVariable( "dipho_rapidity",  &_dipho_rapidity);
        TTHMva_->AddVariable( "dipho_cosphi",  &_dipho_cosphi);
        TTHMva_->AddVariable( "jet1_bdiscriminant ",  &_jet1_bdiscriminant );
        TTHMva_->AddVariable( "jet2_bdiscriminant",  &_jet2_bdiscriminant);
        TTHMva_->AddVariable( "jet3_bdiscriminant",  &_jet3_bdiscriminant);
        TTHMva_->AddVariable( "jet1_pt",  &_jet1_pt);
        TTHMva_->AddVariable( "jet2_pt",  &_jet2_pt);
        TTHMva_->AddVariable( "jet3_pt",  &_jet3_pt);
        TTHMva_->AddVariable( "jet1_eta",  &_jet1_eta);
        TTHMva_->AddVariable( "jet2_eta",  &_jet2_eta);
        TTHMva_->AddVariable( "jet3_eta",  &_jet3_eta);
        TTHMva_->AddVariable( "min_id",  &_min_id);
        TTHMva_->AddVariable( "max_id",  &_max_id);
        TTHMva_->AddVariable( "nbjets",  &_nbjets);
        TTHMva_->AddVariable( "ncentraljets",  &_ncentraljets);
        TTHMva_->AddVariable( "bjet1_pt",  &_bjet1_pt);
        TTHMva_->AddVariable( "dr_leadpho_fwdjet",  &_dr_leadpho_fwdjet);
        TTHMva_->AddVariable( "dr_subleadpho_fwdjet",  &_dr_subleadpho_fwdjet);
        TTHMva_->AddVariable( "dr_bjet_fwdjet",  &_dr_bjet_fwdjet);
        TTHMva_->BookMVA( "BDT", TTHMVAweightfile_.fullPath() );

        BkgMva_.reset( new TMVA::Reader( "!Color:!Silent" ) );
        BkgMva_->AddVariable( "lead_ptoM", &_lead_ptoM );
        BkgMva_->AddVariable( "sublead_ptoM", &_sublead_ptoM );
        BkgMva_->AddVariable( "leadPixelSeed", &_leadPixelSeed );
        BkgMva_->AddVariable( "subleadPixelSeed", &_subleadPixelSeed );
        BkgMva_->AddVariable( "bjet1_csv", &_bjet1_csv );
        BkgMva_->AddVariable( "bjet2_csv", &_bjet2_csv );
        BkgMva_->AddVariable( "n_jets", &_n_jets );
        BkgMva_->AddVariable( "dipho_ptoM", &_dipho_ptoM );
        BkgMva_->AddVariable( "deltaR_gg", &_deltaR_gg );
        BkgMva_->AddVariable( "dipho_rapidity", &_dipho_rapidity );
        BkgMva_->AddVariable( "jet1_bdiscriminant", &_jet1_bdiscriminant );
        BkgMva_->AddVariable( "jet2_bdiscriminant", &_jet2_bdiscriminant );
        BkgMva_->AddVariable( "jet3_bdiscriminant", &_jet3_bdiscriminant );
        BkgMva_->AddVariable( "jet1_pt", &_jet1_pt );
        BkgMva_->AddVariable( "jet2_pt", &_jet2_pt );
        BkgMva_->AddVariable( "jet3_pt", &_jet3_pt );
        BkgMva_->AddVariable( "jet1_eta", &_jet1_eta );
        BkgMva_->AddVariable( "jet2_eta", &_jet2_eta );
        BkgMva_->AddVariable( "jet3_eta", &_jet3_eta );
        BkgMva_->AddVariable( "dipho_helicity", &_dipho_helicity );
        BkgMva_->AddVariable( "min_id", &_min_id );
        BkgMva_->AddVariable( "max_id", &_max_id );
        BkgMva_->AddVariable( "nbjets", &_nbjets );
        BkgMva_->AddVariable( "ncentraljets", &_ncentraljets );
        BkgMva_->AddVariable( "bjet1_pt", &_bjet1_pt );
        BkgMva_->AddVariable( "bjet2_pt", &_bjet2_pt );
        BkgMva_->AddVariable( "bjet3_pt", &_bjet3_pt );
        BkgMva_->AddVariable( "bjet1_eta", &_bjet1_eta );
        BkgMva_->AddVariable( "bjet2_eta", &_bjet2_eta );
        BkgMva_->AddVariable( "bjet3_eta", &_bjet3_eta );
        BkgMva_->AddVariable( "bjet3_csv", &_bjet3_csv );
        BkgMva_->AddVariable( "dr_leadpho_fwdjet", &_dr_leadpho_fwdjet );
        BkgMva_->AddVariable( "dr_subleadpho_fwdjet", &_dr_subleadpho_fwdjet );
        BkgMva_->AddVariable( "dr_bjet_fwdjet", &_dr_bjet_fwdjet );
        BkgMva_->BookMVA( "BDT", BkgMVAweightfile_.fullPath() );

        produces<vector<THQHadronicTag> >();
    }

    void THQHadronicTagProducer::produce( Event &evt, const EventSetup & )
    {
        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::Met> > METs;
        evt.getByToken( METToken_, METs );

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        edm::Handle<double>  rho;
        evt.getByToken(rhoTag_,rho);
        double rho_    = *rho;

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );
        assert( diPhotons->size() == mvaResults->size() );

        Handle<View<reco::GenParticle> > genParticles;
        evt.getByToken( genParticleToken_, genParticles );

        std::unique_ptr<vector<THQHadronicTag> > thqHadronicTags( new vector<THQHadronicTag> );

        double idmva1 = 0.;
        double idmva2 = 0.;

        bool hasGoodElec  = false;
        bool hasGoodMuons = false;
        
        for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
            edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

            THQHadronicTag thqHadronicTags_obj( dipho, mvares );
            thqHadronicTags_obj.includeWeights( *dipho );

            idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            
            if( dipho->leadingPhoton()->pt()    < dipho->mass() * leadPhoOverMassThreshold_ ) continue;
            if( dipho->subLeadingPhoton()->pt() < dipho->mass() * subleadPhoOverMassThreshold_ ) continue;
            if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) continue;
            if( mvares->result < MVAThreshold_ ) continue;

            // Lepton
            std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons(  theMuons->ptrs(), 
                                                                            dipho, 
                                                                            vertices->ptrs(), 
                                                                            muonEtaThreshold_, 
                                                                            leptonPtThreshold_,
                                                                            muPFIsoSumRelThreshold_, 
                                                                            deltaRMuonPhoThreshold_, 
                                                                            deltaRMuonPhoThreshold_ 
                                                                            );
            
            std::vector<edm::Ptr<Electron> >goodElectrons = selectStdElectrons( theElectrons->ptrs(), 
                                                                                dipho,
                                                                                vertices->ptrs(), 
                                                                                leptonPtThreshold_, 
                                                                                electronEtaThresholds_,
                                                                                useElectronMVARecipe_,
                                                                                useElectronLooseID_,
                                                                                deltaRPhoElectronThreshold_,
                                                                                DeltaRTrkElec_,
                                                                                deltaMassElectronZThreshold_,
                                                                                rho_, 
                                                                                evt.isRealData() 
                                                                                );

            hasGoodElec  = ( goodElectrons.size() > 0 );
            hasGoodMuons = ( goodMuons.size() > 0 );
                        
            if( !hasGoodElec && !hasGoodMuons ) { continue; }
            //including SFs for leading muon or electron 
            if( goodMuons.size() > 0 ) {
                thqHadronicTags_obj.includeWeights( *goodMuons.at(0) );
            } else if ( goodElectrons.size() > 0 ) {
                thqHadronicTags_obj.includeWeights( *goodElectrons.at(0) );
            }

            //Jets
            std::vector<edm::Ptr<Jet> > tagJets;
            std::vector<edm::Ptr<Jet> > bJets;
            std::vector<edm::Ptr<Jet> > forwardJets;
            unsigned nCentralJets = 0;
            unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();
            for( unsigned int candIndex_outer = 0; candIndex_outer < Jets[jetCollectionIndex]->size() ; candIndex_outer++ ) {

                bool keepJet = true;
                edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );
                if(!thejet->passesJetID  ( flashgg::Tight2017 ) ) { continue; }
                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { keepJet = false; }
                if( thejet->pt() < jetPtThreshold_ ) { keepJet = false; }
                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                dipho->subLeadingPhoton()->superCluster()->phi() );
                
                if( dRPhoLeadJet < deltaRPhoLeadJet_ || dRPhoSubLeadJet < deltaRPhoSubLeadJet_ ) { keepJet = false; }
                if( hasGoodElec ) 
                    for( unsigned int electronIndex = 0; electronIndex < goodElectrons.size(); electronIndex++ ) {
                            Ptr<flashgg::Electron> electron = goodElectrons[electronIndex];
                            float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), electron->eta(), electron->phi() ) ;
                            if( dRJetElectron < deltaRJetMuonThreshold_ ) { keepJet = false; }
                    }
                if( hasGoodMuons ) 
                    for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) {
                            Ptr<flashgg::Muon> muon = goodMuons[muonIndex];
                            float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;
                            if( dRJetMuon < deltaRJetMuonThreshold_ ) { keepJet = false; }
                    }
                if(keepJet) { 
                    tagJets.push_back( thejet );
                    if (thejet->bDiscriminator("pfDeepCSVJetTags:probb") + thejet->bDiscriminator("pfDeepCSVJetTags:probbb") > bDiscriminator_[0]) bJets.push_back( thejet );
                    if (forwardJets.size()==0) forwardJets.push_back(thejet);
                    else if (abs(thejet->eta()) > abs(forwardJets[0]->eta())) {
                        forwardJets.clear();
                        forwardJets.push_back(thejet);
                    }
                    if (abs(thejet->eta() < 1.)) nCentralJets += 1;
                }
            }

            if (tagJets.size() > jetsNumberThreshold_) continue;

            //MET
            Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
            if ( theMET->getCorPt() < METThreshold_ ) continue;

            //MVAs
            //double ele_pt = -999.;
            //double mu_pt  = -999.;

            //if (goodElectrons.size() > 0) ele_pt = goodElectrons[0]->pt();
            //if (goodMuons.size() > 0)     mu_pt  = goodMuons[0]->pt();

            TLorentzVector diphoP4(dipho->px(), dipho->py(), dipho->pz(), dipho->energy());

            _lead_ptoM        = dipho->leadingPhoton()->eta();
            _sublead_ptoM     = dipho->subLeadingPhoton()->eta();
            _leadEta          = dipho->leadingPhoton()->phi();
            _subleadEta       = dipho->subLeadingPhoton()->phi();
            _leadPhi          = dipho->leadingPhoton()->pt() / dipho->mass();
            _subleadPhi       = dipho->subLeadingPhoton()->pt() / dipho->mass();
            _min_id           = TMath::Min(idmva1, idmva2);
            _max_id           = TMath::Max(idmva1, idmva2);
            _leadPixelSeed    = dipho->leadingPhoton()->hasPixelSeed() > 0.5 ? 0. : 1.;
            _subleadPixelSeed = dipho->subLeadingPhoton()->hasPixelSeed() > 0.5 ? 0. : 1.;

            _deltaR_gg      = deltaR( dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi());
            TLorentzVector pho1, pho2;
            pho1.SetPtEtaPhiE(dipho->leadingPhoton()->pt(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), dipho->leadingPhoton()->energy());
            pho2.SetPtEtaPhiE(dipho->subLeadingPhoton()->pt(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi(), dipho->subLeadingPhoton()->energy());
            _dipho_helicity = helicity(pho1, pho2);
            _dipho_cosphi   = TMath::Cos( deltaPhi(dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->phi()) );
            _dipho_ptoM     = dipho->pt()/dipho->mass();
            _dipho_rapidity = dipho->rapidity();

            _n_jets = tagJets.size();
            if (tagJets.size() > 0) {
                _jet1_bdiscriminant = tagJets[0]->bDiscriminator("pfDeepCSVJetTags:probb") + tagJets[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _jet1_pt = tagJets[0]->pt();
                _jet1_eta = tagJets[0]->eta();
            }
            if (tagJets.size() > 1) {
                _jet2_bdiscriminant = tagJets[1]->bDiscriminator("pfDeepCSVJetTags:probb") + tagJets[1]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _jet2_pt = tagJets[1]->pt();
                _jet2_eta = tagJets[1]->eta();
            }
            if (tagJets.size() > 2) {
                _jet3_bdiscriminant = tagJets[2]->bDiscriminator("pfDeepCSVJetTags:probb") + tagJets[2]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _jet3_pt = tagJets[2]->pt();
                _jet3_eta = tagJets[2]->eta();
            }

            _nbjets = bJets.size();
            if (bJets.size() > 0) {
                _bjet1_csv = bJets[0]->bDiscriminator("pfDeepCSVJetTags:probb") + bJets[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _bjet1_pt  = bJets[0]->pt();
                _bjet1_eta = bJets[0]->eta();
            }
            if (bJets.size() > 1) {
                _bjet2_csv = bJets[1]->bDiscriminator("pfDeepCSVJetTags:probb") + bJets[1]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _bjet2_pt  = bJets[1]->pt();
                _bjet2_eta = bJets[1]->eta();
            }
            if (bJets.size() > 2) {
                _bjet3_csv = bJets[2]->bDiscriminator("pfDeepCSVJetTags:probb") + bJets[2]->bDiscriminator("pfDeepCSVJetTags:probbb");
                _bjet3_pt  = bJets[2]->pt();
                _bjet3_eta = bJets[2]->eta();
            }

            _ncentraljets = nCentralJets;
            if (forwardJets.size()>0) {
                _forward_pt = forwardJets[0]->pt();
                _forward_eta = forwardJets[0]->eta();
                _dr_leadpho_fwdjet = deltaR( dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), forwardJets[0]->eta(), forwardJets[0]->phi() );
                _dr_subleadpho_fwdjet = deltaR( dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi(), forwardJets[0]->eta(), forwardJets[0]->phi() );
                if (bJets.size() > 0) _dr_bjet_fwdjet = deltaR( bJets[0]->eta(), bJets[0]->phi(), forwardJets[0]->eta(), forwardJets[0]->phi() );
            }

            float tthMva = TTHMva_->EvaluateMVA( "BDT" );
            float tthMvaTform = 1. / ( 1. + exp( 0.5*log( 2./(tthMva+1.) - 1 ) ) );
            float bkgMva = BkgMva_->EvaluateMVA( "BDT" );
            float bkgMvaTform = 1. / ( 1. + exp( 0.5*log( 2./(bkgMva+1.) - 1 ) ) );

            // Categorization with both MVAs
            bool isSelected = (tthMvaTform > tthMvaThreshold_ && bkgMvaTform > bkgMvaThreshold_);

            if( isSelected ) {
                thqHadronicTags_obj.setStage1recoTag( DiPhotonTagBase::stage1recoTag::RECO_THQ_HAD );
                thqHadronicTags_obj.setJets( tagJets );
                thqHadronicTags_obj.setMuons( goodMuons );
                thqHadronicTags_obj.setElectrons( goodElectrons );
                thqHadronicTags_obj.setDiPhotonIndex( diphoIndex );
                thqHadronicTags_obj.setSystLabel( systLabel_ );
                thqHadronicTags_obj.setMET( theMET );
                thqHadronicTags->push_back( thqHadronicTags_obj );
            }
        }
        evt.put( std::move( thqHadronicTags ) );
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

