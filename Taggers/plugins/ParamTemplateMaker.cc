// ---------------------------------------------------------------------
// E.Scott    08/2016
//
// Makes initial templates for prompt-fake parameterisation
// ---------------------------------------------------------------------

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"

#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/VBFMVAResult.h"
#include "flashgg/DataFormats/interface/VBFDiPhoDiJetMVAResult.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CLHEP/Random/RandFlat.h"


using namespace std;
using namespace edm;
using namespace flashgg;

class ParamTemplateMaker : public edm::EDAnalyzer
{
public:
    explicit ParamTemplateMaker( const edm::ParameterSet & );
    ~ParamTemplateMaker();

    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


private:

    edm::Service<TFileService> fs_;

    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    EDGetTokenT< edm::View<reco::GenJet> >               genJetToken_;
    EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    EDGetTokenT< View<reco::Vertex> >                    vertexToken_;
    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;

    bool        debug_;

    TH1F* hFakeGenJetRatio;
    TH1F* hBarrelLowTemplateIDMVA;
    TH1F* hBarrelHighTemplateIDMVA;
    TH1F* hEndcapLowTemplateIDMVA;
    TH1F* hEndcapHighTemplateIDMVA;
    TH1F* hSigmarvLowEB;
    TH1F* hSigmarvHighEB;
    TH1F* hSigmarvLowEE;
    TH1F* hSigmarvHighEE;
    TH1F* hFakeVtxprob;
    //TH1F* hCorrectPt;
    //TH1F* hWrongPt;
};

ParamTemplateMaker::ParamTemplateMaker( const edm::ParameterSet &iConfig ):
    genJetToken_( consumes<View<reco::GenJet> >( iConfig.getUntrackedParameter<InputTag> ( "GenJetTag", InputTag( "slimmedGenJets" ) ) ) ),
    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    debug_( iConfig.getUntrackedParameter<bool>( "debug", false ) )
{
    cout << "Inside constructor of ParamTemplateMaker" << endl;

    hFakeGenJetRatio = fs_->make<TH1F>( "hFakeGenJetRatio","fakePhoton energy / energy nearest genJet ",48,0.,1.2 );
    hBarrelLowTemplateIDMVA = fs_->make<TH1F>( "hBarrelLowTemplateIDMVA", "IDMVA for low FakeGenJetRatio, in barrel", 20,-1.,1. );
    hBarrelHighTemplateIDMVA = fs_->make<TH1F>( "hBarrelHighTemplateIDMVA","IDMVA for high FakeGenJetRatio, in barrel",20,-1.,1. );
    hEndcapLowTemplateIDMVA = fs_->make<TH1F>( "hEndcapLowTemplateIDMVA", "IDMVA for low FakeGenJetRatio, in endcap", 20,-1.,1. );
    hEndcapHighTemplateIDMVA = fs_->make<TH1F>( "hEndcapHighTemplateIDMVA","IDMVA for high FakeGenJetRatio, in endcap",20,-1.,1. );
    hSigmarvLowEB = fs_->make<TH1F>( "hSigmarvLowEB", "Sigmarv for low FakeGenJetRatio, in barrel", 20, 0., 0.05 );
    hSigmarvHighEB = fs_->make<TH1F>( "hSigmarvHighEB", "Sigmarv for high FakeGenJetRatio, in barrel", 20, 0., 0.05 );
    hSigmarvLowEE = fs_->make<TH1F>( "hSigmarvLowER", "Sigmarv for low FakeGenJetRatio, in endcap", 20, 0., 0.05 );
    hSigmarvHighEE = fs_->make<TH1F>( "hSigmarvHighEE", "Sigmarv for high FakeGenJetRatio, in endcap", 20, 0., 0.05 );
    hFakeVtxprob = fs_->make<TH1F>( "hFakeVtxprob", "Fake vtxprob", 200, 0., 1.0 );
    
    // this can't be done here; need to think exactly how to implement elegantly
    //hCorrectPt = fs_->make<TH1F>(  );
    //hWrongPt = fs_->make<TH1F>(  );
}

ParamTemplateMaker::~ParamTemplateMaker()
{
}


void
ParamTemplateMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{
    if( debug_ ) cout << "Inside ParamTemplateMaker analyze method" << endl;

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );

    Handle<View<reco::GenJet> > genJets;
    iEvent.getByToken( genJetToken_, genJets );

    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
    iEvent.getByToken( mvaResultToken_, mvaResults );
    
    size_t diPhotonsSize = diPhotons->size();

    //unsigned int actualDiphoIndex = randomDiphoIndex->Integer( diPhotonsSize );
    //unsigned int actualDiphoIndex = CLHEP::RandFlat::shootInt( &engine, diPhotonsSize ); // to be parameterised

    //cout << "actualDiphoIndex = " << actualDiphoIndex << endl;
    //cout << "diPhotonsSize = " << diPhotonsSize << endl;

    for( unsigned int diphoIndex = 0; diphoIndex < diPhotonsSize; diphoIndex++ ) {

        edm::Ptr<flashgg::DiPhotonMVAResult> dipho_mvares = mvaResults->ptrAt( diphoIndex );

        //--------------------------------------------------------------------------------------------------
        // Begin analysis of prompt-fake events // Ed
        //--------------------------------------------------------------------------------------------------
        auto printDipho = diPhotons->ptrAt( diphoIndex );
        // need to figure out how to add this properly
        //float weight = printDipho->weight();
        float weight = 1.;

        if( printDipho->mass() < 100 || printDipho->mass() > 180 ) continue;
        /*if( genJets->size()    < 2 ) continue;
        if( genPhotons->size() < 2 ) continue;
        if( gens->size()       < 2 ) continue;*/

        auto printLeadPho    = printDipho->leadingPhoton();
        auto printSubLeadPho = printDipho->subLeadingPhoton();
        flashgg::Photon::mcMatch_t leadMatchType    = printLeadPho->genMatchType();
        //cout << "lead match type is " << leadMatchType << endl;
        flashgg::Photon::mcMatch_t subleadMatchType = printSubLeadPho->genMatchType();
        //cout << "sublead match type is " << subleadMatchType << endl;
        bool eventIsPromptFake = ( (!(leadMatchType==1 && subleadMatchType==1)) && (leadMatchType==1 || subleadMatchType==1) );
        auto promptPhoton      = printLeadPho;
        auto fakePhoton        = printSubLeadPho;
        //bool reversed          = false;
        //float promptIDMVA      = printDipho->leadPhotonId();
        float fakeIDMVA        = printDipho->subLeadPhotonId();
        if( leadMatchType != 1 ) {
            promptPhoton = printSubLeadPho;
            fakePhoton   = printLeadPho;
            //reversed     = true;
            //promptIDMVA  = printDipho->subLeadPhotonId();
            fakeIDMVA    = printDipho->leadPhotonId();
        }
        
        auto  genJetNearestPrompt     = genJets->ptrAt(0);
        auto  genJetNearestFake       = genJets->ptrAt(1);
        int   promptNumGenJets        = 0;
        int   fakeNumGenJets          = 0;

        if( eventIsPromptFake ) {
            //cout << "printDipho nConv = " << printDipho->nConv() << endl;
            //cout << "promptPhoton has conv tracks? = " << promptPhoton->hasConversionTracks() << endl;
            //cout << "fakePhoton   has conv tracks? = " << fakePhoton->hasConversionTracks()   << endl << endl;
            float promptEta   = promptPhoton->eta();
            //float promptEta   = promptPhoton->superCluster()->eta();
            float promptPhi   = promptPhoton->phi();
            float fakeEta     = fakePhoton->eta();
            //float fakeEta     = fakePhoton->superCluster()->eta();
            float fakePhi     = fakePhoton->phi();
            float minDrPrompt = 10.0;
            float minDrFake   = 10.0;

            // Print nearest genJet and whether it is within dR of 0.4
            minDrPrompt = 10.0;
            minDrFake   = 10.0;
            for( uint genJetIndex = 0; genJetIndex < genJets->size(); genJetIndex++ ) {
                float tempEta  = genJets->ptrAt( genJetIndex )->eta();
                float tempPhi  = genJets->ptrAt( genJetIndex )->phi();
                float drPrompt = deltaR( tempEta, tempPhi, promptEta, promptPhi );
                float drFake   = deltaR( tempEta, tempPhi, fakeEta,   fakePhi   );
                //float drPrompt = sqrt( (tempEta-promptEta)*(tempEta-promptEta) + (tempPhi-promptPhi)*(tempPhi-promptPhi) );
                //float drFake   = sqrt( (tempEta-fakeEta)*(tempEta-fakeEta)     + (tempPhi-fakePhi)*(tempPhi-fakePhi)     );
                if( drPrompt < minDrPrompt ) {
                    minDrPrompt = drPrompt;
                    genJetNearestPrompt = genJets->ptrAt( genJetIndex );
                    if ( drPrompt < 0.4 ) promptNumGenJets++;
                }
                else if( drPrompt < 0.4 ) promptNumGenJets++;
                if( drFake < minDrFake ) {
                    minDrFake = drFake;
                    genJetNearestFake = genJets->ptrAt( genJetIndex );
                    if ( drFake < 0.4 ) fakeNumGenJets++;
                }
                else if( drFake < 0.4 ) fakeNumGenJets++;
            } // end of genJet matching

            float fakeGenJetRatio = fakePhoton->energy() / genJetNearestFake->energy();
            hFakeGenJetRatio->Fill( fakeGenJetRatio, weight );
            if(      abs( fakeEta ) < 1.5 && fakeGenJetRatio < 0.8 ) hBarrelLowTemplateIDMVA->Fill(  fakeIDMVA, weight );
            else if( abs( fakeEta ) < 1.5 && fakeGenJetRatio < 1.2 ) hBarrelHighTemplateIDMVA->Fill( fakeIDMVA, weight );
            else if( abs( fakeEta ) < 2.5 && fakeGenJetRatio < 0.8 ) hEndcapLowTemplateIDMVA->Fill(  fakeIDMVA, weight );
            else if( abs( fakeEta ) < 2.5 && fakeGenJetRatio < 1.2 ) hEndcapHighTemplateIDMVA->Fill( fakeIDMVA, weight );

            float sigmarv = dipho_mvares->sigmarv;
            if(      abs( fakeEta ) < 1.5 && fakeGenJetRatio < 0.8 ) hSigmarvLowEB->Fill(  sigmarv, weight );
            else if( abs( fakeEta ) < 1.5 && fakeGenJetRatio < 1.2 ) hSigmarvHighEB->Fill( sigmarv, weight );
            else if( abs( fakeEta ) < 2.5 && fakeGenJetRatio < 0.8 ) hSigmarvLowEE->Fill(  sigmarv, weight );
            else if( abs( fakeEta ) < 2.5 && fakeGenJetRatio < 1.2 ) hSigmarvHighEE->Fill( sigmarv, weight );

            float vtxprob = dipho_mvares->vtxprob;
            hFakeVtxprob->Fill( vtxprob, weight );

        //--------------------------------------------------------------------------------------------------
        } // End of prompt-fake events */ // Ed
        //--------------------------------------------------------------------------------------------------

    } //  end of loop over the diphotons
    //cout << "Exiting JetValTreeMaker analyze method" << endl;
}

void
ParamTemplateMaker::beginJob()
{
}

void ParamTemplateMaker::endJob()
{
}

void ParamTemplateMaker::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}

typedef ParamTemplateMaker FlashggParamTemplateMaker;
DEFINE_FWK_MODULE( FlashggParamTemplateMaker );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
