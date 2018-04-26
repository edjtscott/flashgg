#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

using namespace flashgg;

DiPhotonTagBase::DiPhotonTagBase()
{
    category_number_ = -1;
    isGold_ = -1;
}

DiPhotonTagBase::DiPhotonTagBase( edm::Ptr<flashgg::DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvaRes )
    : DiPhotonTagBase::DiPhotonTagBase( diPho, *mvaRes ) {}

DiPhotonTagBase::DiPhotonTagBase( edm::Ptr<flashgg::DiPhotonCandidate> diPho, DiPhotonMVAResult mvaRes )
{
    mva_result_ = mvaRes;
    category_number_ = -1;
    dipho_ = diPho;
    isGold_ = -1;
    stage1KinematicLabel_ = "";
}

DiPhotonTagBase::~DiPhotonTagBase()
{
}


bool DiPhotonTagBase::operator <( const DiPhotonTagBase &b ) const
{
    // For choosing which of two tags OF THE SAME TYPE are of higher priority
    // Comparison of different tags not currently supported - is it ever needed?
    // Overloading may be appropriate if different tags have different priorities

    if( categoryNumber() == b.categoryNumber() ) {
        return ( sumPt() < b.sumPt() );
    } else {
        return ( categoryNumber() < b.categoryNumber() );
    }
}

void DiPhotonTagBase::setIsGold( int runNumber ) {
    // Below is the subtraction between
    // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt
    // and
    // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
    //
    // They were last changed on 18 December and still the most recent (as of 28 January)
    // I have confirmed that none of these runs are in the Gold at all, so checking runNumber without lumiSection suffices
    // Of course the result may be wrong for events not even in the Silver JSON
    // I am hard coding this to make dumping and comparing events faster, because this issue now needs intensive study
    //
    /*
      {"256729": [[1, 331], [346, 598], [600, 755], [758, 760], [765, 1165], [1167, 1292], [1295, 1327], [1329, 1732]],
 "256734": [[1, 57], [60, 213]],
 "257394": [[41, 72]],
 "257395": [[1, 13]],
 "257396": [[1, 216]],
 "257397": [[1, 119]],
 "257399": [[1, 271]],
 "257400": [[1, 291], [295, 819], [1011, 1418]],
 "257487": [[50, 102], [104, 202], [204, 1124]],
 "257490": [[1, 591]],
 "257822": [[1, 719], [721, 1389]],
 "257823": [[1, 171]],
 "258443": [[1, 291]]}
    */
    isGold_ = 1;
    if ( runNumber == 256729 ) { isGold_ = 0; }
    if ( runNumber == 256734 ) { isGold_ = 0; }
    if ( (runNumber >= 257394) && (runNumber <= 257397) ) { isGold_ = 0; }
    if ( (runNumber >= 257399) && (runNumber <=257400)) { isGold_ = 0; }
    if ( runNumber == 257487 ) { isGold_ = 0; }
    if ( runNumber == 257490 ) { isGold_ = 0; }
    if ( (runNumber >= 257822) && (runNumber <=257823)) { isGold_ = 0; }
    if ( runNumber == 258443 ) { isGold_ = 0; }
}

float DiPhotonTagBase::ggHweightCentralised( std::string weightName ) const {
    float weightCentral  = centralWeight();
    float weightOriginal = tagTruth()->HTXSggHweight( weightName );
    return weightCentral * weightOriginal;
}

void DiPhotonTagBase::computeStage1Kinematics( const edm::Handle<edm::View<flashgg::Jet> >& jets, float ptV, float lepeta1, float lepphi1, float lepeta2, float lepphi2 ) {
    stage1KinematicLabel_ = "LOGICERROR";
    stage1recoTag_ = DiPhotonTagBase::LOGICERROR;
    float ptH = this->diPhoton()->pt();
    unsigned int nJ = 0;
    float dEta = 0.;
    float mjj = 0.;
    float ptHjj = 0.;
    float mvaScore = this->diPhotonMVA().mvaValue();
    edm::Ptr<flashgg::Jet> j0;
    edm::Ptr<flashgg::Jet> j1;
    for ( unsigned int i = 0 ; i < jets->size(); i++ ) {
        edm::Ptr<flashgg::Jet> jet = jets->ptrAt(i);

        //        std::cout << " Jet " << i << " pt=" << jet->pt() << " eta=" << jet->eta() << std::endl;

        if ( jet->pt() < 30. ) continue;
        if ( fabs(jet->eta()) > 4.7 ) continue;
        bool _useJetID = true; // temporary
        std::string _JetIDLevel = "Tight"; // temporary
        float _drJetPhoton = 0.4; // Temporary
        float _drJetLepton = 0.4; // Temporary
        // Temporary: No PU JET ID
        if( _useJetID ){
            if( _JetIDLevel == "Loose" && !jet->passesJetID  ( flashgg::Loose ) ) continue;
            if( _JetIDLevel == "Tight" && !jet->passesJetID  ( flashgg::Tight ) ) continue;
        }

        // close to lead photon?                                                                                                                                                                                                                          
        float dPhi = deltaPhi( jet->phi(), this->diPhoton()->leadingPhoton()->phi() );
        float dEta = jet->eta() - this->diPhoton()->leadingPhoton()->eta();
        if( sqrt( dPhi * dPhi + dEta * dEta ) < _drJetPhoton ) { continue; }

        // close to sublead photon?                                                                                                                                                                                                                       
        dPhi = deltaPhi( jet->phi(), this->diPhoton()->subLeadingPhoton()->phi() );
        dEta = jet->eta() - this->diPhoton()->subLeadingPhoton()->eta();
        if( sqrt( dPhi * dPhi + dEta * dEta ) < _drJetPhoton ) { continue; }

        // close to lepton1 (if any)
        dPhi = deltaPhi( jet->phi(), lepphi1 );
        dEta = jet->eta() - lepeta1;
        if( sqrt( dPhi * dPhi + dEta * dEta ) < _drJetLepton ) { continue; }

        // close to lepton2 (if any)
        dPhi = deltaPhi( jet->phi(), lepphi2 );
        dEta = jet->eta() - lepeta2;
        if( sqrt( dPhi * dPhi + dEta * dEta ) < _drJetLepton ) { continue; }

        nJ++;

        if ( j0.isNull() ) {
            //            std::cout << " Save jet " << i << " as j0" << std::endl;
            j0 = jet;
        } else if ( j1.isNull() ) {
            //            std::cout << " Save jet " << i << " as j1" << std::endl;
            j1 = jet;
        } else { 
            //            std::cout << " Not saving jet " << i << " - two jets already " << std::endl;
        }
    }
    //    std::cout << " nJ=" << nJ << " ptV=" << ptV << " ptH=" << ptH << std::endl;

    unsigned nlep = 0;
    if (lepphi1 > -998. ) nlep++;
    if (lepphi2 > -998. ) nlep++;
    string nlepstring = std::to_string(nlep)+"LEP";

    if ( nJ >= 2 ) {
        dEta = fabs( j0->eta() - j1->eta() );
        mjj = ( j0->p4() + j1->p4() ).mass();
        ptHjj = ( j0->p4() + j1->p4() + this->diPhoton()->p4() ).pt();
        //        std::cout << " dEta=" << dEta << " mjj=" << mjj << " ptHjj=" << ptHjj << std::endl;
    }
    // have now added two categories for each RECO tag, using the moment diphoton MVA, with boundaries currently hard-coded below..
    if ( ptV < -0.5 ) {
        if (nJ == 0) {
            if (mvaScore > 0.46) {
                stage1KinematicLabel_ = "RECO_0J_Tag0";
                stage1recoTag_ = DiPhotonTagBase::RECO_0J_Tag0;
            }
            else if (mvaScore > -0.2) {
                stage1KinematicLabel_ = "RECO_0J_Tag1";
                stage1recoTag_ = DiPhotonTagBase::RECO_0J_Tag1;
            }
            else { 
                stage1KinematicLabel_ = "NOTAG";
                stage1recoTag_ = DiPhotonTagBase::NOTAG;
            }
        } else if ( nJ == 1 ) {
            if ( ptH > 200 ) {
                if (mvaScore > 0.86) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_GT200";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_GT200;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else if ( ptH > 120. ) {
                if (mvaScore > 0.8) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_120_200_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_120_200_Tag0;
                }
                else if (mvaScore > 0.4) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_120_200_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_120_200_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else if ( ptH > 60. ) {
                if (mvaScore > 0.64) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_60_120_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_60_120_Tag0;
                }
                else if (mvaScore > 0.1) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_60_120_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_60_120_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else {
                if (mvaScore > 0.56) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_0_60_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_0_60_Tag0;
                }
                else if (mvaScore > -0.05) {
                    stage1KinematicLabel_ = "RECO_1J_PTH_0_60_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_1J_PTH_0_60_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            }
        } else { // 2 jets
            if ( ptH > 200 ) {
                if (mvaScore > 0.92) {
                    stage1KinematicLabel_ ="RECO_GE2J_PTH_GT200_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_GT200_Tag0;
                }
                else if (mvaScore > 0.8) {
                    stage1KinematicLabel_ ="RECO_GE2J_PTH_GT200_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_GT200_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else if ( mjj > 400. && dEta > 2.8 ) {
                if ( ptHjj < 25. ) {
                    if (mvaScore > 0.52) {
                        stage1KinematicLabel_ = "RECO_VBFTOPO_JET3VETO_Tag0";
                        stage1recoTag_ = DiPhotonTagBase::RECO_VBFTOPO_JET3VETO_Tag0;
                    }
                    else if (mvaScore > 0.) {
                        stage1KinematicLabel_ = "RECO_VBFTOPO_JET3VETO_Tag1";
                        stage1recoTag_ = DiPhotonTagBase::RECO_VBFTOPO_JET3VETO_Tag1;
                    }
                    else { 
                        stage1KinematicLabel_ = "NOTAG";
                        stage1recoTag_ = DiPhotonTagBase::NOTAG;
                    }
                } else {
                    if (mvaScore > 0.72) {
                        stage1KinematicLabel_ = "RECO_VBFTOPO_JET3_Tag0";
                        stage1recoTag_ = DiPhotonTagBase::RECO_VBFTOPO_JET3_Tag0;
                    }
                    else if (mvaScore > 0.3) {
                        stage1KinematicLabel_ = "RECO_VBFTOPO_JET3_Tag1";
                        stage1recoTag_ = DiPhotonTagBase::RECO_VBFTOPO_JET3_Tag1;
                    }
                    else { 
                        stage1KinematicLabel_ = "NOTAG";
                        stage1recoTag_ = DiPhotonTagBase::NOTAG;
                    }
                }
            } else if ( ptH > 120. ) {
                if (mvaScore > 0.86) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_120_200_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_120_200_Tag0;
                }
                else if (mvaScore > 0.5) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_120_200_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_120_200_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else if ( ptH > 60. ) {
                if (mvaScore > 0.75) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_60_120_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_60_120_Tag0;
                }
                else if (mvaScore > 0.3) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_60_120_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_60_120_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            } else {
                if (mvaScore > 0.56) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_0_60_Tag0";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_0_60_Tag0;
                }
                else if (mvaScore > -0.1) {
                    stage1KinematicLabel_ = "RECO_GE2J_PTH_0_60_Tag1";
                    stage1recoTag_ = DiPhotonTagBase::RECO_GE2J_PTH_0_60_Tag1;
                }
                else { 
                    stage1KinematicLabel_ = "NOTAG";
                    stage1recoTag_ = DiPhotonTagBase::NOTAG;
                }
            }
        }
    } else { // Leptonic vector boson assigned. Leave this up to existing VH tags for now
        stage1KinematicLabel_ = "NOTAG";
        stage1recoTag_ = DiPhotonTagBase::NOTAG;
    }
    //    std::cout << " Final label " << stage1KinematicLabel_ << std::endl;
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

