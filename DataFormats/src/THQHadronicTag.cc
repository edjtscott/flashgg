#include "flashgg/DataFormats/interface/THQHadronicTag.h"
#include "flashgg/DataFormats/interface/Jet.h"

using namespace flashgg;

THQHadronicTag::THQHadronicTag() : DiPhotonTagBase::DiPhotonTagBase() {}

THQHadronicTag::~THQHadronicTag() {}

THQHadronicTag::THQHadronicTag( edm::Ptr<flashgg::DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvaRes,
                                std::vector<edm::Ptr<flashgg::Jet> > theJetVec , std::vector<edm::Ptr<flashgg::Jet> > theBJetVec ) :
    THQHadronicTag::THQHadronicTag( diPho, *mvaRes, theJetVec, theBJetVec ) {}

THQHadronicTag::THQHadronicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares,
                                std::vector<edm::Ptr<flashgg::Jet> > theJetVec , std::vector<edm::Ptr<flashgg::Jet> > theBJetVec ) :
    DiPhotonTagBase::DiPhotonTagBase( dipho, mvares )
{
    theJetVec_ = std::vector<edm::Ptr<flashgg::Jet> >( theJetVec );
    theBJetVec_ = std::vector<edm::Ptr<flashgg::Jet> >( theBJetVec );
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

