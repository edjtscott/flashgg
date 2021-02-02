#include "flashgg/DataFormats/interface/THQHadronicTag.h"
#include <algorithm>

using namespace flashgg;

THQHadronicTag::THQHadronicTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

THQHadronicTag::~THQHadronicTag()
{}


THQHadronicTag::THQHadronicTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : DiPhotonTagBase::DiPhotonTagBase( diPho, *mvares ) {}
THQHadronicTag::THQHadronicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

