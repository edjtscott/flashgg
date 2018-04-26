#ifndef FLASHgg_DiPhotonTagBase_h
#define FLASHgg_DiPhotonTagBase_h

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"

namespace flashgg {

    class DiPhotonTagBase : public WeightedObject
    {
    public:
        enum tag_t { kUndefined = 0, kUntagged, kVBF, kTTHHadronic, kTTHLeptonic, kVHTight, kVHLoose, kVHHadronic, kVHEt,  kZHLeptonic, kWHLeptonic, kVHLeptonicLoose, kVHMet };
        enum stage1recoTag { LOGICERROR = -1, NOTAG = 0, RECO_0J_Tag0, RECO_0J_Tag1, RECO_1J_PTH_0_60_Tag0, RECO_1J_PTH_0_60_Tag1, RECO_1J_PTH_60_120_Tag0, RECO_1J_PTH_60_120_Tag1,
                             RECO_1J_PTH_120_200_Tag0, RECO_1J_PTH_120_200_Tag1, RECO_1J_PTH_GT200, 
                             RECO_GE2J_PTH_0_60_Tag0, RECO_GE2J_PTH_0_60_Tag1, RECO_GE2J_PTH_60_120_Tag0, RECO_GE2J_PTH_60_120_Tag1, RECO_GE2J_PTH_120_200_Tag0, RECO_GE2J_PTH_120_200_Tag1, 
                             RECO_GE2J_PTH_GT200_Tag0, RECO_GE2J_PTH_GT200_Tag1, RECO_VBFTOPO_JET3VETO_Tag0, RECO_VBFTOPO_JET3VETO_Tag1, RECO_VBFTOPO_JET3_Tag0, RECO_VBFTOPO_JET3_Tag1, 
                             RECO_WHLEP, RECO_ZHLEP, RECO_VHLEPLOOSE, RECO_VHMET, RECO_VHHAD, RECO_TTH_LEP, RECO_TTH_HAD };

        DiPhotonTagBase();
        virtual ~DiPhotonTagBase(); 
        DiPhotonTagBase( edm::Ptr<DiPhotonCandidate>, DiPhotonMVAResult );
        DiPhotonTagBase( edm::Ptr<DiPhotonCandidate>, edm::Ptr<DiPhotonMVAResult> );
        const edm::Ptr<DiPhotonCandidate> diPhoton() const { return dipho_; }

        const flashgg::Photon *leadingPhoton() const { return dipho_->leadingPhoton(); }
        const flashgg::Photon *subLeadingPhoton() const { return dipho_->subLeadingPhoton(); }
        const flashgg::SinglePhotonView *leadingView() const { return dipho_->leadingView(); }
        const flashgg::SinglePhotonView *subLeadingView() const { return dipho_->subLeadingView(); }

        const DiPhotonMVAResult diPhotonMVA() const { return mva_result_; }
        int diPhotonIndex() const {return diPhotonIndex_;}
        void setDiPhotonIndex( int i ) { diPhotonIndex_ = i; }
        float sumPt() const { return this->diPhoton()->sumPt() ;}
        bool operator <( const DiPhotonTagBase &b ) const;
        operator int() const { return categoryNumber(); }
        virtual DiPhotonTagBase *clone() const { return ( new DiPhotonTagBase( *this ) ); }
        void setCategoryNumber( int value ) { category_number_ = value; }
        int categoryNumber() const { return category_number_; }
        void setTagTruth( const edm::Ptr<TagTruthBase> value ) { truth_ = value; }
        const edm::Ptr<TagTruthBase> tagTruth() const { return truth_; }
        void setSystLabel( const std::string label ) { systLabel_ = label; }
        std::string systLabel() const { return systLabel_; }
        bool hasSyst( const string &label ) const { return ( systLabel_ == label );}
        void setIsGold ( int runNumber );
        void setIsGoldMC( bool isGold ) { isGold_ = isGold; }
        bool isGold() const { return isGold_; }
        virtual DiPhotonTagBase::tag_t tagEnum() const { return DiPhotonTagBase::kUndefined; }
        int stage1recoEnum() const { return stage1recoTag_; }
        unsigned nOtherTags() const { 
            assert(otherTagTypes_.size() == otherTagCategories_.size());
            assert(otherTagTypes_.size() == otherTagIndices_.size());
            return otherTagTypes_.size(); 
        }
        void addOtherTag( const DiPhotonTagBase& other ) { 
            otherTagTypes_.push_back(other.tagEnum());
            otherTagCategories_.push_back(other.categoryNumber());
            otherTagIndices_.push_back(other.diPhotonIndex());
        }
        void addOtherTags( std::vector<std::tuple<DiPhotonTagBase::tag_t,int,int> > others ) { 
            for (unsigned i = 0 ; i < others.size() ; i++) {
                otherTagTypes_.push_back(std::get<0>(others[i]));
                otherTagCategories_.push_back(std::get<1>(others[i]));
                otherTagIndices_.push_back(std::get<2>(others[i]));
            }
        }
        DiPhotonTagBase::tag_t otherTagType( unsigned i ) const { return otherTagTypes_[i]; }
        int otherTagCategory( unsigned i ) const { return otherTagCategories_[i]; }
        int otherTagDiPhotonIndex ( unsigned i ) const { return otherTagIndices_[i]; }
        float ggHweightCentralised( std::string weightName ) const;

        void computeStage1Kinematics( const edm::Handle<edm::View<flashgg::Jet> > & jets, float ptV = -1., float lepphi1 = -999., float lepeta1 = -999., float lepphi2 = -999., float lepeta2 = -999. );
        string stage1KinematicLabel() const { return stage1KinematicLabel_; }
        void setStage1KinematicLabel( const string label ) { stage1KinematicLabel_ = label; }
        void setStage1recoTag( const int tag ) { stage1recoTag_ = tag; }
    private:
        DiPhotonMVAResult mva_result_;
        int category_number_;
        int diPhotonIndex_;
        edm::Ptr<DiPhotonCandidate> dipho_;
        edm::Ptr<TagTruthBase> truth_;
        string systLabel_;
        bool isGold_;
        //        std::vector<std::tuple<DiPhotonTagBase::tag_t,int,int> > otherTags_; // (type,category,diphoton index) 
        std::vector<DiPhotonTagBase::tag_t> otherTagTypes_;
        std::vector<int> otherTagCategories_;
        std::vector<int> otherTagIndices_;
        // weights for 2017 ggH WG1 uncertainty scheme
        // with central object weight applied, unlike TagTruthBase version
        // order: THU_ggH_Mu, THU_ggH_Res, THU_ggH_Mig01, THU_ggH_Mig12, THU_ggH_VBF2j, THU_ggH_VBF3j, THU_ggH_PT60, THU_ggH_PT120, THU_ggH_qmtop
        std::map<std::string,float> ggHweightsCentralised_;

        string stage1KinematicLabel_;
        int stage1recoTag_;

    };

}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

