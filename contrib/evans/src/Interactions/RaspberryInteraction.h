//
// Created by josh evans on 3 December 2024.
// Name is temporary (aka, we're going to keep saying we should get around to renaming it, then not do that)
// code will be pilfered (sorry, "adapted") from both romano/src/Interactions/PatchyShapeInteraction (on oxDNA_torsion repo)
// and a little from Lorenzo Rovigatti's interactions
//

#ifndef RASPBERRYINTERACTION_H
#define RASPBERRYINTERACTION_H

#include "Interactions/BaseInteraction.h"
#include <unordered_map>

// constants from Flavio Romano's PatchyShapeInteraction
// i'm not 100% sure why r and rc aren't 1
#define PLEXCL_S   1.0f
#define PLEXCL_R   0.9053f
#define PLEXCL_B   677.505671539f
// cutoff for repulsive interaction. if r ^ 2 < ((r1+r2) * PLEXCL_RC) ^ 2, will calculate
#define PLEXCL_RC  0.99888f
#define PLEXCL_EPS 2.0f

// todo: make this dynamic or smth!!!! or at least warn when it will cause issues
#define PATCHY_CUTOFF 0.18f


// hash unordered pairs
struct UnorderedPairHash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        // Ensure that (a, b) and (b, a) have the same hash value
        int a = std::min(p.first, p.second);
        int b = std::max(p.first, p.second);
        // Use a simple hash combination technique
        return std::hash<int>()(a) ^ (std::hash<int>()(b) << 1);
    }
};


// Custom equality function for unordered pairs
struct UnorderedPairEqual {
    bool operator()(const std::pair<int, int>& p1, const std::pair<int, int>& p2) const {
        // Since the pair is unordered, (a, b) == (b, a)
        return (p1.first == p2.first && p1.second == p2.second) ||
               (p1.first == p2.second && p1.second == p2.first);
    }
};

/**
* @brief interaction between patchy particles with multiple repulsive points
* we aren't going to use classes for particles and patches to better do forward compatibility with CUDA
*/
class RaspberryInteraction : public BaseInteraction {
protected:
    /**
     * Patch type members:
     * - patch unique ID
     * - patch position
     * - patch alignment (a1 vector). expected to be normalized
     * - patch color (used w/ interaction matrix to find interaction strength)
     * - state variable (for future use)
     * - activation variable (for future use)
     * - patch poly-T spacer length (in nucleotides) (for future use)
     * - patch sticky sequence (for future use)
     */
    using Patch = std::tuple<
            int, //
            LR_vector, // position
            LR_vector, // orientation
            int, // color
            float, //strength
            int, //state variable
            int, //activation variable
            number, // polyT (aka sigma)
            std::string // sticky sequence
            >;
#define PPATCH_TYPE_ID 0
#define PPATCH_POS 1
#define PPATCH_ORI 2
#define PPATCH_COLOR 3
#define PPATCH_STATE 5
#define PPATCH_INT_DIST 7
    /**
     * get in loser we're doing native multidentate patches
     * this is mostly an organizational thing
     * tbd
     */
    using PatchGroupType = std::tuple<std::vector< int, Patch&>>;

    /**
     * Repulsion point (not implemented yet)
     * - repulsion point position, repulsion distance
     */
    using RepulsionPoint = std::tuple<int, LR_vector, number>;
#define REPULSION_ID 0
#define REPULSION_COORDS 1
#define REPULSION_DIST 2

    /**
     * - particle type
     * - number of instances
     * - list of patch ids
     * - list of repulsion point ids
     * - state size
     */
    using ParticleType = std::tuple<int,
                                    int,
                                    std::vector<int>,
                                    std::vector<int>,
                                    int>;
#define PTYPE_ID  0
#define PTYPE_INST_COUNT 1
#define PTYPE_PATCH_IDS  2
#define PTYPE_REP_PTS 3
#define PTYPE_STATE_SIZE 4

    // repulsion points
    std::vector<RepulsionPoint> m_RepulsionPoints;
    // patch types
    std::vector<Patch> m_PatchesTypes;
    // particle types
    std::vector<ParticleType> m_ParticleTypes;
    // should be length _N, each value is a particle type ID in the above m_ParticleTypes
    std::vector<int> m_ParticleList;

    // runtime variables
    std::vector<std::vector<int>> m_ParticleStates; // todo

    // bonds. this is tricky
    // use this name here for clarity
    using ParticlePatch = std::pair<int,int>;

    // list of lists
    // each item in the outer list of particles, the inner list is of patches on the particles
    std::vector<std::vector<ParticlePatch>> m_PatchyBonds;

    // patch color interactions

//    std::unordered_map<std::pair<int, int>, number,  UnorderedPairHash, UnorderedPairEqual> m_PatchColorInteractions;

// todo impl option to specify potential version in input file
    // lorenzo version
    // a PatchPatch object describes how two patches interact
    // there is - by design - very redundant, to minimize required calculations at runtime
    // order: sigma_ss, rcut_ss a_part, b_part, eps
//    using PatchPatch = std::tuple<number, number, number, number, number>;
//#define PP_INT_RCUT_SS  0
//#define PP_INT_SIGMA_SS 1
//#define PP_INT_A_PART   2
//#define PP_INT_B_PART   3
//#define PP_INT_EPS      4

    using PatchPatch = std::tuple<number, number, number, number >;
#define PP_INT_EPS          0
#define PP_INT_ALPHA_POW    1
#define PP_MAX_DIST_SQR     2
#define PP_E_CUT            3

    // i've gone back and forth a LOT about how to work these, settled on this method, for now
    // i don't think this hash or equal function are very fast
    // in this case for speed i am using patch type ids as my hash, color should not be discussed
    // outside initialization
    std::unordered_map<std::pair<int,int>, PatchPatch, UnorderedPairHash, UnorderedPairEqual> m_PatchPatchInteractions;

    number m_nPatchyBondEnergyCutoff;
    number m_nDefaultAlpha;
    bool _has_read_top; // flag to avoid double-reading the top file
public:
    RaspberryInteraction();
    virtual ~RaspberryInteraction();

    virtual void get_settings(input_file &inp);
    // initializes constants for the interaction
    virtual void init();
    // allocate particles
    virtual void allocate_particles(std::vector<BaseParticle *> &particles);

    // utility functions
    int numParticles() const;
    LR_vector getParticlePatchPosition(BaseParticle* p, int patch_idx) const;
    LR_vector getParticlePatchAlign(BaseParticle* p, int patch_idx) const;
    LR_vector getParticleInteractionSitePosition(BaseParticle* p, int int_site_idx) const;
    const RaspberryInteraction::Patch& getParticlePatchType(BaseParticle* p, int patch_idx) const;
//    number get_r_max_sqr(const int &intSite1, const int &intSite2) const;
//    number get_r_sum(const int &intSite1, const int &intSite2) const;

    // pair interaction functions
    virtual number pair_interaction(BaseParticle *p,
                                    BaseParticle *q,
                                    bool compute_r = true,
                                    bool update_forces = false);
    virtual number pair_interaction_bonded(BaseParticle *p,
                                           BaseParticle *q,
                                           bool compute_r = true,
                                           bool update_forces = false);
    virtual number pair_interaction_nonbonded(BaseParticle *p,
                                              BaseParticle *q,
                                              bool compute_r = true,
                                              bool update_forces = false);

    virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
    virtual int get_N_from_topology();
    virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

    int getPatchBondEnergyCutoff() const {
        return m_nPatchyBondEnergyCutoff;
    }

//    virtual void begin_energy_computation() ;

    // interaction functions
    number repulsive_pt_interaction(BaseParticle *p, BaseParticle *q, bool update_forces);

    number patchy_pt_interaction(BaseParticle *p, BaseParticle *q, bool update_forces);

    void readPatchString(const std::string &patch_line);

    // better to pass refs to objects for speed reasons
    bool patches_can_interact(BaseParticle *p, BaseParticle *q,
                              int ppatch_idx, int qpatch_idx) const;
    bool patch_is_active(BaseParticle* p, const Patch& patch_type) const;
    bool patch_types_interact(const Patch &ppatch_type, const Patch &qpatch_type) const;
//    number patch_types_eps(const Patch &ppatch_type, const Patch &qpatch_type) const;


    // methods for handling locking
    bool is_bound_to(int p, int ppatch_idx, int q, int qpatch_idx) const;
    const ParticlePatch& patch_bound_to(BaseParticle* p, int patch_idx) const;
    bool patch_bound(BaseParticle* p, int patch_idx) const;
    void set_bound_to(int p, int ppatch_idx, int q, int qpatch_idx);
    void clear_bound_to(int p, int ppatch_idx);

    static number compute_energy(number patch_dist_sqr, number alpha_exp, number &r8b10) ;
};

std::string readLineNoComment(std::istream& inp);
LR_vector parseVector(const std::string& sz);



#endif //RASPBERRYINTERACTION_H
