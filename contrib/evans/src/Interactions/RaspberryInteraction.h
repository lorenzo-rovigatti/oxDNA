//
// Created by josh evans on 3 December 2024.
// Name is temporary (aka, we're going to keep saying we should get around to renaming it, then not do that)
// code will be pilfered (sorry, "adapted") from both romano/src/Interactions/PatchyShapeInteraction (on oxDNA_torsion repo)
// and a little from Lorenzo Rovigatti's interactions
//

#ifndef ESTIMATE_TM_PY_RASPBERRYINTERACTION_H
#define ESTIMATE_TM_PY_RASPBERRYINTERACTION_H

#include "Interactions/BaseInteraction.h"
#include <unordered_map>

// constants from Flavio Romano's PatchyShapeInteraction
// i'm not 100% sure why r and rc aren't 1
#define PLEXCL_S   1.0f
#define PLEXCL_R   0.9053f
#define PLEXCL_B   677.505671539f
#define PLEXCL_RC  0.99888f
#define PLEXCL_EPS 2.0f

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
            int,
            LR_vector, // position
            LR_vector, // orientation
            int, // color
            int, //state variable
            int, //activation variable
            unsigned int, // polyT
            std::string // sticky sequence
            >;

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
     * particle type
     * - number of instances
     * - list of repulsion point ids
     * - list of patch ids
     */
    using ParticleType = std::tuple<int, std::vector<int>, std::vector<int> >;
#define PTYPE_ID  0
#define PTYPE_REP_PTS 1
#define PTYPE_PATCH_IDS  2


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
    std::unordered_map<std::pair<int, int>, number, UnorderedPairHash, UnorderedPairEqual> m_RSums;

public:
    enum {
        PATCHY, // patch-patch interaction
        SPHERICAL // spherical repulsive interaction
    };
    RaspberryInteraction();
    virtual ~RaspberryInteraction();

    virtual void get_settings(input_file &inp);
    // initializes constants for the interaction
    virtual void init();
    // allocate particles
    virtual void allocate_particles(std::vector<BaseParticle *> &particles);

    // utility functions
    int numParticles() const;
    LR_vector getParticlePatchPosition(BaseParticle* p, int patch_idx);
    LR_vector getParticlePatchAlign(BaseParticle* p, int patch_idx);
    LR_vector getParticleInteractionSitePosition(BaseParticle* p, int int_site_idx);
    number get_r_max_sqr(const int &intSite1, const int &intSite2) const;
    number get_r_sum(const int &intSite1, const int &intSite2) const;
    number get_int_sites_b(const int &intSite1, const int &intSite2) const;

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
    virtual void check_input_sanity(std::vector<BaseParticle *> &particles);

//    virtual void begin_energy_computation() ;

    // interaction functions
    number repulsive_pt_interaction(BaseParticle *p, BaseParticle *q, bool update_forces);
    number patchy_kf_interaction(BaseParticle* p, BaseParticle* q, bool compute_r = true, bool update_forces = false);
    number patchy_lj_interaction(BaseParticle* p, BaseParticle* q, bool compute_r = true, bool update_forces = false);
    number patchy_pt_like_interaction(BaseParticle* p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
};

std::string readLineNoComment(std::istream& inp);
LR_vector parseVector(const std::string& sz);

#endif //ESTIMATE_TM_PY_RASPBERRYINTERACTION_H
