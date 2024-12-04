//
// Created by josh evans on 3 December 2024
// Name is temporary (aka, we're going to keep saying we should get around to renaming it, then not do that)
// code will be pilfered (sorry, "adapted") from both romano/src/Interactions/PatchyShapeInteraction (on oxDNA_torsion repo)
// and a little from Lorenzo Rovigatti's interactions
//

#include <sstream>
#include "RaspberryInteraction.h"

RaspberryInteraction::RaspberryInteraction()  : BaseInteraction(){
    // we actually don't want to add interactions to map here, since for raspberry particles
    // this will depend somewhat on inputs
}

RaspberryInteraction::~RaspberryInteraction() {

}

/**
 * function to allocate particles
 */
void RaspberryInteraction::init() {

}

void RaspberryInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
    int i_type = 0;
    int type_count = 0;
    for (int i = 0; i < particles.size(); i++){
        m_ParticleList[i] = i_type;
        type_count++;
        if (type_count == std::get<0>(m_ParticleTypes[i_type])){
            i_type++;
            type_count = 0;
        }
    }
}

int RaspberryInteraction::numParticles() const {
    return m_ParticleList.size();
}

// interaction sites in particle positions

/***
 * retrieves the position of patch with the specified index on the specified particle
 * @param p
 * @param patch_idx INDEX of the patch WITHIN THE PARTICLE TYPE
 * @return
 */
LR_vector RaspberryInteraction::getParticlePatchPosition(BaseParticle *p, int patch_idx) {
    // patch positions are listed first
    return p->int_centers[patch_idx];
}

/**
 * retrieves the alignment vector of patch with the specified position on the specified particle
 * @param p
 * @param patch_idx INDEX of the patch WITHIN THE PARTICLE TYPE
 * @return
 */
LR_vector RaspberryInteraction::getParticlePatchAlign(BaseParticle *p, int patch_idx) {
    ParticleType* particleType = &m_ParticleTypes[m_ParticleList[p->get_index()]];
    // patch orientations are listed after patch positions
    return p->int_centers[patch_idx + std::get<1>(*particleType).size() ];
}

/**
 *
 * @param p
 * @param int_site_idx INDEX of the interaction sitw WITHIN THE PARTICLE TYPE
 * @return
 */
LR_vector RaspberryInteraction::getParticleInteractionSitePosition(BaseParticle *p, int int_site_idx) {
    ParticleType* particleType = &m_ParticleTypes[m_ParticleList[p->get_index()]];
    // interaction sites are listed after the patch geometry
    return p->int_centers[int_site_idx + 2 * std::get<1>(*particleType).size() ];
}

number
RaspberryInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
    return 0;
}

number RaspberryInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r,
                                                     bool update_forces) {
    // going off Lorenzo's example we are treating all particles
    return 0.;
}

number RaspberryInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r,
                                                        bool update_forces) {
    return 0;
}

/**
 * reads the topology file
 * @param N_strands
 * @param particles
 */
void RaspberryInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
    // open topology file
    std::ifstream topology(_topology_filename, std::ios::in);

    // resize list of particle type IDs to match particle pointers
    m_ParticleList.resize(particles.size());
    // set number of "strands"
    *N_strands = m_ParticleList.size();
    // read header
    std::string first_line = readLineNoComment(topology);
    // skip first line entirely (for now)

    // let's read all lines first
    std::vector<std::string> particle_type_lines;
    std::vector<std::string> patch_type_lines;
    std::vector<std::string> patch_group_lines;
    std::vector<std::string> repulsion_pt_lines;
    std::vector<std::string> signal_passing_operations;

    std::string sz;
    while (std::getline(topology, sz)){
        // read line
        std::getline(topology, sz);
        sz = Utils::trim(sz);
        // if line is not blank or a comment
        if (sz.length() > 0 && sz[0] != '#') {
            // if line is patch descriptor
            if (sz.substr(0, 2) == "iP") {
                patch_type_lines.push_back(sz);
            }
                // particle type ("corpuscule") descriptor
            else if (sz.substr(0, 2) == "iC") {
                patch_type_lines.push_back(sz);
            }
                // operation descriptor
            else if (sz.substr(0, 2) == "iO") {
                signal_passing_operations.push_back(sz);
            } else if (sz.substr(0, 2) == "iG") {
                patch_group_lines.push_back(sz);
            }
            else if (sz.substr(0,2) == "iR"){
                repulsion_pt_lines.push_back(sz);
            } else {
                throw oxDNAException("Malformed topology! Line `" + sz + "` does not clearly convey information");
            }
        }
    }
    // resize patch type vector
    m_PatchesTypes.resize(patch_type_lines.size());
    // read patches
    for (int i = 0; i < patch_type_lines.size(); i++){
        std::stringstream ss(patch_type_lines[i].substr(2, patch_type_lines[i].size()-2));
        int iTypeID, iColor, iState, iActivation;
        std::string vecStr, oriStr;
        LR_vector position, orientation;

        // read patch type ID
        if (!(ss >> iTypeID)){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }

        // read patch position
        if (!(ss >> vecStr)){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }
        try {
            position = parseVector(vecStr);
        } catch (oxDNAException &e){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }

        // read patch orientation
        if (!(ss >> oriStr)){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }
        try {
            orientation = parseVector(oriStr);
        } catch (oxDNAException &e){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }

        // read color
        if (!(ss >> iColor)){
            throw oxDNAException("Invalid patch type str `" + patch_type_lines[i] + "`!");
        }

        // following values are optional!
        if (!(ss >> iState && ss >> iActivation)){
            iState = 0;
            iActivation =0;
        }
        if (iTypeID >= m_PatchesTypes.size()){
            throw oxDNAException("Invalid patch type ID %d", iTypeID);
        }
        // TODO: polyTs, sequence?
        m_PatchesTypes[iTypeID] = {iTypeID, position, orientation, iColor, iState, iActivation, 0, ""};
    }

    // read patch groups
    for (int i = 0; i < patch_group_lines.size(); i++){
        std::stringstream ss(patch_group_lines[i].substr(2, patch_group_lines[i].size()-2));
        std::vector<int> patches;
        int p;
        while (ss >> p){
            patches.push_back(p);
        }
        // todo: if patch groups do anything put that here
    }

    // read repulsion points
    m_RepulsionPoints.resize(repulsion_pt_lines.size());
    for (int i = 0; i < repulsion_pt_lines.size(); i++){
        std::stringstream ss(repulsion_pt_lines[i].substr(2, repulsion_pt_lines[i].size()-2));
        LR_vector position;
        number r;
        std::string pos_str;
        if (!(ss >> pos_str >> r)){
            throw oxDNAException("Invalid repulsion point type str `" + repulsion_pt_lines[i] + "`!");
        }
        try {
            position = parseVector(pos_str);
        } catch (oxDNAException &e){
            throw oxDNAException("Invalid repulsion point type str `" + repulsion_pt_lines[i] + "`!");
        }
        m_RepulsionPoints[i] = {i, position, r};
    }

    // read particle types
    m_ParticleTypes.resize(particle_type_lines.size());
    for (int i = 0; i < particle_type_lines.size(); i++){
        std::stringstream ss(particle_type_lines[i].substr(2, particle_type_lines[i].size()-2));
        int iParticleType;
        std::string patch_id_strs, interaction_pt_id_strs;
        if (!(ss >> iParticleType >> std::get<0>(m_ParticleTypes[iParticleType]) >> patch_id_strs >> interaction_pt_id_strs)){
            throw oxDNAException("Invalid particle type str `" + particle_type_lines[i] + "`!");
        }
        if (iParticleType >= m_ParticleTypes.size()){
            throw oxDNAException("Invalid particle type ID %d", iParticleType);
        }
        std::vector<std::string> patch_id_strs_list = Utils::split(patch_id_strs,',');
        std::vector<std::string> int_pt_strs_list = Utils::split(interaction_pt_id_strs,',');
        // process patch IDs
        for (std::string& sz : patch_id_strs_list){
            std::get<1>(m_ParticleTypes[iParticleType]).push_back(std::stoi(sz));
        }
        // process interaction points
        for (std::string &sz : int_pt_strs_list){
            std::get<2>(m_ParticleTypes[iParticleType]).push_back(std::stoi(sz));
        }
        // TODO: operations
    }

    // todo: read operations

    // close topology file
    topology.close();
    // TODO: CUDA version of this func will need to cache vs. max patches
}

void RaspberryInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {
    // ????
}

/**
 * a lot of this is pilfered from DetailedPatchySwapInteraction::_patchy_two_body_point
 * and DetailedPatchySwapInteraction::_spherical_patchy_two_body in contrib/rovigatti
 * @param p
 * @param q
 * @param update_forces
 * @return
 */
number RaspberryInteraction::repulsive_pt_interaction(BaseParticle *p, BaseParticle *q, bool update_forces) {
    int p_type = m_ParticleList[p->get_index()];
    int q_type = m_ParticleList[q->get_index()];

    number energy = 0.;

    // TODO: make this combinatorics problem less shit
    for (int pp = 0; pp < std::get<PTYPE_REP_PTS>(m_ParticleTypes[p_type]).size(); pp++){
        for (int qq = 0; qq < std::get<PTYPE_REP_PTS>(m_ParticleTypes[q_type]).size(); qq++) {
            // find interaction site positions for p & q
            LR_vector ppos = getParticleInteractionSitePosition(p, pp);
            LR_vector qpos = getParticleInteractionSitePosition(q, qq);
            // repulsive spheres have variable sizes,
            // we precompute the rmaxes-squareds for each pair of interaction types
            number rmax_dist_sqr = get_r_max_sqr(
                    std::get<PTYPE_REP_PTS>(m_ParticleTypes[p_type])[pp],
                    std::get<PTYPE_REP_PTS>(m_ParticleTypes[q_type])[qq]
            );

            LR_vector patch_dist = _computed_r + qpos - ppos;
            number r_patch_sqr = patch_dist.norm();
            // if the distance between the two interaction points is greater than the cutoff
            if (r_patch_sqr < rmax_dist_sqr) {
                // unlike in other impls our model does not assume repulsive spheres have radius=0.5

                number lj_rfac = get_r_lj_offset(
                        std::get<PTYPE_REP_PTS>(m_ParticleTypes[p_type])[pp],
                        std::get<PTYPE_REP_PTS>(m_ParticleTypes[q_type])[qq]
                );

                // this is a bit faster than calling r.norm()
                number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
                number energy = (number) 0;


                sigma *= (2*lj_rfac);
                rstar *= (2*lj_rfac);
                b /= SQR(2*lj_rfac);
                rc *= (2*lj_rfac);

                if(rnorm < SQR(rc)) {
                    if(rnorm > SQR(rstar)) {
                        number rmod = sqrt(rnorm);
                        number rrc = rmod - rc;
                        energy = PLEXCL_EPS * b * SQR(rrc);
                        if(update_forces) force = -r * (2 * PLEXCL_EPS * b * rrc / rmod);
                    }
                    else {
                        number tmp = SQR(sigma) / rnorm;
                        number lj_part = tmp * tmp * tmp;
                        energy = 4 * PLEXCL_EPS * (SQR(lj_part) - lj_part);
                        if(update_forces) force = -r * (24 * PLEXCL_EPS * (lj_part - 2*SQR(lj_part)) / rnorm);
                    }
                }

                if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

                return energy;

                // check that the center-center distance isn't less than the sum of the radii; this should not be allowed
                assert((r_patch_sqr - int_site_rad_sum_sqr) > 0 );

                // very much stole the force calculations from Lorenzo, because there is no way i get those right

                // inverse square of the radius
                number inv_sqr_r = 1. / (r_patch_sqr - int_site_rad_sum_sqr);
                // compute 1 / r ^ 6 by computing (1/(r ^ 2))^3
                number lj_part = inv_sqr_r * inv_sqr_r * inv_sqr_r;
                // energy from lennard-jones potential
                number e_lj = 4 * (SQR(lj_part) - lj_part) + 1.0;

                // update energy
                energy += e_lj;

                number ir2 = 1. / sqr_r;
                number lj_part = ir2 * ir2 * ir2;
                energy = 4 * (SQR(lj_part) - lj_part) + 1.0 - _spherical_attraction_strength - _spherical_E_cut;
                if(update_forces) {
                    LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);
                    p->force -= force;
                    q->force += force;

                    _update_stress_tensor(p->pos, -force);
                    _update_stress_tensor(p->pos + _computed_r, force);
                }

                if (update_forces) {
                    // force is the derivitive of energy w/ respect to time
                    //
                    LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);

                    // compute torque from force vectors
                    LR_vector p_torque = p->orientationT * ppos.cross(tmp_force);
                    LR_vector q_torque = q->orientationT * qpos.cross(tmp_force);

                    // update forces for particles
                    p->force -= tmp_force;
                    q->force += tmp_force;

                    // update torque for particles
                    p->torque -= p_torque;
                    q->torque += q_torque;

                    // no idea what this does
                    _update_stress_tensor(p->pos, -tmp_force);
                    _update_stress_tensor(p->pos + _computed_r, tmp_force);
                }
            }
        }
    }
    return energy;
}

number RaspberryInteraction::patchy_kf_interaction(BaseParticle *p, BaseParticle *q,
                                                   bool compute_r, bool update_forces) {
    return 0;
}

number RaspberryInteraction::patchy_lj_interaction(BaseParticle *p, BaseParticle *q,
                                                   bool compute_r, bool update_forces) {
    return 0;
}

number RaspberryInteraction::patchy_pt_like_interaction(BaseParticle *p, BaseParticle *q,
                                                        bool compute_r, bool update_forces) {
    return 0;
}


/**
 * retrieves the "lennard jones scale factor" (named after Dr. John Scalefactor)
 * the default lennard-jones potential assumes a energy y-intercept at r=1
 * to offset that we need to add the sum of the radii
 * @param intSite1
 * @param intSite2
 * @return
 */
number RaspberryInteraction::get_r_lj_offset(int intSite1, int intSite2) const {
    return 0;
}

/**
 * retrieves the matimum distance at which two repulsive interaction point types (that's a mouthful!)
 * will interact
 * @param intSite1 type ID of first interaction point
 * @param intSite2 type ID of second interaction point
 * @return
 */
number RaspberryInteraction::get_r_max_sqr(int intSite1, int intSite2) const {
    return 0;
}


std::string readLineNoComment(std::istream& inp){
    std::string sz;
    do {
        // read first line
        std::getline(inp, sz);
        sz = Utils::trim(sz);

    } while (sz.length() > 0 && sz[0] == '#'); // skip comments
    return sz;
}

LR_vector parseVector(const std::string& vec){
    number tmpf[3];
    int tmpi = sscanf(vec.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
    if(tmpi != 3)
        throw oxDNAException ("Could not parse vector %s. Aborting", vec.c_str());
    return {tmpf[0], tmpf[1], tmpf[2]};
}

extern "C" RaspberryInteraction* make_RaspberryInteraction() {
    return new RaspberryInteraction();
}