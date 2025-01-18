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

void RaspberryInteraction::get_settings(input_file &inp) {
    BaseInteraction::get_settings(inp);
}


void RaspberryInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
    int i_type = 0;
    int type_count = 0;
    for (int i = 0; i < particles.size(); i++){
        // assign particle type
        m_ParticleList[i] = i_type;

        particles[i] = new BaseParticle();
        particles[i]->index = i;
        particles[i]->strand_id = i;
        particles[i]->type = i_type;

        // ordering of interaction centers is important!
        // first we do patches, as pairs of position, a1
        // then repuslive spheres!
        int num_int_centers = std::get<PTYPE_REP_PTS>(m_ParticleTypes[i_type]).size() + 2*std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type]).size();
        particles[i]->int_centers.resize(num_int_centers);
        // use same indexer for repulsion points & patches
        int j = 0;

        for (; j < 2*std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type]).size(); j += 2){
            // two "int centers" for each patch: position and a1
            int idx = j /2 ;
            assert(idx < std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type]).size());
            int iPatch = std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type])[idx];
            assert(iPatch < m_PatchesTypes.size());
            // assign position from patch info
            particles[i]->int_centers[j] = std::get<1>(m_PatchesTypes[iPatch]);
            // assign a1 from patch info
            particles[i]->int_centers[j+1] = std::get<2>(m_PatchesTypes[iPatch]);
        }

        for (; j < num_int_centers; j++){
            // get index of info in m_RepulsionPoints
            int idx = j - 2*std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type]).size();
            int iRepulsionPoint = std::get<PTYPE_REP_PTS>(m_ParticleTypes[i_type])[idx];
            assert(iRepulsionPoint < m_RepulsionPoints.size());
            // set int center to repulsion point position
            particles[i]->int_centers[j] = std::get<1>(m_RepulsionPoints[iRepulsionPoint]);
        }



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
    return p->int_centers[patch_idx + std::get<PTYPE_REP_PTS>(*particleType).size() ];
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
    return p->int_centers[int_site_idx + 2 * std::get<PTYPE_REP_PTS>(*particleType).size() ];
}

number RaspberryInteraction::pair_interaction(BaseParticle *p,
                                              BaseParticle *q,
                                              bool compute_r,
                                              bool update_forces) {
    return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

number RaspberryInteraction::pair_interaction_bonded(BaseParticle *p,
                                                     BaseParticle *q,
                                                     bool compute_r,
                                                     bool update_forces) {
    // going off Lorenzo's example we are treating all particles as nonbonded
    return 0.;
}

number RaspberryInteraction::pair_interaction_nonbonded(BaseParticle *p,
                                                        BaseParticle *q,
                                                        bool compute_r,
                                                        bool update_forces) {
    number e = 0.;
    if(compute_r) {
        _computed_r = _box->min_image(p->pos, q->pos);
    }

    e += repulsive_pt_interaction(p, q, update_forces);

    return e;
}

int RaspberryInteraction::get_N_from_topology() {
    std::ifstream topology(_topology_filename, std::ios::in);
    if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
    std::string sz, _;
    int N = 0;
    int n;

    // topology file should only be a few dozen lines, we can just skim it
    while (std::getline(topology, sz)) {
        if (sz.length() > 2 && sz.substr(0,2) == "iC"){
            // yuck way to do this but i don't care
            std::stringstream ss(sz);
            ss >> _ >> _ >> n;
            N += n;
        }
    }
    return N;
}

/**
 * reads the topology file
 * @param N_strands
 * @param particles
 */
void RaspberryInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
    // open topology file
    std::ifstream topology(_topology_filename, std::ios::in);
    if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);

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
        sz = Utils::trim(sz);
        // if line is not blank or a comment
        if (sz.length() > 0 && sz[0] != '#') {
            // if line is patch descriptor
            if (sz.substr(0, 2) == "iP") {
                patch_type_lines.push_back(sz);
            }
            // particle type ("corpuscule") descriptor
            else if (sz.substr(0, 2) == "iC") {
                particle_type_lines.push_back(sz);
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
        readPatchString(patch_type_lines[i]);
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
//    // cache sums of all possible repulsion point pairs
//    // todo: warn if this is too big?
//    for (int i = 0; i < m_RepulsionPoints.size(); i++){
//        for (int j = i; j < m_RepulsionPoints.size(); j++){
//            m_RSums[{i, j}] = std::get<2>(m_RepulsionPoints[i]) + std::get<2>(m_RepulsionPoints[j]);
//        }
//    }

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
            int patch = std::stoi(sz);
            assert(patch < m_PatchesTypes.size());
            std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[iParticleType]).push_back(patch);
        }
        // process interaction points
        for (std::string &sz : int_pt_strs_list){
            int repulsionPt = std::stoi(sz);
            assert(repulsionPt < m_RepulsionPoints.size());
            std::get<PTYPE_REP_PTS>(m_ParticleTypes[iParticleType]).push_back(repulsionPt);
        }
        // TODO: operations
    }

    // todo: read operations

    // close topology file
    topology.close();
    // TODO: CUDA version of this func will need to cache vs. max patches


    allocate_particles(particles);
}

void RaspberryInteraction::readPatchString(const std::string& patch_line) {
    std::string patch_infostr = patch_line.substr(3);
    std::stringstream ss(patch_infostr);
    int iTypeID, iColor, iState, iActivation;
    float fStrength;
    std::string vecStr, oriStr;
    LR_vector position, orientation;

    // order of line should be UID color strength position a1 statevar activationvar

    // read patch type ID
    if (!(ss >> iTypeID)){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }

    if (!(ss >> fStrength)){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }
    if (fStrength != 1.0){
        OX_LOG(Logger::LOG_WARNING, "Strength is required for compatability with Patchy Helix Bundle but i don't currently use it");
    }

    // read color
    if (!(ss >> iColor)){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }

    // read patch position
    if (!(ss >> vecStr)){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }
    try {
        position = parseVector(vecStr);
    } catch (oxDNAException &e){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }

    // read patch orientation
    if (!(ss >> oriStr)){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
    }
    try {
        orientation = parseVector(oriStr);
    } catch (oxDNAException &e){
        throw oxDNAException("Invalid patch type str `" + patch_line + "`!");
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
    m_PatchesTypes[iTypeID] = {iTypeID, position, orientation, iColor, fStrength, iState, iActivation, 0, ""};
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

    number e_lj;
    number energy = 0.;
    number sigma, rstar, b, rc, rsum;
    LR_vector force;
    LR_vector p_torque, q_torque;
    int ppidx, qqidx;
    number r_patch;
    // TODO: make this combinatorics problem less shit
    for (int pp = 0; pp < std::get<PTYPE_REP_PTS>(m_ParticleTypes[p_type]).size(); pp++){
        for (int qq = 0; qq < std::get<PTYPE_REP_PTS>(m_ParticleTypes[q_type]).size(); qq++) {
            // find interaction site positions for p & q
            LR_vector ppos = getParticleInteractionSitePosition(p, pp);
            LR_vector qpos = getParticleInteractionSitePosition(q, qq);
            // repulsive spheres have variable sizes,
            // we precompute the rmaxes-squareds for each pair of interaction types
            ppidx = std::get<PTYPE_REP_PTS>(m_ParticleTypes[p_type])[pp];
            qqidx = std::get<PTYPE_REP_PTS>(m_ParticleTypes[q_type])[qq];
            // lookup r-max squared
            // todo: probably possible to precompute these to save a little time
            number rmax_dist_sqr = SQR(std::get<2>(m_RepulsionPoints[ppidx]) + std::get<2>(m_RepulsionPoints[ppidx])) * 1.2;
//            number rmax_dist_sqr = get_r_max_sqr(ppidx, qqidx);

            LR_vector rep_pts_dist = _computed_r + qpos - ppos;
            number rep_pts_dist_sqr = rep_pts_dist.norm();
            // if the distance between the two interaction points is greater than the cutoff
                if (rep_pts_dist_sqr < rmax_dist_sqr) {
                // unlike in other impls our model does not assume repulsive spheres have radius=0.5
                // r-factor = the sum of the radii of the repulsive interaction spheres, minus 1
                rsum = std::get<2>(m_RepulsionPoints[ppidx]) + std::get<2>(m_RepulsionPoints[ppidx]);

                // we use a Lennard Jones function with quadatic smoothin
                // sigma
                sigma = PLEXCL_S;

                // radial cutoff
                rc = PLEXCL_RC * rep_pts_dist_sqr;

                // if the normalized radius is less than the interaction cutoff
                if (rep_pts_dist_sqr < SQR(rc)) {
                    // is this even correct??
                    b = PLEXCL_B / rsum;
                    // compute quadratic smoothing cutoff rstar
                    rstar = PLEXCL_R * rsum;
                    // if r is less than rstar, use the lennard-jones equation
                    if(rep_pts_dist_sqr < SQR(rstar)) {
                        // partial lennard-jones value
                        number lj_partial = (sigma / (sqrt(rep_pts_dist_sqr) - rsum - 1));
                        // lj partial = (sigma/r)^6
                        lj_partial = lj_partial * lj_partial * lj_partial;

                        // compute lennard-jones interaction energy
                        e_lj = 8 * (lj_partial * lj_partial - lj_partial);
                        // update energy
                        energy += e_lj;
                        if(update_forces) {
                            // compute forces
                            // i really hope this is correct
                            force = -rep_pts_dist * (24 * PLEXCL_EPS * (lj_partial - 2 * SQR(lj_partial)) / rep_pts_dist_sqr);
                        }
                    }
                    else {
                        // if r > rstar, use quadratic smoothing
                        // Vsmooth = b(xc - x)^2
                        // actual distance value, computed by taking the square root of rsquared
                        r_patch = sqrt(rep_pts_dist_sqr);
                        number rrc = r_patch - rc;
                        energy += PLEXCL_EPS * b * SQR(rrc);
                        if(update_forces) force = -rep_pts_dist * (2 * PLEXCL_EPS * b * rrc / r_patch);
                    }

                    if (update_forces){
                        // update particle force vectors
                        p->force -= force;
                        q->force += force;
                        // compute torques
                        // i grabbed this eqn from Lorenzo's code
                        p_torque = -p->orientationT * ppos.cross(force);
                        q_torque =  q->orientationT * qpos.cross(force);
                        // update particle torque vectors
                        p->torque += p_torque;
                        q->torque += q_torque;
                    }
                }
                // if the radius is greater than the cutoff, we can just ignore it
            }
        }
    }
    return energy;
}

number RaspberryInteraction::patchy_kf_interaction(BaseParticle *p, BaseParticle *q,
                                                   bool compute_r, bool update_forces) {
    throw ;
}

number RaspberryInteraction::patchy_lj_interaction(BaseParticle *p, BaseParticle *q,
                                                   bool compute_r, bool update_forces) {
    return 0;
}

number RaspberryInteraction::patchy_pt_like_interaction(BaseParticle *p, BaseParticle *q,
                                                        bool compute_r, bool update_forces) {
    return 0;
}


///**
// * returns the sum of the two radii of the interaction sites given
// * @param intSite1 type ID of first interaction point
// * @param intSite1 type ID of first interaction point
// * @return
// */
//number RaspberryInteraction::get_r_sum(const int &intSite1, const int &intSite2) const {
//    return m_RSums.at({intSite1, intSite2});
//}
//
///**
// * retrieves the maximum distance at which two repulsive interaction point types (that's a mouthful!)
// * will interact
// * @param intSite1 type ID of first interaction point
// * @param intSite2 type ID of second interaction point
// * @return
// */
//number RaspberryInteraction::get_r_max_sqr(const int &intSite1, const int &intSite2) const {
//    return 1.2 * SQR(get_r_sum(intSite1, intSite2));
//}



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

extern "C" RaspberryInteraction* make_RaspberryPatchyInteraction() {
    return new RaspberryInteraction();
}