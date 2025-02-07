//
// Created by josh evans on 3 December 2024
// Name is temporary (aka, we're going to keep saying we should get around to renaming it, then not do that)
// code will be pilfered (sorry, "adapted") from both romano/src/Interactions/PatchyShapeInteraction (on oxDNA_torsion repo)
// and a little from Lorenzo Rovigatti's interactions
//

#include <sstream>
#include "RaspberryInteraction.h"

#define PATCHY_CUTOFF 0.18f

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
    m_nPatchyBondEnergyCutoff = -0.1;
    getInputNumber(&inp, "PATCHY_bond_energy", &m_nPatchyBondEnergyCutoff, 0);
}


void RaspberryInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
    m_ParticleStates.resize(particles.size());
    m_PatchyBonds.resize(particles.size());
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

        // init interaction runtime variables
        // particle states
        // todo: at some point we will need to be able to load these from somewhere
        m_ParticleStates[i].resize(std::get<PTYPE_STATE_SIZE>(m_ParticleTypes[i_type]));
        m_ParticleStates[i][0] = true;
        for (int ii = 1; ii < m_ParticleStates[i].size(); ii++){
            m_ParticleStates[i][ii] = false;
        }

        m_PatchyBonds[i].resize(std::get<PTYPE_PATCH_IDS>(m_ParticleTypes[i_type]).size());
        // for each patch
        for (int ii = 0; ii < m_PatchyBonds.size(); ii++){
            // ParticlePatch = {-1, -1} means no bond
            m_PatchyBonds[i][ii] = {-1, -1};
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
    e += patchy_pt_interaction(p, q, false);

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

    // todo: load manually defined color interactions
    // define color interactions
    for (int pi = 0; pi < m_PatchesTypes.size(); pi++){
        assert(std::get<0>(m_ParticleTypes[pi]) == pi);
        int pi_color = std::get<PPATCH_COLOR>(m_PatchesTypes[pi]);
        number pi_dist = std::get<PPATCH_INT_DIST>(m_PatchesTypes[pi]);
        for (int pj = 0; pj < m_PatchesTypes.size(); pj++){
            assert(std::get<0>(m_ParticleTypes[pj]) == pj);
            int pj_color = std::get<PPATCH_COLOR>(m_PatchesTypes[pj]);
            // use petr + flavio version of default color interactions
            number pj_dist = std::get<PPATCH_INT_DIST>(m_PatchesTypes[pj]);
            if (((abs(pi_color) > 20) && (pi_color + pj_color == 0)) || ((abs(pi_color) <= 20) && (pi_color == pj_color))){
                // todo: patch strengths? ugh
                number sigma_ss = pi_dist + pj_dist;
                number rcut_ss = 1.5 * sigma_ss;
                number B_ss = 1. / (1. + 4. * SQR(1. - rcut_ss / sigma_ss));
                number a_part = -1. / (B_ss - 1.) / exp(1. / (1. - rcut_ss / sigma_ss));
                number b_part = B_ss * pow(sigma_ss, 4.);
                m_PatchPatchInteractions[{pi, pj}] = {
                    sigma_ss, // sigma_ss
                    1.5 * sigma_ss, // rcut_ss
                    a_part, // a_part
                    b_part, //b_part
                    1. //epsilon
                };
            }
        }
    }
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
//          number rmax_dist_sqr = get_r_max_sqr(ppidx, qqidx);

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

number RaspberryInteraction::patchy_pt_interaction(BaseParticle *p, BaseParticle *q, bool update_forces) {
    /**
     * computes the forces and energies from the interactions of patches on two particles
     * I've mostly copied this code from alloassembly, which is itself a modified version of
     * oxDNA_torsion, which is a modified version of Flavio Romano's PatchyShapeInteraction code
     */
    int c = 0;
    LR_vector tmptorquep(0, 0, 0);
    LR_vector tmptorqueq(0, 0, 0);
    const ParticleType& p_type = m_ParticleTypes[p->type];
    const ParticleType& q_type = m_ParticleTypes[q->type];
    number energy = 0.;

    // iter patches on particle p
    for(int ppatch_idx = 0; ppatch_idx < std::get<PTYPE_PATCH_IDS>(p_type).size(); ppatch_idx++) {
        // important: keep pi (index of patch in particle type) distinct from ppatch_tid (patch type global id)
        // lookup patch type id
        int ppatch_tid = std::get<PTYPE_PATCH_IDS>(p_type)[ppatch_idx];

        LR_vector ppatch = getParticlePatchPosition(p, ppatch_idx);

        // iter patches on particle q
        for(int qpatch_idx = 0; qpatch_idx < std::get<PTYPE_PATCH_IDS>(q_type).size(); qpatch_idx++) {
            // note: do NOT use the position from this, since it's not rotated
            int qpatch_tid = std::get<PTYPE_PATCH_IDS>(p_type)[qpatch_idx];

            number e_2patch = 0.;

            // if patches can interact based on color, and both patches are active and not bound to another patch
            if (patches_can_interact(p, q,
                                     ppatch_idx,
                                     qpatch_idx))
            {
                LR_vector qpatch = getParticlePatchPosition(p, ppatch_idx);
                number eps = std::get<PP_INT_EPS>(m_PatchPatchInteractions[{ppatch_tid, qpatch_tid}]);


                // get displacement vector between patches p and q
                LR_vector patch_dist = _computed_r + qpatch - ppatch;
                // compute distancesquared
                number patch_dist_sqr = patch_dist.norm();
                //LR_vector<number> patch_dist_dir = patch_dist / sqrt(dist);
                //number rdist = sqrt(rnorm);
                //LR_vector<number> r_dist_dir = _computed_r / rdist;

                //printf("Patches %d and %d distance %f  cutoff is: %f,\n",pp->patches[pi].id,qq->patches[pj].id,dist,SQR(PATCHY_CUTOFF));

                if (patch_dist_sqr < PATCHY_CUTOFF) {
                    // wait till now to look up either patch align, since we only expect to get here for one p,q pair tops
                    LR_vector qpatch_a1 = getParticlePatchAlign(q, qpatch_idx);
                    LR_vector ppatch_a1 = getParticlePatchAlign(p, ppatch_idx);



                    number r_p = sqrt(patch_dist_sqr);
                    number sigma_ss = std::get<PP_INT_SIGMA_SS>(m_PatchPatchInteractions[{ppatch_tid, qpatch_tid}]);
                    number rcut_ss = std::get<PP_INT_RCUT_SS>(m_PatchPatchInteractions[{ppatch_tid, qpatch_tid}]);
                    number b_part = std::get<PP_INT_B_PART>(m_PatchPatchInteractions[{ppatch_tid, qpatch_tid}]);
                    number a_part = std::get<PP_INT_A_PART>(m_PatchPatchInteractions[{ppatch_tid, qpatch_tid}]);
                    number exp_part = exp(sigma_ss / (r_p - rcut_ss));
                    number tmp_energy = eps * a_part * exp_part * (b_part / patch_dist_sqr - 1.);

                    e_2patch += tmp_energy;


                    if (update_forces) {
                        // i gotta deal with this
                        number force_mod = eps * a_part * exp_part * (4. * b_part / (SQR(r_p) * r_p)) +
                                           sigma_ss * tmp_energy / SQR(r_p - rcut_ss);
                        LR_vector tmp_force = patch_dist * (force_mod / r_p);

                        LR_vector p_torque = p->orientationT * ppatch.cross(tmp_force);
                        LR_vector q_torque = q->orientationT * qpatch.cross(tmp_force);

                        p->force -= tmp_force;
                        q->force += tmp_force;

                        p->torque -= p_torque;
                        q->torque += q_torque;
                    }

                }
            }
            // if energy is great large enough for bond
            if (e_2patch > getPatchBondEnergyCutoff()) {
                // if bond does not alread exist
                if (!is_bound_to(p->index, ppatch_idx, q->index, qpatch_idx)) {
                    // set bound
                    set_bound_to(p->index, ppatch_idx, q->index, qpatch_idx);
                }
            }
            // if energy is not small enough for a bond, but one exists
            else if (is_bound_to(p->index, ppatch_idx, q->index, qpatch_idx)){
                // if patches are bound, but shouldn't be
                clear_bound_to(p->index, ppatch_idx);
                clear_bound_to(q->index, qpatch_idx);

            }
            energy += e_2patch;
        }
    }
    return energy;
}

/**
 *
 * @param p
 * @param q
 * @param ppatch_idx index of patch to check in p (NOT TYPE ID)
 * @param qpatch_idx  index of patch to check in q (NOT TYPE ID)
 * @return
 */
bool RaspberryInteraction::patches_can_interact(BaseParticle *p,
                                                BaseParticle *q,
                                                int ppatch_idx,
                                                int qpatch_idx) const {
    const Patch ppatch_type = m_PatchesTypes[ppatch_idx];
    const Patch qpatch_type = m_PatchesTypes[qpatch_idx];
    // if patch types can't interact, these two patches can't interact. full stop.
    if (!patch_types_interact(ppatch_type, qpatch_type)){
        return false;
    }
    // if either patch isn't active, return false
    if (!patch_is_active(p, ppatch_type)){
        return false;
    }
    if (!patch_is_active(q, qpatch_type)){
        return false;
    }
    if (patch_bound(p, ppatch_idx) && patch_bound_to(p, ppatch_idx) != std::pair<int,int>(q->index, qpatch_idx))){
        return false;
    }
    // if q patch is already locked (we can assume if it's locked, it's not to p patch)
    if (patch_bound(q, qpatch_idx)){
        return false;
    }
    // todo: find and report asymmetric locks
    return true;
}


bool RaspberryInteraction::patch_is_active(BaseParticle* p, const Patch& patch_type) const {
    int state_var = std::get<PPATCH_STATE>(patch_type);
    if (state_var == 0){
        return true; // var 0 = identity variable
    }
    // if state var < 0, var is virtual
    if (state_var < 0){
        return !m_ParticleStates[p->index][-state_var];
    }
    else {
        return m_ParticleStates[p->index][state_var];
    }
}

bool RaspberryInteraction::patch_types_interact(const RaspberryInteraction::Patch &ppatch_type,
                                                const RaspberryInteraction::Patch &qpatch_type) const {
    return m_PatchPatchInteractions.count(
            {
                std::get<3>(ppatch_type),
                std::get<3>(qpatch_type)
            }) > 0 &&
            std::get<PP_INT_EPS>(m_PatchPatchInteractions.at(
            {
                std::get<3>(ppatch_type),
                std::get<3>(qpatch_type)
            }
        )) != 0.;
}

/**
 * checks if patch is bound to anything
 * @param p
 * @param patch_idx
 * @return
 */
bool RaspberryInteraction::patch_bound(BaseParticle *p, int patch_idx) const {
    return patch_bound_to(p, patch_idx).first != -1;
}

/**
 * @param p
 * @param patch_idx index WITHIN THE PARTICLE TYPE, NOT PATCH TYPE ID
 * @return
 */
const RaspberryInteraction::ParticlePatch& RaspberryInteraction::patch_bound_to(BaseParticle *p,
                                                                                int patch_idx) const {
    return m_PatchyBonds[p->get_index()][patch_idx];
}

bool RaspberryInteraction::is_bound_to(int p, int ppatch_idx, int q, int qpatch_idx) const {
    assert(p < m_PatchyBonds.size());
    assert(ppatch_idx < m_PatchyBonds[p].size());
    assert(p < m_PatchyBonds.size());
    assert(qpatch_idx < m_PatchyBonds[q].size());
    return m_PatchyBonds[p][ppatch_idx] == ParticlePatch(q, qpatch_idx);
}

void RaspberryInteraction::set_bound_to(int p, int ppatch_idx, int q, int qpatch_idx) {
    assert(p < m_PatchyBonds.size());
    assert(ppatch_idx < m_PatchyBonds[p].size());
    assert(p < m_PatchyBonds.size());
    assert(qpatch_idx < m_PatchyBonds[q].size());
    m_PatchyBonds[p][ppatch_idx] = {q,qpatch_idx};
    m_PatchyBonds[q][qpatch_idx] = {p,ppatch_idx};
}

void RaspberryInteraction::clear_bound_to(int p, int ppatch_idx) {
    assert(p < m_PatchyBonds.size());
    assert(ppatch_idx < m_PatchyBonds[p].size());
    m_PatchyBonds[p][ppatch_idx] = {-1,-1};
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