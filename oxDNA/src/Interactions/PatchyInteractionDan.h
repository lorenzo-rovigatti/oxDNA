/*
 * PatchyInteractionDan.h
 *
 *  Created on: 25/jan/2016
 *      Author: lorenzo -> dan
 */

#ifndef PATCHYINTERACTIONDAN_H_
#define PATCHYINTERACTIONDAN_H_

#include "BaseInteraction.h"
#include "../Particles/PatchyParticleDan.h"

//Small value for tolerance when comparing: squares of lengths of projection vectors with 0.0 (if 'truly' zero, i.e. vector perpendicular to projection plane, then of the order of <E-28); squares of lengths of cross products of projection vectors with 0.0 (if 'truly' zero, i.e. projection vectors parallel/antiparallel (or length zero, but we already ruled this out), then of the order of <E-28); dot product of projection vectors with 0.0 (don't know magnitude if 'truly' zero, i.e. projection vectors perpendicular (or length zero, but we already ruled this out), because we already showed projection vectors parallel/antiparallel)
#define TOLERANCE_PROJ 1.0E-16
/*Old code (DT)
//Small value for tolerance when comparing squares of lengths of projection vectors with 0.0 (if 'truly' zero, i.e. vector perpendicular to projection plane, then of the order of <E-28)
//Small value for tolerance when comparing squares of lengths of cross products of projection vectors with 0.0 (if truly zero, then of the order of E-28) 
#define TOLERANCE_PROD 1.0E-16*/

//Default value of _rcut (in BaseInteraction.h)
#define DEFAULT_RCUT 3.0
//number _default_rcut = (number) 3.0;
//Default value of _tor_flag (0 = OFF, 1 = ON)
#define DEFAULT_TOR_FLAG 0

//Value assigned to V_tor if one/both patch's reference vector is parallel or antiparallel with interparticle vector
#define V_TOR_PARALLEL 0.5

/**
@brief Manages the interaction between simple patchy particles (as described in http://scitation.aip.org/content/aip/journal/jcp/131/17/10.1063/1.3243581)
 *
 * This interaction is selected with
 * interaction_type = patchyDan
 *
 */

template <typename number>
class PatchyInteractionDan : public BaseInteraction<number, PatchyInteractionDan<number> > {
protected:
        /*//Flag to denote when initialisation of pointers is complete
        bool _initialised;*/

        //Flag denoting whether or not interactions are torsional
        bool _tor_flag;

        //Number of particles
        int _N_particles;
        //Number of particle types
        int _N_particle_types;
	//Number of patch types
        int _N_patch_types;
        //Number of particles of each particle type
        int *_N_particles_of_type;
        //Particle type of each particle
        int *_particle_type_of;
	//Number of patches per particle for each particle type
	int *_N_patches_type;
	//Number of patches on each particle
        int *_patches_on_particle;
	//Patch vector of each patch on each particle type
        LR_vector<number> **_patch_vectors_type;
	//Patch vector of each patch on each particle
        LR_vector<number> **_patch_vectors_particle;
	//Patch type of each patch on each particle type
        int **_patch_type_of;
        //sigma_ang for each patch type
        number *_sigma_ang_patch;
        //epsilon for each patch type interacting with each patch type
        number **_epsilon_patch;

	//Reference vector for each patch on each particle type
        LR_vector<number> **_ref_vectors_type;
	//Reference vector for each patch on each particle
        LR_vector<number> **_ref_vectors_particle;
        //sigma_tor for each patch type interacting with each patch type
        number **_sigma_tor_patch;
	//Number of equivalent offset angles for each patch type-patch type pair
	int **_number_offset_angles;
	//Offset angle(s) for each patch type interacting with each patch type
        number ***_offset_angle_patch;

        //Repulsive LJ interaction at the cut-off
	number _LJ_cut;

	//Piecemeal potential crossover distance
	number _r_crossover;

	//Error-checking: count number of times pref and qref are approximately 0
	//int count_pref, count_qref, count_cross;

	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r - from p to q
	 * @param update_forces - always false for MC
	 * @return
	 */

	inline number _patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);

public:
//DT - bond type (irrelevant for me as all bonds are the same type; I use 4 = hydrogen bond as this allows me to use Observable/HBList & HBEnergy)
	enum {
		PATCHY = 4
	};

//Constructor and destructor
	PatchyInteractionDan();
	virtual ~PatchyInteractionDan();

//All in main PatchyInteractionDan.cpp file
        virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);
	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	/*16-06-07virtual void generate_random_configuration(BaseParticle<number> **particles, int N, number box_side);*/
	virtual void check_input_sanity(BaseParticle<number> **particles, int N);

 };

template<typename number>
number PatchyInteractionDan<number>::_patchy_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
        //printf("PID.h, _patchy_interaction\n");

        /*printf("PID.h (0) p->index (%d), q->index (%d), r-vector %f, %f, %f\n", p->index, q->index, r->x, r->y, r->z);
	  fflush(stdout);*/

	//if(p->index>1) throw oxDNAException("[In PatchyInteractionDan::_patchy_interaction] p->index>0. Aborting\n");

        //Interaction energy; initialise to zero
        number energy = (number) 0.0;

        /*for (int particle = 0; particle < _N_particles; particle++) {
            printf("particle %d, patches %d, particle type %d\n", particle, _patches_on_particle[particle], _particle_type_of[particle]);
          }*/

        //Square of length of vector between particle centres (LR_vector->norm is actually square of real norm)
	number sqr_r_dist = r->norm();

	//EVAprintf("%d, %d, %.16f, %.16f, %.16f, %.16lf, %.16f, %.16f, %.16lf, %.16f, %.16f, %.16lf, ", p->index, q->index, p->pos.x, p->pos.y, p->pos.z, q->pos.x, q->pos.y, q->pos.z, r->x, r->y, r->z, sqr_r_dist);

	//If interparticle distance is greater than cut-off, potential is zero
	if (sqr_r_dist >= this->_sqr_rcut) {
	  //EVAprintf("\n, , %.16lf\n", energy);
	  /*printf("PID.h (1) sqr_r_dist (%f) > this->_sqr_rcut (%f); energy 0.0\n", sqr_r_dist, this->_sqr_rcut);
	    fflush(stdout);*/

	  return energy;

	}

        //Length of vector between particle centres
	number r_dist = sqrt(sqr_r_dist);
	//EVAprintf("%.16lf, ", r_dist);

	//For computational efficiency, calculate part of LJ potential between particles
	//(sigma_LJ would be numerator; power of 3 (not 6) because sqr_r_dist is length squared)
	number V_LJ_part =  pow((1 / sqr_r_dist), 3);
	//LJ potential between particles; epsilon would be * 4
	number V_LJ = 4 * (SQR(V_LJ_part) - V_LJ_part);
	//EVAprintf("%.16lf, ", V_LJ);
	//Shift potential based on cutoff (cut-and-shift)
	V_LJ -= _LJ_cut;
	//EVAprintf("%.16lf\n", V_LJ);

	/*Old code - think it is irrelevant
	//Update forces is always false
	if(update_forces) {
		LR_vector<number> force = *r * (PATCHY_POWER * part / sqr_r_dist);
		p->force -= force;
		q->force += force;
	}

	//??Some sort of counter?
	int c = 0;

	LR_vector<number> tmptorquep(0, 0, 0);
	LR_vector<number> tmptorqueq(0, 0, 0);*/

	//Piecemeal potential; length units of sigma_LJ
	if(r_dist < _r_crossover) {
	  
	  energy = V_LJ;
	  //EVAprintf("\n, , %.16lf\n", energy);
	  /*printf("PID.h (2) r_dist (%f) < _r_crossover (%f); energy %f\n", r_dist, _r_crossover, energy);
	    fflush(stdout);*/

	  return energy;
	  
	}

	//Variables (explained when used)
	number V_ang_p, V_ang_q, V_ang, V_ang_epsilon, V_tor, V_ang_V_tor_epsilon;
	number angle_r_p, angle_r_q;
	/*May be used later
	  int max_patch_index[2] = {-1, -1};*/
	LR_vector<number> ppatch, qpatch, pref, qref;
	
	//Current maximum of patch-patch modulation to LJ potential (now: start at 0.0 for the case where one/both particles have no patches; formerly: start at -1.0 rather than 0.0 to make sure that first patch pair is checked)
	number max_V_ang_V_tor_epsilon = (number) 0.0;
	//FOR COMPARISONnumber max_V_tor = (number) 0.0;
	
	//Recast BaseParticle's as PatchyParticleDan's, so can access members specific to PatchyParticleDan (in particular, _ref_vectors)
	PatchyParticleDan<number> * pp = (PatchyParticleDan<number> *) p;
	PatchyParticleDan<number> * qq = (PatchyParticleDan<number> *) q;

	//Iterate through patches on each particle
	for(int p_patch = 0; p_patch < pp->N_int_centers; p_patch++) {
	  for(int q_patch = 0; q_patch < qq->N_int_centers; q_patch++) {
	    //EVAprintf(", , %d, %d, ", p_patch, q_patch);
	    //If epsilon_patch is less than or equal to the current max_V_ang_V_tor_epsilon, we cannot exceed max_V_ang_V_tor_epsilon (as V_ang and V_tor are <= 1) for this pair of patches, so we skip to the next patch pair (next iteration)
	    if (_epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]] <= max_V_ang_V_tor_epsilon) {
	      
	      /*printf("PID.h (3) _epsilon_patch[_patch_type_of[pp->type %d][p_patch %d] %d][_patch_type_of[qq->type %d][q_patch %d] %d] (%f) < max_V_ang_V_tor_epsilon (%f); energy %f\n", pp->type, p_patch, _patch_type_of[pp->type][p_patch], qq->type, q_patch, _patch_type_of[qq->type][q_patch], _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]], max_V_ang_V_tor_epsilon, energy);
		fflush(stdout);*/
	      
	      continue;
	      
	    }

	    //Define patch vectors and reference vectors; ppatch and qpatch are vectors to patches from centres of particles, pref and qref are reference vectors for patches
	    ppatch = pp->int_centers[p_patch];
	    qpatch = qq->int_centers[q_patch];
	    pref = pp->_ref_vectors[p_patch];
	    qref = qq->_ref_vectors[q_patch];
	    //EVAprintf("%.16lf, %.16f, %.16f, %.16f, %.16f, %.16f, ", ppatch.x, ppatch.y, ppatch.z, qpatch.x, qpatch.y, qpatch.z);
	    //EVAprintf("%.16lf, %.16f, %.16f, %.16f, %.16f, %.16f, ", pref.x, pref.y, pref.z, qref.x, qref.y, qref.z);
	    
	    /*Old code//Vector between patches, from p's patch to q's patch
	      interpatch_vector = *r + qpatch - ppatch;*/
	    
	    //Calculate angle between r and patch vector on p; second '*' in numerator is dot product
	    angle_r_p = LRACOS((*r * ppatch) / (r_dist * ppatch.module()));
	    //Calculate angle between REVERSE of r and patch vector on q (r is from p to q)
	    angle_r_q = LRACOS(((*r * -1.0) * qpatch) / (r_dist * qpatch.module()));
	    //EVAprintf("%.16lf, %.16f, ", angle_r_p, angle_r_q);
	    
	    //Calculate G_ij and G_ji terms in patch-patch interaction
	    V_ang_p = exp(-SQR(angle_r_p) / (2.0 * SQR(_sigma_ang_patch[_patch_type_of[pp->type][p_patch]])));
	    V_ang_q = exp(-SQR(angle_r_q) / (2.0 * SQR(_sigma_ang_patch[_patch_type_of[qq->type][q_patch]])));
	    //EVAprintf("%.16lf, %.16f, ", V_ang_p, V_ang_q);
	    
	    //Calculate V_ang and V_ang_epsilon
	    V_ang = V_ang_p * V_ang_q;
	    V_ang_epsilon = V_ang * _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]];
	    //EVAprintf("%.16lf, ", V_ang);
	    
	    /*Old code (DT)
	    //Calculate V_ang_epsilon
	    V_ang_V_tor_epsilon = V_ang_p * V_ang_q * _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]];
	    */
	    
	    /*printf("PID.h (4):");
	    //printf(" pp->index %d\n", pp->index);
	    //printf(" pp->type %d\n", pp->type);
	    //printf(" pp->N_int_centers %d\n", pp->N_int_centers);
	    printf(" p_patch %d\n", p_patch);
	    //printf(" _patch_type_of[pp->type][p_patch] %d\n", _patch_type_of[pp->type][p_patch]);
	    printf(" pp->_base_patch_vectors[p_patch] .x %f .y %f .z %f\n", pp->_base_patch_vectors[p_patch].x, pp->_base_patch_vectors[p_patch].y, pp->_base_patch_vectors[p_patch].z);
	    printf(" pp->int_centers[p_patch] .x %f .y %f .z %f\n", pp->int_centers[p_patch].x, pp->int_centers[p_patch].y, pp->int_centers[p_patch].z);
	    //printf(" ppatch .x %f .y %f .z %f\n", ppatch.x, ppatch.y, ppatch.z);
	    //printf(" pp->_base_ref_vectors[p_patch] .x %f .y %f .z %f\n", pp->_base_ref_vectors[p_patch].x, pp->_base_ref_vectors[p_patch].y, pp->_base_ref_vectors[p_patch].z);
	    //printf(" pp->_ref_vectors[p_patch] .x %f .y %f .z %f\n", pp->_ref_vectors[p_patch].x, pp->_ref_vectors[p_patch].y, pp->_ref_vectors[p_patch].z);
	    //printf(" pref .x %f .y %f .z %f\n", pref.x, pref.y, pref.z);
	    //printf(" qq->index %d\n", qq->index);
	    //printf(" qq->type %d\n", qq->type);
	    //printf(" qq->N_int_centers %d\n", qq->N_int_centers);
	    //printf(" q_patch %d\n", q_patch);
	    //printf(" _patch_type_of[qq->type][q_patch] %d\n", _patch_type_of[qq->type][q_patch]);
	    //printf(" qq->int_centers[q_patch] .x %f .y %f .z %f\n", qq->int_centers[q_patch].x, qq->int_centers[q_patch].y, qq->int_centers[q_patch].z);
	    //printf(" qpatch .x %f .y %f .z %f\n", qpatch.x, qpatch.y, qpatch.z);
	    //printf(" qq->_ref_vectors[q_patch] .x %f .y %f .z %f\n", qq->_ref_vectors[q_patch].x, qq->_ref_vectors[q_patch].y, qq->_ref_vectors[q_patch].z);
	    //printf(" qref .x %f .y %f .z %f\n", qref.x, qref.y, qref.z);
	    fflush(stdout);*/
	    
	    //printf("PIh1 _sigma_ang_patch[_patch_type_of %d[pp->type %d][p_patch %d] %f\n", _patch_type_of[pp->type][p_patch], pp->type, p_patch, _sigma_ang_patch[_patch_type_of[pp->type][p_patch]]);
	    //printf("PIh2 _sigma_ang_patch[_patch_type_of %d[qq->type %d][q_patch %d] %f\n", _patch_type_of[qq->type][q_patch], qq->type, q_patch, _sigma_ang_patch[_patch_type_of[qq->type][q_patch]]);
	    
	    //If V_ang_epsilon is less than the current max_V_ang_V_tor_epsilon, we cannot exceed max_V_ang_V_tor_epsilon (as V_tor <= 1) for this pair of patches, so we skip to the next pair (iteration)
	    if (V_ang_epsilon < max_V_ang_V_tor_epsilon) {
	      
	      /*printf("PID.h (5) V_ang_V_tor_epsilon (%f) < max_V_ang_V_tor_epsilon (%f)\n", V_ang_V_tor_epsilon, max_V_ang_V_tor_epsilon);
		fflush(stdout);*/
	      
	      continue;
	      
	    }
	    
	    if (_tor_flag == true) {
	      
	      //Variables (explained when used)
	      number dot_proj_pref_proj_qref, sign, angle_diff, sqr_angle_diff, angle_tor; // V_tor_offset;
	      LR_vector<number> proj_pref, proj_qref, cross_proj_pref_proj_qref;
	      
	      //Current minimum of square of angle_diff (difference between angle_tor and offset angle), across equivalent offset angles
	      number min_sqr_angle_diff = M_PI*M_PI;
	      
	      /*//Current maximum of patch-patch torsional modulation to LJ potential (maximum across equivalent offset angles)
		number max_V_tor_offset = (number) 0.0;*/
	      
	      /*printf("PID.h (6) _tor_flag==true:");
		printf(" pp->index %d\n", pp->index);
		printf(" pp->type %d\n", pp->type);
		printf(" p_patch %d\n", p_patch);
		printf(" _patch_type_of[pp->type][p_patch] %d\n", _patch_type_of[pp->type][p_patch]);
		printf(" qq->index %d\n", qq->index);
		printf(" qq->type %d\n", qq->type);
		printf(" q_patch %d\n", q_patch);
		printf(" _patch_type_of[qq->type][q_patch] %d\n", _patch_type_of[qq->type][q_patch]);
		fflush(stdout);*/
	      
	      //Projections of particle reference vectors onto the plane perpendicular to the interparticle vector, r ('*' between vectors is dot product)
	      proj_pref = pref - (*r * ((pref * *r) / sqr_r_dist));
	      
	      //As they are calculated, check if one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r). If so, print a warning statement, and assign a default value to V_tor. If a reference vector is parallel or antiparallel with r, it is likely the patch vector is roughly perpendicular with r (because the reference vectors should be roughly perpendicular with the patch vectors). Thus, it is likely the patch-patch interaction is weak (V_ang will be small)
	      //Note that if one or both of the projected vectors is zero, their cross product and dot product will also be zero. Note also that if one or both of the projected vectors is zero, we cannot calculate angle_tor using the method below
	      //Previously, we skipped this pair of patches; this is unsatisfactory if there is only one pair of patches, or the interaction between all other patch pairs is very weak. Given, as stated above, it is likely V_ang is small, the default V_tor should not have a significant effect
	      
	      if (proj_pref.norm() < TOLERANCE_PROJ) {
		
		//if ((abs(proj_pref.x) < TOLERANCE) && (abs(proj_pref.y) < TOLERANCE) && (abs(proj_pref.z) < TOLERANCE)) {
		//throw oxDNAException("Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because former projected vector is zero. Aborting\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);
		/*Full warning statement (reduced below):
		  printf("WARNING: Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because former projected vector is zero.\nSkipping this pair of patches.\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);*/
		
		//printf("WARNING: Projection of reference vector onto r is 0. [Projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle pp (%d), for r (%.16lf %.16lf %.16lf).]\nSkipping this pair of patches.\n", proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z);
		
		//count_pref++;
		//EVAprintf("\n");
		
		V_tor = V_TOR_PARALLEL;

		//continue;
		
	      } else {
		
		proj_qref = qref - (*r * ((qref * *r) / sqr_r_dist));
		
		if (proj_qref.norm() < TOLERANCE_PROJ) {
		  
		  //} else if ((abs(proj_qref.x) < TOLERANCE) && (abs(proj_qref.y) < TOLERANCE) && (abs(proj_qref.z) < TOLERANCE)) {
		  //throw oxDNAException("Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because latter projected vector is zero. Aborting\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);
		  
		  /*Full warning statement (reduced below):
		    printf("WARNING: Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because latter projected vector is zero.\nSkipping this pair of patches.\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);*/
		  
		  //printf("WARNING: Projection of reference vector onto r is 0. [Projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle qq (%d), for r (%.16lf %.16lf %.16lf).]\nSkipping this pair of patches.\n", proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index, r->x, r->y, r->z);
		
		  //count_qref++;
		  //EVAprintf("\n");
		  
		  V_tor = V_TOR_PARALLEL;

		  //continue;
		  
		} else {
		  
		  //If neither of the projected vectors is zero, proceed
		  
		  //EVAprintf("%.16lf, %.16f, %.16f, %.16f, %.16f, %.16f, ", proj_pref.x, proj_pref.y, proj_pref.z, proj_qref.x, proj_qref.y, proj_qref.z);
		  
		  /*printf("PID.h (7) proj_pref.x (%.16lf) .y (%.16lf) .z (%.16lf), proj_qref.x (%.16lf) .y (%.16lf) .z (%.16lf)\n", proj_pref.x, proj_pref.y, proj_pref.z, proj_qref.x, proj_qref.y, proj_qref.z);
		    fflush(stdout);*/
		  
		  //Angle between these projected vectors (numerator is dot product)
		  angle_tor = LRACOS((proj_pref * proj_qref) / (proj_qref.module() * proj_pref.module()));
		  //EVAprintf("%.16lf, ", angle_tor);
		  
		  /*printf("PID.h (8) angle_tor (%f)\n", angle_tor);
		    fflush(stdout);*/
		  
		  //Cross product of these projected vectors
		  cross_proj_pref_proj_qref = proj_pref.cross(proj_qref);
		  
		  /*printf("PID.h (9) cross_proj_pref_proj_qref.x (%.16lf), .y (%.16lf), .z (%.16lf)\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z);
		    fflush(stdout);*/
		  
		  //If the cross product is zero, we cannot determine whether it is parallel or antiparallel with r (which we need to do to calculate angle_tor). So, we work through the possible cases
		  if (cross_proj_pref_proj_qref.norm() < TOLERANCE_PROJ) {
		    //if ((abs(cross_proj_pref_proj_qref.x) < TOLERANCE) && (abs(cross_proj_pref_proj_qref.y) < TOLERANCE) && (abs(cross_proj_pref_proj_qref.z) < TOLERANCE)) {
		    
		    /*printf("PID.h (10) (abs(cross_proj_pref_proj_qref.x) (%.16lf) < TOLERANCE) && (abs(cross_proj_pref_proj_qref.y) (%.16lf) < TOLERANCE) && (abs(cross_proj_pref_proj_qref.z) (%.16lf) < TOLERANCE)\n", abs(cross_proj_pref_proj_qref.x), abs(cross_proj_pref_proj_qref.y), abs(cross_proj_pref_proj_qref.z));
		      fflush(stdout);*/
		    
		    //count_cross++;
		    
		    //To begin, calculate dot product of projected vectors
		    dot_proj_pref_proj_qref = proj_pref * proj_qref;
		    
		    //Case 1: cross product is zero because projected vectors are parallel (hence dot product is positive)
		    //In this case, the angle between the projected vectors is zero
		    if (dot_proj_pref_proj_qref > TOLERANCE_PROJ) {
		      /*printf("PID.h (11) dot_proj_pref_proj_qref (%.16lf) > TOLERANCE)\n", dot_proj_pref_proj_qref);
			fflush(stdout);*/
		      
		      angle_tor = 0.0;
		      //EVAprintf("L, G, n/a, ");
		      //printf("cross_proj_pref_proj_qref = 0 and dot_proj_pref_proj_qref > 0, for vector r (%f, %f, %f), from particle p %d to particle q %d, cross product (%f, %f, %f), projection (onto plane perpendicular to r) of reference vector (%f, %f, %f) of patch %d on p, and projection of reference vector (%f, %f, %f) of patch %d on q.\n", r->x, r->y, r->z, pp->index, qq->index, cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, p_patch, proj_qref.x, proj_qref.y, proj_qref.z, q_patch);
		      
		      //Case 2: cross product is zero because projected vectors are antiparallel (hence dot product is negative)
		      //In this case, the angle between the projected vectors is pi
		    } else if (dot_proj_pref_proj_qref < -TOLERANCE_PROJ) {
		      /*printf("PID.h (12) dot_proj_pref_proj_qref (%.16lf) < -TOLERANCE)\n", dot_proj_pref_proj_qref);
			fflush(stdout);*/
		      
		      angle_tor = M_PI;
		      //EVAprintf("L, L, n/a, ");
		      
		      //Case 3: cross product is zero, but not because projected vectors are parallel or antiparallel, and also not because one or both of the projected vectors is zero (we checked this earlier). Thus, cross product is zero for unknown reason, and so - for now - we abort.
		      
		      /*Old code
		      //Case 3a: cross product is zero because one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r)
		      //In this case, we print a warning statement, and skip this pair of patches. If a reference vector is parallel or antiparallel with r, it is likely the patch vector is roughly perpendicular with r (because the reference vectors should be roughly perpendicular with the patch vectors). Thus, it is likely the patch-patch interaction is weak (V_ang will be small).
		      
		      //Case 3b: cross product is zero, but not because one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r). The reason is unknown, and so - for now - we abort.
		      */
		      
		    } else {
		      
		      /*printf("PID.h (13) TOLERANCE > dot_proj_pref_proj_qref (%.16lf) > -TOLERANCE)\n", dot_proj_pref_proj_qref);
			printf(" pp->type %d\n", pp->type);
			printf(" pp->N_int_centers %d\n", pp->N_int_centers);
			printf(" pp->index %d\n", pp->index);
			printf(" p_patch %d\n", p_patch);
			printf(" _patch_type_of[pp->type][p_patch] %d\n", _patch_type_of[pp->type][p_patch]);
			printf(" qq->type %d\n", qq->type);
			printf(" qq->N_int_centers %d\n", qq->N_int_centers);
			printf(" qq->index %d\n", qq->index);
			printf(" q_patch %d\n", q_patch);
			printf(" _patch_type_of[qq->type][q_patch] %d\n", _patch_type_of[qq->type][q_patch]);
			fflush(stdout);*/
		      
		      //printf("count_pref %d, count_qref %d, count_cross %d\n", count_pref, count_qref, count_cross);
		      
		      throw oxDNAException("Cross product (%.16lf %.16lf %.16lf) of {projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle pp (%d) onto the plane perpendicular to r (%.16lf %.16lf %.16lf), with projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle qq (%d) onto same plane} is zero, for unknown reason. Aborting\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);
		      
		    }
		    
		  } else {
		    
		    //Cross product dotted with interparticle vector
		    sign = cross_proj_pref_proj_qref.x * r->x + cross_proj_pref_proj_qref.y * r->y + cross_proj_pref_proj_qref.z * r->z;
		    
		    //Adjust angle_tor (if necessary), so it is between 0 and 360 degrees, measured clockwise (positive) from proj_pref when looking along r
		    if (sign > 0.0) {
		      //Do nothing
		      //EVAprintf("G, n/a, G, ");
		    } else if (sign < 0.0) {
		      //Adjust
		      angle_tor = (2.0*M_PI) - angle_tor;
		      //EVAprintf("G, n/a, L, ");
		    } else {
		      throw oxDNAException("Sign (sign of %.16lf) of dot product of vector r (%.16lf, %.16lf, %.16lf) (from particle pp %d to particle qq %d) with vector (%.16lf, %.16lf, %.16lf) that is cross product of projection (onto plane perpendicular to r) of reference vector (%.16lf, %.16lf, %.16lf) of patch %d on p and projection of reference vector (%.16lf, %.16lf, %.16lf) of patch %d on qq is exactly 0.0. Reason is unknown: should have been caught earlier, when checked if cross product was zero. Aborting\n", sign, r->x, r->y, r->z, pp->index, qq->index, cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, p_patch, proj_qref.x, proj_qref.y, proj_qref.z, q_patch);
		    }
		    
		  }
		  //EVAprintf("%.16lf, ", angle_tor);
		  
		  //Try all possible offset angles for each pair of patches; choose that which gives the maximum V_tor (minimum sqr_angle_diff)
		  for (int offset_angle = 0; offset_angle < _number_offset_angles[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]; offset_angle++) {
		    
		    //Difference between angle_tor and offset angle 
		    angle_diff = angle_tor - _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]][offset_angle];
		    //Adjust difference so that it is between -PI and PI
		    if (angle_diff < -M_PI) {
		      angle_diff = angle_diff + 2.0*M_PI;
		    } else if (angle_diff > M_PI) {
		      angle_diff = angle_diff - 2.0*M_PI;
		    }
		    //EVAprintf("%.16lf, %.16lf, ", _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]][offset_angle], angle_diff);
		    
		    /*Old code (DT)
		    //Calculate V_tor term in patch-patch interaction
		    V_tor_offset = exp(-SQR(angle_diff) / (2.0 * SQR(_sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]])));
		    
		    if (V_tor_offset > max_V_tor_offset) max_V_tor_offset = V_tor_offset;
		    
		    }
		    
		    V_tor = max_V_tor_offset;
		    */
		    
		    //Square of angle_diff
		    sqr_angle_diff = SQR(angle_diff);
		    
		    if (sqr_angle_diff < min_sqr_angle_diff) min_sqr_angle_diff = sqr_angle_diff;
		    
		  }
		  
		  //Calculate V_tor term in patch-patch interaction
		  V_tor = exp(-min_sqr_angle_diff / (2.0 * SQR(_sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]])));
		  //EVAprintf("%.16lf, %.16lf, ", min_sqr_angle_diff, V_tor);
		  
		  //FOR COMPARISONif (V_tor > max_V_tor) max_V_tor = V_tor;
		  
		  /*Old code (DT)
		  //Calculate V_tor term in patch-patch interaction
		  //Numerator is square of angle between projected vectors, measured clockwise from proj_pref when looking along r (clockwise +, anticlockwise -)
		  if (sign > 0.0) {
		  
		  V_tor = exp(-SQR(angle_tor - _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]) / (2 * SQR(_sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]])));
		  
		  } else if (sign < 0.0) {
		  
		  V_tor = exp(-SQR((2*M_PI -angle_tor) - _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]) / (2 * SQR(_sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]])));
		  
		  } else {
		  
		  throw oxDNAException("Sign (sign of %f) of dot product of vector r (%f, %f, %f) (from particle p %d to particle q %d) with vector (%f, %f, %f) that is cross product of projection (onto plane perpendicular to r) of reference vector (%f, %f, %f) of patch %d on p and projection of reference vector (%f, %f, %f) of patch %d on q is 0.0. Aborting\n", sign, r->x, r->y, r->z, pp->index, qq->index, cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, p_patch, proj_qref.x, proj_qref.y, proj_qref.z, q_patch);
		  
		  }*/
		  
		  //printf("PIh3 _sigma_tor_patch[_patch_type_of %d[pp->type %d][p_patch %d][_patch_type_of %d[qq->type %d][q_patch %d]] %f\n", _patch_type_of[pp->type][p_patch], pp->type, p_patch, _patch_type_of[qq->type][q_patch], qq->type, q_patch, _sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]);
		  //printf("PIh4 _offset_angle_patch[_patch_type_of %d[pp->type %d][p_patch %d][_patch_type_of %d[qq->type %d][q_patch %d]] %f\n", _patch_type_of[pp->type][p_patch], pp->type, p_patch, _patch_type_of[qq->type][q_patch], qq->type, q_patch, _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]);

		  //Close 'if' statements for 
		}
	      }

	      //If tor_flag = 0, essentially ignore V_tor term
	    } else {
	      
	      V_tor = 1.0;
	      
	    }
	    
	    //Calculate patch-patch modulation term
	    V_ang_V_tor_epsilon = V_ang_epsilon * V_tor;
	    
	    //printf("PIh5 _epsilon_patch[_patch_type_of %d[pp->type %d][p_patch %d][_patch_type_of %d[qq->type %d][q_patch %d]] %f\n", _patch_type_of[pp->type][p_patch], pp->type, p_patch, _patch_type_of[qq->type][q_patch], qq->type, q_patch, _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]);
	    
	    //If current patch-patch interaction is stronger than all others considered, record it as strongest interaction
	    if (V_ang_V_tor_epsilon > max_V_ang_V_tor_epsilon) {
	      
	      max_V_ang_V_tor_epsilon = V_ang_V_tor_epsilon;
	      //EVAprintf("%.16lf, %.16lf\n", V_ang_V_tor_epsilon, max_V_ang_V_tor_epsilon);

	      /*printf("%d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, MAX\n", p->index, q->index, p_patch, q_patch, pp->type, qq->type, _patch_type_of[pp->type][p_patch], _patch_type_of[qq->type][q_patch], p->pos.x, p->pos.y, p->pos.z, q->pos.x, q->pos.y, q->pos.z, r_dist, V_LJ, _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]], angle_r_p, angle_r_q, V_ang, angle_tor, V_tor, V_ang_V_tor_epsilon);
		fflush(stdout);*/

	      /*if (((_patch_type_of[pp->type][p_patch] == 3) && (_patch_type_of[qq->type][q_patch] == 4)) || ((_patch_type_of[pp->type][p_patch] == 4) && (_patch_type_of[qq->type][q_patch] == 3))) {
	      //if ((_patch_type_of[pp->type][p_patch] == 7) && (_patch_type_of[qq->type][q_patch] == 7)) {
		printf("%d, %d, %d, %d, %d, %d, %d, %d, %f, %f\n", p->index, q->index, p_patch, q_patch, pp->type, qq->type, _patch_type_of[pp->type][p_patch], _patch_type_of[qq->type][q_patch], angle_tor, V_ang_V_tor_epsilon);
		}*/
		  
	      /*May be used later
		max_patch_index[0] = p_patch;
		max_patch_index[1] = q_patch;*/

	    }/* else {

	      printf("%d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, --\n", p->index, q->index, p_patch, q_patch, pp->type, qq->type, _patch_type_of[pp->type][p_patch], _patch_type_of[qq->type][q_patch], p->pos.x, p->pos.y, p->pos.z, q->pos.x, q->pos.y, q->pos.z, r_dist, V_LJ, _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]], angle_r_p, angle_r_q, V_ang, angle_tor, V_tor, V_ang_V_tor_epsilon);
	      fflush(stdout);

	      }*/

	  }
	}

	//printf("%d, %d, %.16lf\n", pp->index, qq->index, max_V_tor);

	energy = V_LJ * max_V_ang_V_tor_epsilon;
	//EVAprintf(", %.16lf, %.16lf\n", max_V_ang_V_tor_epsilon, energy);
	
	return energy;
	
};

#endif /* PATCHYINTERACTIONDAN_H_ */
