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

//Small value for tolerance when comparing norms (squares of lengths) of projection vectors with 0.0 (if 'truly' zero, i.e. vector perpendicular to projection plane, then of the order of <E-28)
#define TOLERANCE_PROJ 1.0E-10
//Small value for tolerance when comparing: norms (squares of lengths) of cross products of projection vectors with 0.0 (if 'truly' zero, i.e. projection vectors parallel/antiparallel (or projection vectors of length zero, but we already ruled this out), then of the order of <E-28); dot product of projection vectors with 0.0 (don't know magnitude if 'truly' zero, i.e. projection vectors perpendicular (or projection vectors of length zero, but we already ruled this out), because we will have already showed projection vectors parallel/antiparallel)
#define TOLERANCE_PROD 1.0E-16

//Small value cutoff for G_ij or G_ji; if angle_r_p or angle_r_q is such that G_ij or G_ji is less than the cutoff (for the given sigma_ang), we can avoid some expensive calculations; for sigma_ang = 0.3, corresponds to angle_r_p or angle_r_q of 127.8 degrees
#define G_CUTOFF 1.0E-12
//Could do the same for large G_ij or G_ji, but (for sigma_ang = 0.3) even for G = 0.999999, angle_r_p or angle_r_q is 0.024 degrees, which will occur rarely, so this will not give a significant saving
//Could do the same for small or large V_tor, but (for sigma_tor = 0.6) even for V_tor = 0.999999, angle_tor is 0.00085 degrees, which will occur rarely, so this will not give a significant saving; and even for angle_tor = +/-180 degrees, V_tor = 1.1E-6, which is not negligible

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

class PatchyInteractionDan: public BaseInteraction {
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
	LR_vector **_patch_vectors_type;
	//Patch vector of each patch on each particle
	LR_vector **_patch_vectors_particle;
	//Patch type of each patch on each particle type
	int **_patch_type_of;
	//sigma_ang for each patch type
	number *_sigma_ang_patch;
	//epsilon for each patch type interacting with each patch type
	number **_epsilon_patch;

	//Reference vector for each patch on each particle type
	LR_vector **_ref_vectors_type;
	//Reference vector for each patch on each particle
	LR_vector **_ref_vectors_particle;
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

	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r - from p to q
	 * @param update_forces - always false for MC
	 * @return
	 */

	inline number _patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

public:
	//DT - bond type (irrelevant for me as all bonds are the same type; I use 4 = hydrogen bond as this allows me to use Observable/HBList & HBEnergy)
	enum {
		PATCHY = 4
	};

	// Constructor and destructor
	PatchyInteractionDan();
	virtual ~PatchyInteractionDan();

	// All in main PatchyInteractionDan.cpp file
	virtual void get_settings(input_file &inp);
	virtual void init();
	virtual void read_topology(int *N_strands, std::vector<BaseParticle*> &particles);
	virtual void allocate_particles(std::vector<BaseParticle*> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	/*16-06-07virtual void generate_random_configuration(std::vector<BaseParticle *> &particles, number box_side);*/
	virtual void check_input_sanity(std::vector<BaseParticle*> &particles);

};

number PatchyInteractionDan::_patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.f;
	//Square of length of vector between particle centres (LR_vector->norm is actually square of real norm)
	number sqr_r_dist = _computed_r.norm();

	//If interparticle distance is greater than cut-off, potential is zero
	if(sqr_r_dist >= _sqr_rcut) {
		return energy;

	}

	//Length of vector between particle centres
	number r_dist = sqrt(sqr_r_dist);
	//EVAprintf("%.16lf, ", r_dist);

	//For computational efficiency, calculate part of LJ potential between particles
	//(sigma_LJ would be numerator; power of 3 (not 6) because sqr_r_dist is length squared)
	number V_LJ_part = pow((1 / sqr_r_dist), 3);
	//LJ potential between particles; epsilon would be * 4
	number V_LJ = 4 * (SQR(V_LJ_part) - V_LJ_part);
	//Shift potential based on cutoff (cut-and-shift)
	V_LJ -= _LJ_cut;

	//Piecemeal potential; length units of sigma_LJ
	if(r_dist < _r_crossover) {
		energy = V_LJ;

		return energy;

	}

	//Variables (explained when used)
	number V_tor;

	//Current maximum of patch-patch modulation to LJ potential (now: start at 0.0 for the case where one/both particles have no patches; formerly: start at -1.0 rather than 0.0 to make sure that first patch pair is checked)
	number max_V_ang_V_tor_epsilon = (number) 0.0;
	//FOR COMPARISONnumber max_V_tor = (number) 0.0;

	//Recast BaseParticle's as PatchyParticleDan's, so can access members specific to PatchyParticleDan (in particular, _ref_vectors)
	PatchyParticleDan *pp = (PatchyParticleDan*) p;
	PatchyParticleDan *qq = (PatchyParticleDan*) q;

	//Calculate cutoff values for cos of angle_r_p and angle_r_q, for each patch type; if angle_r_p or angle_r_q is such that G_ij or G_ji is less than cos_ang_cutoffs, we can avoid some expensive calculations
	number min2_ln_Gcutoff = -2.0 * log(G_CUTOFF);
	number *cos_ang_cutoffs = new number[_N_patch_types];
	for(int patch_type = 0; patch_type < _N_patch_types; patch_type++) {
		cos_ang_cutoffs[patch_type] = cos(sqrt( SQR(_sigma_ang_patch[patch_type]) * min2_ln_Gcutoff));
	}

	//Iterate through patches on q, and calculate G_ji term in patch-patch interaction (V_ang_q)
	number *V_ang_qs = new number[qq->N_int_centers()];
	for(uint q_patch = 0; q_patch < qq->N_int_centers(); q_patch++) {
		//Define patch vector; qpatch is vector to patch from centre of particle
		LR_vector qpatch = qq->int_centers[q_patch];
		//Calculate cosine of angle between REVERSE of r and patch vector on q (r is from p to q); second '*' in numerator is dot product
		number cos_angle_r_q = -(_computed_r * qpatch) / (r_dist * qpatch.module());
		//If cos_angle_r_q is is less than cos_angle_cutoffs (for the given sigma_ang), we can avoid expensive acos and exp calculations
		if(cos_angle_r_q < cos_ang_cutoffs[_patch_type_of[qq->type][q_patch]]) {
			V_ang_qs[q_patch] = (number) 0.0;
		}
		else {
			number angle_r_q = LRACOS(cos_angle_r_q);
			V_ang_qs[q_patch] = exp(-SQR(angle_r_q) / (2.0 * SQR(_sigma_ang_patch[_patch_type_of[qq->type][q_patch]])));
		}
	}

	//Iterate through patches on p, performing operations which relate only to p
	for(uint p_patch = 0; p_patch < pp->N_int_centers(); p_patch++) {

		//Define patch vector and reference vector (ppatch is vector to patch from centre of particle, pref is reference vector for patch), and projection of particle reference vector onto the plane perpendicular to the interparticle vector, r ('*' between vectors is dot product)
		LR_vector ppatch = pp->int_centers[p_patch];
		LR_vector pref;
		LR_vector proj_pref;
		if(_tor_flag == true) {
			pref = pp->_ref_vectors[p_patch];
			proj_pref = pref - (_computed_r * ((pref * _computed_r) / sqr_r_dist));
		}

		//As for q_patches above
		number V_ang_p;
		number cos_angle_r_p = (_computed_r * ppatch) / (r_dist * ppatch.module());
		if(cos_angle_r_p < cos_ang_cutoffs[_patch_type_of[pp->type][p_patch]]) {
			V_ang_p = (number) 0.0;
		}
		else {
			number angle_r_p = LRACOS(cos_angle_r_p);
			V_ang_p = exp(-SQR(angle_r_p) / (2.0 * SQR(_sigma_ang_patch[_patch_type_of[pp->type][p_patch]])));
		}

		//Iterate through patches on q
		for(uint q_patch = 0; q_patch < qq->N_int_centers(); q_patch++) {
			//EVAprintf(", , %d, %d, ", p_patch, q_patch);
			//If epsilon_patch is less than or equal to the current max_V_ang_V_tor_epsilon, we cannot exceed max_V_ang_V_tor_epsilon (as V_ang and V_tor are <= 1) for this pair of patches, so we skip to the next patch pair (next iteration)
			if(_epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]] <= max_V_ang_V_tor_epsilon) {

				/*printf("PID.h (3) _epsilon_patch[_patch_type_of[pp->type %d][p_patch %d] %d][_patch_type_of[qq->type %d][q_patch %d] %d] (%f) < max_V_ang_V_tor_epsilon (%f); energy %f\n", pp->type, p_patch, _patch_type_of[pp->type][p_patch], qq->type, q_patch, _patch_type_of[qq->type][q_patch], _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]], max_V_ang_V_tor_epsilon, energy);
				 fflush(stdout);*/

				continue;

			}

			//Calculate V_ang and V_ang_epsilon
			number V_ang = V_ang_p * V_ang_qs[q_patch];
			number V_ang_epsilon = V_ang * _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]];
			//If V_ang_epsilon is less than the current max_V_ang_V_tor_epsilon, we cannot exceed max_V_ang_V_tor_epsilon (as V_tor <= 1) for this pair of patches, so we skip to the next pair (iteration)
			if(V_ang_epsilon < max_V_ang_V_tor_epsilon) {
				continue;

			}

			if(_tor_flag == true) {
				//Current minimum of square of angle_diff (difference between angle_tor and offset angle), across equivalent offset angles
				number min_sqr_angle_diff = M_PI * M_PI;

				/*//Current maximum of patch-patch torsional modulation to LJ potential (maximum across equivalent offset angles)
				 number max_V_tor_offset = (number) 0.0;*/

				//Check if one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r). If so, print a warning statement, and assign a default value to V_tor. If a reference vector is parallel or antiparallel with r, it is likely the patch vector is roughly perpendicular with r (because the reference vectors should be roughly perpendicular with the patch vectors). Thus, it is likely the patch-patch interaction is weak (V_ang will be small)
				//Note that if one or both of the projected vectors is zero, their cross product and dot product will also be zero. Note also that if one or both of the projected vectors is zero, we cannot calculate angle_tor using the method below
				//Previously, we skipped this pair of patches; this is unsatisfactory if there is only one pair of patches, or the interaction between all other patch pairs is very weak. Given, as stated above, it is likely V_ang is small, the default V_tor should not have a significant effect
				if(proj_pref.norm() < TOLERANCE_PROJ) {

					//if ((abs(proj_pref.x) < TOLERANCE) && (abs(proj_pref.y) < TOLERANCE) && (abs(proj_pref.z) < TOLERANCE)) {
					//throw oxDNAException("Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because former projected vector is zero. Aborting\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);
					/*Full warning statement (reduced below):
					 printf("WARNING: Cross product (%f %f %f) of {projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle p (%d) onto the plane perpendicular to r (%f %f %f), with projection (%f %f %f) of reference vector (%f %f %f) of patch %d on particle q (%d) onto same plane} is zero, because former projected vector is zero.\nSkipping this pair of patches.\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);*/

					//printf("WARNING: Projection of reference vector onto r is 0. [Projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle pp (%d), for r (%.16lf %.16lf %.16lf).]\nSkipping this pair of patches.\n", proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, r->x, r->y, r->z);
					//count_pref++;
					//EVAprintf("\n");
					V_tor = V_TOR_PARALLEL;

					//continue;

				}
				else {

					//Define reference vector for patch on q, and projection of particle reference vector onto the plane perpendicular to the interparticle vector, r ('*' between vectors is dot product)
					LR_vector qref = qq->_ref_vectors[q_patch];
					LR_vector proj_qref = qref - (_computed_r * ((qref * _computed_r) / sqr_r_dist));

					/*printf(" qref .x %f .y %f .z %f\n", qref.x, qref.y, qref.z);
					 printf(" proj_qref .x %f .y %f .z %f\n", proj_qref.x, proj_qref.y, proj_qref.z);*/

					if(proj_qref.norm() < TOLERANCE_PROJ) {
						V_tor = V_TOR_PARALLEL;

						//continue;

					}
					else {
						number angle_tor;

						//Cross product of these projected vectors
						LR_vector cross_proj_pref_proj_qref = proj_pref.cross(proj_qref);

						//If the cross product is zero, we cannot determine whether it is parallel or antiparallel with r (which we need to do to calculate angle_tor). So, we work through the possible cases
						if(cross_proj_pref_proj_qref.norm() < TOLERANCE_PROD) {
							//To begin, calculate dot product of projected vectors
							number dot_proj_pref_proj_qref = proj_pref * proj_qref;

							//Case 1: cross product is zero because projected vectors are parallel (hence dot product is positive)
							//In this case, the angle between the projected vectors is zero
							if(dot_proj_pref_proj_qref > TOLERANCE_PROD) {

								angle_tor = 0.0;
								//Case 2: cross product is zero because projected vectors are antiparallel (hence dot product is negative)
								//In this case, the angle between the projected vectors is pi
							}
							else if(dot_proj_pref_proj_qref < -TOLERANCE_PROD) {
								angle_tor = M_PI;

								//Case 3: cross product is zero, but not because projected vectors are parallel or antiparallel, and also not because one or both of the projected vectors is zero (we checked this earlier). Thus, cross product is zero for unknown reason, and so - for now - we abort.

								/*Old code
								 //Case 3a: cross product is zero because one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r)
								 //In this case, we print a warning statement, and skip this pair of patches. If a reference vector is parallel or antiparallel with r, it is likely the patch vector is roughly perpendicular with r (because the reference vectors should be roughly perpendicular with the patch vectors). Thus, it is likely the patch-patch interaction is weak (V_ang will be small).

								 //Case 3b: cross product is zero, but not because one or both of the projected vectors is zero (which occurs when one or both reference vectors is parallel or antiparallel with r). The reason is unknown, and so - for now - we abort.
								 */

							}
							else {
								throw oxDNAException("Cross product (%.16lf %.16lf %.16lf) of {projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle pp (%d) onto the plane perpendicular to r (%.16lf %.16lf %.16lf), with projection (%.16lf %.16lf %.16lf) of reference vector (%.16lf %.16lf %.16lf) of patch %d on particle qq (%d) onto same plane} is zero, for unknown reason. Aborting\n", cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, pref.x, pref.y, pref.z, p_patch, pp->index, _computed_r.x, _computed_r.y, _computed_r.z, proj_qref.x, proj_qref.y, proj_qref.z, qref.x, qref.y, qref.z, q_patch, qq->index);
							}

						}
						else {

							//Angle between projected vectors (numerator is dot product)
							angle_tor = LRACOS((proj_pref * proj_qref) / (proj_qref.module() * proj_pref.module()));

							//Cross product dotted with interparticle vector
							number sign = cross_proj_pref_proj_qref.x * _computed_r.x + cross_proj_pref_proj_qref.y * _computed_r.y + cross_proj_pref_proj_qref.z * _computed_r.z;

							//Adjust angle_tor (if necessary), so it is between 0 and 360 degrees, measured clockwise (positive) from proj_pref when looking along r
							if(sign > 0.0) {
								//Do nothing
							}
							else if(sign < 0.0) {
								//Adjust
								angle_tor = (2.0 * M_PI) - angle_tor;
							}
							else {
								throw oxDNAException("Sign (sign of %.16lf) of dot product of vector r (%.16lf, %.16lf, %.16lf) (from particle pp %d to particle qq %d) with vector (%.16lf, %.16lf, %.16lf) that is cross product of projection (onto plane perpendicular to r) of reference vector (%.16lf, %.16lf, %.16lf) of patch %d on p and projection of reference vector (%.16lf, %.16lf, %.16lf) of patch %d on qq is exactly 0.0. Reason is unknown: should have been caught earlier, when checked if cross product was zero. Aborting\n", sign, _computed_r.x, _computed_r.y, _computed_r.z, pp->index, qq->index, cross_proj_pref_proj_qref.x, cross_proj_pref_proj_qref.y, cross_proj_pref_proj_qref.z, proj_pref.x, proj_pref.y, proj_pref.z, p_patch, proj_qref.x, proj_qref.y, proj_qref.z, q_patch);
							}
						}

						//Try all possible offset angles for each pair of patches; choose that which gives the maximum V_tor (minimum sqr_angle_diff)
						for(int offset_angle = 0; offset_angle < _number_offset_angles[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]; offset_angle++) {

							//Difference between angle_tor and offset angle
							number angle_diff = angle_tor - _offset_angle_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]][offset_angle];
							//Adjust difference so that it is between -PI and PI
							if(angle_diff < -M_PI) {
								angle_diff = angle_diff + 2.0 * M_PI;
							}
							else if(angle_diff > M_PI) {
								angle_diff = angle_diff - 2.0 * M_PI;
							}
							//Square of angle_diff
							number sqr_angle_diff = SQR(angle_diff);

							if(sqr_angle_diff < min_sqr_angle_diff)
								min_sqr_angle_diff = sqr_angle_diff;

						}

						//Calculate V_tor term in patch-patch interaction
						V_tor = exp(-min_sqr_angle_diff / (2.0 * SQR(_sigma_tor_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]])));
					}
				}

				//If tor_flag = 0, essentially ignore V_tor term
			}
			else {

				V_tor = 1.0;

			}

			//Calculate patch-patch modulation term
			number V_ang_V_tor_epsilon = V_ang_epsilon * V_tor;

			//printf("PIh5 _epsilon_patch[_patch_type_of %d[pp->type %d][p_patch %d][_patch_type_of %d[qq->type %d][q_patch %d]] %f\n", _patch_type_of[pp->type][p_patch], pp->type, p_patch, _patch_type_of[qq->type][q_patch], qq->type, q_patch, _epsilon_patch[_patch_type_of[pp->type][p_patch]][_patch_type_of[qq->type][q_patch]]);

			//If current patch-patch interaction is stronger than all others considered, record it as strongest interaction
			if(V_ang_V_tor_epsilon > max_V_ang_V_tor_epsilon) {

				max_V_ang_V_tor_epsilon = V_ang_V_tor_epsilon;
			}
		}
	}

	energy = V_LJ * max_V_ang_V_tor_epsilon;

	delete[] V_ang_qs;
	delete[] cos_ang_cutoffs;

	return energy;
}
;

#endif /* PATCHYINTERACTIONDAN_H_ */
