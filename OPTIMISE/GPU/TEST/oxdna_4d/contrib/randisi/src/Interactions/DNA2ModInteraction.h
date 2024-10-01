/**
 * @brief Adds some customizability features to the oxDNA2 potential represented in the DNA2Interaction
 *
 * Ferdinando Nov 2017
 * 
 * To use this interaction, set
 *
 * interaction_type = DNA2_mod 
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
[mod_nucleotide_group_<N> = <nucleotide_index_list> (defaults to none) (list of nucleotides belonging to some group (currently only N=1 and N=2 are supported). The list has the indices separated by commas or dashes as in the rest of the code, such that a list 10,11-13 translates to nucleotides 10, 11, 12, 13, as elsewhere in the oxDNA code).
[mod_hb_multiplier_<N> = <float> (defaults to 1) (HB interaction multiplier applied to all nucleotides in group N. Notice that this is in addition to any other multiplier applied by the argument hb_multiplier (see DNAInteraction.h))]
[mod_stacking_roll_<N> = <float> (defaults to 0) (when computing the stacking interaction between all nucleotides of group N and their 5' neighbour, the stack vector will be rotated by this angle (in degrees) around the base vector, so that it should introduce a roll.)]
[mod_stacking_tilt_<N> = <float> (defaults to 0) (when computing the stacking interaction between all nucleotides of group N and their 5' neighbour, the stack vector will be rotated by this angle (in degrees) around the vector normal to the base vector and the stacking vector, so that it should introduce a tilt.)]
[mod_stacking_multiplier_<N> = <float> (defaults to 0) (when computing the stacking interaction between all nucleotides of group N and their 5' neighbour, the energy, forces and torques are going to be multiplied by this factor. Defaults to 1.)]

 @endverbatim
 */

#ifndef DNA2MOD_INTERACTION_H
#define DNA2MOD_INTERACTION_H

#include "DNAInteraction.h"
#include "DNA2Interaction.h"

template<typename number>
class DNA2ModInteraction: public DNA2Interaction<number> {
	protected:
		number *_a_hb_multiplier, *_a_stacking_roll, *_a_stacking_r_roll, *_a_stacking_tilt, *_a_stacking_multiplier;
// Horrible copy-paste that should be substituted with a vector as needed
	private:
		std::vector<int> _vec_group_1;
		number _hb_multiplier_1;
		number _mod_stacking_roll_1;
		number _mod_stacking_r_roll_1;
		number _mod_stacking_tilt_1;
		number _mod_stacking_multiplier_1;

		std::vector<int> _vec_group_2;
		number _hb_multiplier_2;
		number _mod_stacking_roll_2;
		number _mod_stacking_r_roll_2;
		number _mod_stacking_tilt_2;
		number _mod_stacking_multiplier_2;

		static std::vector<int> unsafeGetParticlesFromString(std::string particle_string, char const * identifier);
		virtual bool _is_particle_in_group(BaseParticle<number> *particle, std::vector<int> group){
			return std::find(group.begin(), group.end(), particle->index) != group.end();
		}
		virtual bool _is_integer_in_group(int integer, std::vector<int> group){
			return std::find(group.begin(), group.end(), integer) != group.end();
		}
		// get the matrix that rotates the vectors in order to change the stacking
		inline LR_matrix<number> _get_rotation_matrix(number const roll, number const tilt){
			// first we introduce roll, so we rotate around the hb vector (first component)
			LR_matrix<number> R_roll(1.,0.,0., 0.,cos(roll),-sin(roll), 0.,sin(roll),cos(roll));
			//then we introduce the tilt, so that we rotate around the second vector
			LR_matrix<number> R_tilt(cos(tilt),0.,-sin(tilt), 0.,1.,0., sin(tilt),0.,cos(tilt));
			return R_roll * R_tilt;
			
		}


	public:
		DNA2ModInteraction();
		virtual ~DNA2ModInteraction();
		virtual void get_settings(input_file &inp);
		virtual inline number get_hb_multiplier(BaseParticle<number> *p, BaseParticle<number> *q){ 
			number hb_multiplier = (abs(p->btype) >= 300 && abs(q->btype) >= 300) ? this->_hb_multiplier : 1.f; 
			bool p_in_group_1 = _is_particle_in_group(p, _vec_group_1);
			bool q_in_group_1 = _is_particle_in_group(q, _vec_group_1);

			bool p_in_group_2 = _is_particle_in_group(p, _vec_group_2);
			bool q_in_group_2 = _is_particle_in_group(q, _vec_group_2);

			if (p_in_group_1 or q_in_group_1) hb_multiplier *= _hb_multiplier_1;
			else if(p_in_group_2 or q_in_group_2) hb_multiplier *= _hb_multiplier_2;
			
			return hb_multiplier;
		}
		const virtual inline LR_vector<number> get_rotated_T_v3(BaseParticle<number> *p){
			bool p_in_group_1 = _is_particle_in_group(p, _vec_group_1);
			bool p_in_group_2 = _is_particle_in_group(p, _vec_group_2);
			number angle_1 = p_in_group_1 ? _mod_stacking_roll_1 : (p_in_group_2 ? _mod_stacking_roll_2 : 0);
			number angle_2 = p_in_group_1 ? _mod_stacking_tilt_1 : (p_in_group_2 ? _mod_stacking_tilt_2 : 0);
			LR_vector<number> rotated = rotateVectorAroundVersor(p->orientationT.v3, p->orientationT.v1, angle_1);
			if ((p_in_group_1 or p_in_group_2) and false){ 
				printf("angle_1 %g, angle_2 %g\n",angle_1, angle_2);
				printf("rotated.v3 = %g %g %g\n",rotated.x, rotated.y, rotated.z);
				rotated = rotateVectorAroundVersor(rotated, p->orientationT.v2, angle_2);
				printf("rotated.v3 = %g %g %g\n",rotated.x, rotated.y, rotated.z);
			}
			return rotated;
	}
	
	LR_vector<number> rotateVectorAroundVersor(const LR_vector<number> vector, const LR_vector<number> versor, const number angle);
	LR_matrix<number> rotationMatrixAroundVersorByAngle(const LR_vector<number> vector, const number angle);

	protected:
	virtual number _stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	virtual number _hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);



};
	// Rotate a vector around a versor (TODO: propose to add this to defs.h)
template<typename number>
LR_vector<number> DNA2ModInteraction<number>::rotateVectorAroundVersor(const LR_vector<number> vector, const LR_vector<number> versor, const number angle){
	/* According to section 5.2 of this webpage http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ ,
	   the image of rotating a vector (x,y,z) around a versor (u,v,w) is given by
	        / u(ux + vy + wz)(1 - cos th) + x cos th + (- wy + vz) sin th \
	        | v(ux + vy + wz)(1 - cos th) + y cos th + (+ wx - uz) sin th |
	        \ w(ux + vy + wz)(1 - cos th) + z cos th + (- vx + uy) sin th /
	   Note that the first parenthesis in every component contains the scalar product of the vectors,
	   and the last parenthesys in every component contains a component of the cross product versor ^ vector.
	   The derivation that they show on the webpage seems sound, and I've checked it with Mathematica in about
     1 hour of uninterrupted swearing. */
	number costh = cos(angle);
	number sinth = sin(angle);
	number scalar = vector*versor;
	LR_vector<number> cross = versor.cross(vector);
	return LR_vector<number>( 
		versor.x * scalar * (1. - costh) + vector.x*costh + cross.x * sinth,
		versor.y * scalar * (1. - costh) + vector.y*costh + cross.y * sinth,
		versor.z * scalar * (1. - costh) + vector.z*costh + cross.z * sinth
		  );  
}

	//Create the rotation matrix around a vector
template<typename number>
LR_matrix<number> DNA2ModInteraction<number>::rotationMatrixAroundVersorByAngle(const LR_vector<number> versor, const number angle){
	//Implementing Rodrigues rotation formula, as described e.g. on https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
	LR_matrix<number> K(0.,-versor.z,versor.y, versor.z,0.,-versor.x, -versor.y,versor.x,0.);
	LR_matrix<number> identity(1.,0.,0., 0.,1.,0., 0.,0.,1.);
	
	return identity + K * sin(angle) + K * K * (1.-cos(angle));
}

extern "C" BaseInteraction<float> *make_DNA2ModInteraction_float();
extern "C" BaseInteraction<double> *make_DNA2ModInteraction_double(); 
#endif /* DNA2MOD_INTERACTION_H */
