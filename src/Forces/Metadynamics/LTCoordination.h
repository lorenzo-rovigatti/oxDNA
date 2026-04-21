/*
 * LTCoordination.h
 *
 * Created on: 10/16/2025
 *      Author: Lorenzo
*/

#ifndef LTCOORDINATION_H
#define LTCOORDINATION_H

#include "../BaseForce.h"
#include "meta_utils.h"

#include <unordered_map>
#include <string>

class LTCoordination: public BaseForce {
public:
    LTCoordination();

    virtual ~LTCoordination();

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector force(llint step, LR_vector &pos) override;
    LR_vector torque(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

    // all pairs of particles involved in the coordination number calculation
    // (i.e. all pairs of particles that are part of the same hydrogen-bond-based order parameter defined in the OP file)
    std::vector<std::pair<BaseParticle *, BaseParticle *>> all_pairs;

    // mapping from particle pointer to the index of the pair it belongs to in the all_pairs vector
    std::unordered_map<BaseParticle *, int> particle_to_pair_index;

    number coord_min; // minimum coordination number, defaults to 0.0
	number coord_max; // maximum coordination number, defaults to the number of pairs defined in the OP file
	int N_grid;
	number d_coord; // grid spacing
    number delta_F = 1e-6; // displacement for finite difference calculation of the force
    number delta_T = 1e-6; // displacement for finite difference calculation of the torque

    meta::CoordSettings settings;
    std::vector<number> potential_grid;

private:
    LR_vector _dcoord_dpos(std::pair<BaseParticle*, BaseParticle*> &pair);
    number _dcoord_dr(std::pair<BaseParticle*, BaseParticle*> &pair);
    LR_vector _dcoord_dtheta();

    LR_matrix _rot_matrices[3][2]; // rotation matrices for the 3 possible rotation axes and 2 possible directions of rotation (positive or negative)
};

#endif /* LTCOORDINATION_H */
