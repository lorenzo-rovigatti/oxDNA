/*
 * LTCoordination.h
 *
 * Created on: 10/16/2025
 *      Author: Lorenzo
*/

#ifndef LTCOORDINATION_H
#define LTCOORDINATION_H

#include "../BaseForce.h"

#include <unordered_map>

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
    number d0 = 0.4; // optimal distance for coordination, defaults to 1.2 in internal units
    number r0 = 0.5; // width of the switching function, defaults to 0.5 in internal units
    int n = 6; // exponent of the switching function, defaults to 6
    std::vector<number> potential_grid;

private:
    LR_vector _distance(std::pair<BaseParticle *, BaseParticle *> &pair);
    number _coordination();
};

#endif /* LTCOORDINATION_H */
