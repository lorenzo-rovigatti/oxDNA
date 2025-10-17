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

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

    std::vector<std::pair<BaseParticle *, BaseParticle *>> all_pairs;

    // unfortunately forces in oxDNA are passed vectors rather than particles,
    // so we need to keep a map from positions to pair of particles
    std::unordered_map<LR_vector *, int> pos_to_pair_index;

    number coord_min; // minimum coordination number, defaults to 0.0
	number coord_max; // maximum coordination number, defaults to the number of pairs defined in the OP file
	int N_grid;
	number d_coord; // grid spacing
    number d0; // optimal distance for coordination, defaults to 1.2 in internal units
    number r0; // width of the switching function, defaults to 1.2 in internal units
    int n; // exponent of the switching function, defaults to 6
    std::vector<number> potential_grid;

private:
    number _coordination();
};

#endif /* LTCOORDINATION_H */
