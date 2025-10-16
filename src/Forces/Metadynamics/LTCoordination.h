/*
 * LTCoordination.h
 *
 * Created on: 10/16/2025
 *      Author: Lorenzo
*/

#ifndef LTCOORDINATION_H
#define LTCOORDINATION_H

#include "../BaseForce.h"

class LTCoordination: public BaseForce {
public:
    LTCoordination();

    ~LTCoordination();

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

    std::vector<BaseParticle *> p1_ptr;
    std::vector<BaseParticle *> p2_ptr;

    number xmin;
	number xmax;
	int N_grid;
	number dX;
    std::vector<number> potential_grid;

private:
    number _coordination();
};

#endif /* LTCOORDINATION_H */
