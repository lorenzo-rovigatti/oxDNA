/**
 * @file    Morse.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef Morse_H_
#define Morse_H_

#include "BaseForce.h"

class Morse: public BaseForce {
private:
	int _particle;
	int _ref_id;

public:
	BaseParticle * _p_ptr;
	number _r0;
	bool PBC;
	number _a;
	number _D;

	Morse();
	virtual ~Morse() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
    void graph_data();

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // Morse_H
