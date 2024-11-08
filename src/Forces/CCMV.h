/**
 * @file    CCMV.h
 * @date    10/06/24
 * @author  Manuel Micheloni
 *
 */

#ifndef CCMV_H_
#define CCMV_H_

#include "BaseForce.h"

class CCMV: public BaseForce {
public:
	/// center of the sphere
	LR_vector _center = LR_vector(0., 0., 0.);

	number _radius;
	
    number _A_Vittorio = -7.77913022e-02;
	number _B_Vittorio = 4.59994536;
	number _C_Vittorio = -1.08939804e+01;
	number _D_Vittorio = 7.68423122;
	number _E_Vittorio = -2.79192908;
	number _F_Vittorio = 6.06976204e-01;
	number _G_Vittorio = -8.29137297e-02;
	number _H_Vittorio = 7.15232653e-03;
	number _I_Vittorio = -3.77388959e-04;
    number _L_Vittorio = 1.11053305e-05;
	number _M_Vittorio = -1.39522295e-07;

	number _cutoff;
	number _V_amplt;

	CCMV();
	virtual ~CCMV() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // CCMV_H_
