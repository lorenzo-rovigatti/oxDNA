/*
 * SaltExtrapolation.h
 *
 *  Created on: Jul 11, 2014
 *      Author: Flavio
 */

#ifndef SALTEXTRAPOLATION_H_
#define SALTEXTRAPOLATION_H_

#include "BaseObservable.h"
#include "../Interactions/InteractionFactory.h"

#include "../Utilities/OrderParameters.h"
#include "../Utilities/Histogram.h"
#include "../Utilities/Weights.h"

/**
 * @brief Observable that prints out a histogram extrapolated with a proper treatment of the temperature dependence in the Debye potential and also extrapolates to different salt conditions
 *
 * @verbatim
salts = <float>, <float>, ... (list of salt concentration to extrapolate to)
temps = <T>, <T>, ... (list of temperatures to extrapolate to, separated with commas. Temperatures can be specified in reduced units, Kelvin, Celsius as 0.10105, 30C, 30c, 30 c, 303.15 k, 303.15K, 303.15k)
[op_file = <string> (order parameter file. If not found, it will use the one from the input file)]
[weights_file = <string> (weights file. If not found, the one from the input file will be used.)]
@endverbatim
 */
template<typename number>
class SaltExtrapolation: public BaseObservable<number> {
protected:
	std::vector<number> _salts;
	std::vector<number> _temps;
	std::vector<std::vector< IBaseInteraction<number> * > > _interactions;
	IBaseInteraction<number> * _ref_interaction;
	OrderParameters _op;
	Weights _weights;
	std::string _op_file, _weights_file;
	std::vector<std::vector <std::vector <long double> > > _hists;
	number _sim_temp;
	bool _skip_zeros;

public:
	SaltExtrapolation();
	virtual ~SaltExtrapolation();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	void init (ConfigInfo<number> &Info);
};

#endif /* SALTEXTRAPOLATION_H_ */

