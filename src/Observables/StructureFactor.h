/*
 * StructureFactor.h
 *
 *  Created on: Nov 21, 2013
 *      Author: Flavio
 */

#ifndef STRUCTUREFACTOR_H_
#define STRUCTUREFACTOR_H_

#include "BaseObservable.h"
#include <sstream>
#include <list>

/**
 * @brief Outputs the structure factor S(q).
 *
 * This is a classic observable from liquid simulations.
 *
 * the parameters are as follows:
@verbatim
max_q = <float> (maximum q to consider)
[type = <int> (particle species to consider. Defaults to -1, which means "all particles")]
@endverbatim
 *
 * example input file section:
<pre>
analysis_data_output_1 = {
name = sq.dat
print_every = 1
only_last = yes
col_1 = {
type = Sq
max_q = 6.
}
}
</pre>
 *
 *
 */
template<typename number>
class StructureFactor : public BaseObservable<number> {
private:
	long int _nconf;
	number _max_q;
	std::vector<long double> _sq, _sq_cos, _sq_sin;
	std::list<LR_vector<double> > _qs;
	int _type;
	int _max_qs_in_interval;
	number _max_qs_delta;
	bool _always_reset;

public:
	StructureFactor();
	virtual ~StructureFactor();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
	virtual std::string get_output_string(llint curr_step);
};

#endif /* STRUCTUREFACTOR_H_ */
