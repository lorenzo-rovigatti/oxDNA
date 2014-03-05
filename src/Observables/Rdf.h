/*
 * Rdf.h
 *
 *  Created on: Nov 21, 2013
 *      Author: Flavio
 */

#ifndef RDF_H_
#define RDF_H_

#include "BaseObservable.h"
#include <sstream>

/**
 * @brief Outputs the radial distribution function g(r).
 *
 * This is a classic observable from liquid simulations. The g(r) is normalized
 * in the standard way (goes to 1 at large r).
 * It can be made 1- or 2-dimensional by specifying the axes. If they are not
 * specified, a 3-D system is assumed.
 *
 * the parameters are as follows:
@verbatim
max_value = <float> (maximum r to consider)
bin_size = <float> (bin size for the g(r))
[axes = <string> (Possible values: x, y, z, xy, yx, zy, yz, xz, zx. Those are the axes to consider in the computation. Mind that the normalization always assumes 3D sytems for the time being.)]
@endverbatim
 *
 * example input file section:
<pre>
analysis_data_output_1 = {
name = gdr.dat
print_every = 1
only_last = yes
col_1 = {
type = rdf
max_value = 6.
bin_size = 0.02
\#axes = xy # if enabled, this would ignore the z coordinate in the computation
}
}
</pre>
 *
 *
 */
template<typename number>
class Rdf : public BaseObservable<number> {
private:
	long int _nconf;
	number _max_value;
	number _bin_size;
	int _nbins;
	std::vector<long double> _profile;
	LR_vector<number> _mask;
public:
	Rdf();
	virtual ~Rdf();

	virtual std::string get_output_string(llint curr_step);
	void get_settings (input_file &my_inp, input_file &sim_inp);
};

#endif /* RDF_H_ */
