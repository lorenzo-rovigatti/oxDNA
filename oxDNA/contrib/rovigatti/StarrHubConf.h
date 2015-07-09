/*
 * StarrHubConf.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef STARRCONF_H_
#define STARRCONF_H_

#include "Configurations/Configuration.h"

/**
 * @brief Prints a configuration in the Starr format.
 *
 * @verbatim
[in_box = <bool> (if true all the positions are brought back between -L/2 and L/2. Defaults to false)]
@endverbatim
 */
template<typename number>
class StarrHubConf: public Configuration<number>  {
protected:
	int _N_per_strand;
	int _N_strands_per_tetramer;
	int _N_per_tetramer;
	int _N_tetramers;

	bool _print_bonds;
	number _dt;

	std::vector<LR_vector<number> > _tetra_poss;
	std::vector<LR_vector<number> > _tetra_vels;

	std::vector<std::map<int, int> > _tetra_bonds;

	virtual std::string _headers(llint step);
	std::string _configuration(llint step);

public:
	StarrHubConf();
	virtual ~StarrHubConf();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_StarrHubConf_float() { return new StarrHubConf<float>(); }
extern "C" BaseObservable<double> *make_StarrHubConf_double() { return new StarrHubConf<double>(); }

#endif /* STARRCONF_H_ */
