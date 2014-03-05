/*
 * Configuration.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "../BaseObservable.h"
#include "../TotalEnergy.h"

/**
 * @brief Prints ascii configurations. It can be extended to provide further
 * output types (e.g. for visualisation).
 *
 * The supported syntax is (optional values are between [])
 * @verbatim
[back_in_box = <bool> (if true the particle positions will be brought back in the box, defaults to false)]
[show = <int>,<int>,... (list of comma-separated particle indexes whose positions will be put into the final configuration)]
[hide = <int>,<int>,... (list of comma-separated particle indexes whose positions won't be put into the final configuration)]
[reduced = <bool> (if true only the strand centres of mass will be printed, defaults to false)]
 @endverbatim
 * Note that you cannot put both the 'show' and 'hide' keys in the input file.
 * If you do so, the 'hide' key will not be considered.
 */
template<typename number>
class Configuration: public BaseObservable<number>  {
protected:
	bool _back_in_box;
	bool _reduced;
	std::set<int> _visible_particles;
	std::set<int> _hidden_particles;
	std::map<int, LR_vector<number> > _strands_cdm;

	/**
	 * @brief Returns the configuration header(s)
	 *
	 * @param step
	 * @return
	 */
	virtual std::string _headers(llint step);

	/**
	 * @brief Returns the portion of output for the given particle
	 *
	 * @param p
	 * @return
	 */
	virtual std::string _particle(BaseParticle<number> *p);

	/**
	 * @brief Returns the configuration output for the whole system. It does not comprise the headers.
	 *
	 * @param step
	 * @return
	 */
	virtual std::string _configuration(llint step);

	/**
	 * @brief A auxiliary function that fills _strands_cdm, a vector with the c.o.m. of the strands;
	 */
	void _fill_strands_cdm ();

	TotalEnergy<number> _tot_energy;

public:
	Configuration();
	virtual ~Configuration();

	virtual void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
	std::string get_output_string(llint curr_step);
};

#endif /* CONFIGURATION_H_ */

