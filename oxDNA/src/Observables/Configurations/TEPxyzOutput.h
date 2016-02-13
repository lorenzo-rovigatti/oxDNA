/*
 * TEPxyzOutput.h
 *
 *  Created on: 19/jun/2015
 *      Author: Ferdinando (after TcLOutput.h, written by Flavio)
 */

#ifndef TEPXYZOUTPUT_H_
#define TEPXYZOUTPUT_H_

#include "../Configurations/Configuration.h"
#include "../../Particles/TEPParticle.h"

/**
 * @brief Prints a xyz configuration for TEP systems. The configuration can be loaded from the command 
 * line in VMD with source "filename.xyz"
 *
 * The supported syntax is (optional values are between [])
@verbatim
[back_in_box = <bool>  (Default: true; if true the particle positions will be brought back in the box)]
[show = <int>,<int>,... (Default: all particles; list of comma-separated indexes of the particles that will be shown. Other particles will not appear)]
[hide = <int>,<int>,... (Default: no particles; list of comma-separated indexes of particles that will not be shown)]
[ref_particle = <int> (Default: -1, no action; The particle with the id specified (starting from 0) is set at the centre of the box. Overriden if ref_strands is specified. Ignored if negative or too large for the system.)]
[ref_strand = <int> (Default: -1, no action; The strand with the id specified, starting from 1, is set at the centre of the box. Ignored if negative or too large for the system.)]
@endverbatim
 * Note that you cannot put both the 'show' and 'hide' keys in the input file.
 * If you do so, the 'hide' key will not be considered.
 *
 * example section to add to an input file:
<pre>
data_output_1 = {
name = gino.xyz
only_last = true 
print_every = 1e4
col_1 = {
  type = xyz_configuration
  ref_strand = 1
  print_labels = true
  }
}
</pre>
 *
 *
 * this will print on the file "gino.xyz" the configuration every 10^4
 * simulation steps, printing out the strand labels and centering the box on
 * the center of mass of strand 1
 */
template<typename number>
class TEPxyzOutput: public Configuration<number>  {
protected:
	bool _back_in_box;
	bool _print_labels;
	int _ref_particle_id;
	int _ref_strand_id;
	int _resolution;
	int *_bead_types;
	
	int _weak_bead_1_index;
	int _weak_bead_2_index;
	std::string _particles_type1_string;
	std::string _particles_type2_string;
	number _core_radius, _side_radius, _front_radius, _backback_radius;
	number _side_shift, _front_shift;
	
	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);

public:
	TEPxyzOutput();
	virtual ~TEPxyzOutput();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init(ConfigInfo<number> &config_info);
};

#endif /* TEPXYZOUTPUT_H_ */
