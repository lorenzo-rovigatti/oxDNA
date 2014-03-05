/*
 * tcl.h
 *
 *  Created on: 25/ott/2013
 *      Author: Flavio
 */

#ifndef TCLOUTPUT_H_
#define TCLOUTPUT_H_

#include "../Configurations/Configuration.h"
#include "../../Particles/DNANucleotide.h"

/**
 * @brief Prints a tcl configuration for DNA systems. It still is not M/M grooving aware. The configuration can be loaded from the command 
 * line in VMD with source "filename.tcl"
 *
 * The supported syntax is (optional values are between [])
@verbatim
[back_in_box = <bool>  (Default: true; if true the particle positions will be brought back in the box)]
[show = <int>,<int>,... (Default: all particles; list of comma-separated indexes of the particles that will be shown. Other particles will not appear)]
[hide = <int>,<int>,... (Default: no particles; list of comma-separated indexes of particles that will not be shown)]
[print_labels = <bool> (Default: false; if true labels with the strand id are printed next to one end of the strand.)]
[resolution = <int> (Default: 20; resolution set in the tcl file.)]
[ref_particle = <int> (Default: -1, no action; The nucleotide with the id specified (starting from 0) is set at the centre of the box. Overriden if ref_strands is specified. Ignored if negative or too large for the system.)]
[ref_strand = <int> (Default: -1, no action; The strand with the id specified, starting from 1, is set at the centre of the box. Ignored if negative or too large for the system.)]
@endverbatim
 * Note that you cannot put both the 'show' and 'hide' keys in the input file.
 * If you do so, the 'hide' key will not be considered.
 *
 * example section to add to an input file:
<pre>
data_output_1 = {
name = gino.tcl
only_last = true 
print_every = 1e4
col_1 = {
  type = tcl_configuration
  ref_strand = 1
  print_labels = true
  }
}
</pre>
 *
 *
 * this will print on the file "gino.tcl" the configuration every 10^4
 * simulation steps, printing out the strand labels and centering the box on
 * the center of mass of strand 1
 */
template<typename number>
class TclOutput: public Configuration<number>  {
protected:
	bool _back_in_box;
	bool _print_labels;
	int _ref_particle_id;
	int _ref_strand_id;
	int _resolution;
	number _back_radius, _base_radius, _backbase_radius, _backback_radius;
	
	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle<number> *p);
	virtual std::string _configuration(llint step);

public:
	TclOutput();
	virtual ~TclOutput();

	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* TCLOUTPUT_H_ */
