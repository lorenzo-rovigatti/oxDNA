/*
 * ChimeraOutput.h
 *
 *  Created on: 18/jan/2014
 *      Author: C. Matek
 */

#ifndef CHIMERAOUTPUT_H_
#define CHIMERAOUTPUT_H_

#include "../Configurations/Configuration.h"
#include "../../Particles/DNANucleotide.h"

/**
 * @brief Prints the commands that are needed to obtain the nice representation in UCSF chimera
 * for DNA systems
 *
 * this observable is selected with type = chimera_script
 *
 @verbatim
 [colour_by_sequece = <bool> (Default: false; whether to coulour the bases according to the base type (A, C, G, T)]
 @endverbatim
 *
 */

class ChimeraOutput: public Configuration {
protected:
	bool _colour_by_seq;

	virtual std::string _headers(llint step);
	virtual std::string _configuration(llint step);

public:
	ChimeraOutput();
	virtual ~ChimeraOutput();

	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* CHIMERAOUTPUT_H_ */
