/*
 * JordanOutput.h
 *
 *  Created on: 23/mar/2016
 *      Author: Flavio
 */

#ifndef JORDANOUTPUT_H_
#define JORDANOUTPUT_H_

#include "../Configurations/Configuration.h"

/**
 * @brief Prints a file with a line per each particle in the system. The line contains position (x, y, z) and the position of each patch.
 *
 * example section to add to an input file:
 <pre>
 data_output_1 = {
 name = caca.jordan
 only_last = false
 print_every = 1e4
 col_1 = {
 type = jordan_conf
 }
 }
 </pre>
 *
 *
 * this will print on the file "caca.jordan" the configuration every 10^4
 * simulation steps.
 * Since only_last is set to true, every 10^4 steps, the file will be augmented
 */

class JordanOutput: public Configuration {
protected:
	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle *p);
	virtual std::string _configuration(llint step);

public:
	JordanOutput();
	virtual ~JordanOutput();

	void get_settings(input_file &my_inp, input_file &sim_inp);
};

#endif /* JORDANOUTPUT_H_ */
