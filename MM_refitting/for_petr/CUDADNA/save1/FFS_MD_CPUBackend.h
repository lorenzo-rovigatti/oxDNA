/*
 * MD_CPUBackend.h
 *
 *  Created on: 03/set/2010
 *      Author: sulc
 */

#ifndef FFS_MD_CPUBACKEND_H_
#define FFS_MD_CPUBACKEND_H_

#include "MD_CPUBackend.h"
#include "OrderParameters.h"

#include <vector>
#include <algorithm>

#define BOND_CUTOFF (-0.1f)

//this class is responsible for processing an OP

struct parsed_expression {
public:
	//allowed symbols: >=, <= , < >, ==
	enum {
		GEQ = 0, LEQ = 1, LT = 2, GT = 3, EQ = 4
	};

	int expression_type; //1 for distance, 0 for hydrogen bonds
	//int expression_operator;
	int compare_type; // GEQ, LT, etc..

	int value_to_compare_with_hb;
	double value_to_compare_with_dist;

	int parameter_index;

	parsed_expression() {
		expression_type = -1;
	}
	bool parse_expression(const char *expression, OrderParameters *op);

	parsed_expression(const char *expression, OrderParameters *op) {
		parse_expression(expression, op);
	}

	bool eval_expression(OrderParameters *op);

};

/*
 * Example of input file options for FFS:
 order_parameters_file = op.txt
 ffs_file = ffs.txt

 example of ffs_file:
 {
 action = stop_and
 condition1 = params_a > 1
 condition2 = params_b >= 4
 condition3 = params_c < 4
 }

 */
template<typename number>
class FFS_MD_CPUBackend: public MD_CPUBackend<number> {
protected:

	OrderParameters _op;
	std::vector<parsed_expression> _conditions;
	int _ffs_type;
	char _order_parameters_file[512];
	char _ffs_file[512];
	char _state_str[1024];

	number _ffs_particle_particle_interaction(Particle<number> *p,
			Particle<number> *q);
	void _ffs_compute_forces(void);

public:
	FFS_MD_CPUBackend(IOManager *IO);
	virtual ~FFS_MD_CPUBackend();

	void print_order_parameters(void);
	void init_ffs_from_file(const char *fname);
	bool check_stop_conditions(void);

	virtual void get_settings(input_file &inp) {
		MD_CPUBackend<number>::get_settings(inp);
		getInputString(&inp, "order_parameters_file", _order_parameters_file,
				1);
		getInputString(&inp, "ffs_file", _ffs_file, 1);
	}

	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);
	
	void sim_step(llint cur_step);
	char * get_op_state_str(void);
	virtual void print_energy(llint curr_step);

};

#endif /* FFS_MD_CPUBACKEND_H_ */
