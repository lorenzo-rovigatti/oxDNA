/**
 * @file    FFS_MD_CPUBackend.h
 * @date    03/set/2013
 * @author  sulc
 *
 *
 */

#ifndef FFS_MD_CPUBACKEND_H_
#define FFS_MD_CPUBACKEND_H_

#include "MD_CPUBackend.h"
#include "../Utilities/OrderParameters.h"

#include <vector>
#include <algorithm>

#define MAX_BOND_CUTOFF (0.f)

/**
 * @brief This class implements an MD algorithm, and stops if conditions are met. The conditions are specified
 * in a file (say ffs.txt, and one has to put ffs_file = ffs.txt into input file)
 * The format of the ffs file is the following, for example:
 * condition1 = {
 *   native_bonds == 4
 *   all_bonds >= 4
 * }
 * condition2 = {
 *  all_bonds <= 0
 * }
 *
 * The program stops when condition1 OR condition2 are satisfied
 * For a conditiond to be satisified, all subcondistions enclosed in {} have to be true
 * Optionally, you can choose to use custom cutoff in the definition of order parametersfor hydrogen bonds, for example, the order paramter file can look like this:
 *  {
 *   order_parameter = bond
 *   name = hb1
 *   cutoff = 0
 *   pair1 = 0, 607
 *   pair2 = 1, 606
 *  }
 *  The cutoff is by default -0.1 as defined HB_CUTOFF in OrderParameters.h
 *  The maximum value you can use is 0.
 *  Note that cutoff option is relevant only for FFS backend. It is ignored by all other Backends and Observables
 *
 @verbatim
 backend = CPU (For CPU FFS)
 sim_type = FFS_MD (This must be set for an FFS simulation)
 @endverbatim
 */

struct parsed_expression {
public:
	/// allowed symbols: >=, <= , < >, ==
	enum {
		GEQ = 0, LEQ = 1, LT = 2, GT = 3, EQ = 4
	};

	/// 1 for distance, 0 for hydrogen bonds
	int expression_type;
	/// GEQ, LT, etc..
	int compare_type;

	int value_to_compare_with_hb;
	double value_to_compare_with_dist;

	int parameter_index;

	parsed_expression() :
					expression_type(-1) {
	}
	bool parse_expression(const char *expression, OrderParameters *op);

	parsed_expression(const char *expression, OrderParameters *op) {
		parse_expression(expression, op);
	}

	bool eval_expression(OrderParameters *op);

};

struct parsed_condition {
	std::vector<parsed_expression> all_expressions;

	bool parse_condition(const char *expression, OrderParameters *op);

	parsed_condition(void) {
	}
	parsed_condition(const char *expression, OrderParameters *op) {
		parse_condition(expression, op);
	}
	bool eval_condition(OrderParameters *op);
};

/**
 * Example of input file options for FFS:
 * order_parameters_file = op.txt
 * ffs_file = ffs.txt
 *
 * example of ffs_file:
 * {
 * action = stop_and
 * condition1 = params_a > 1
 * condition2 = params_b >= 4
 * condition3 = params_c < 4
 * }
 *
 */

class FFS_MD_CPUBackend: public MD_CPUBackend {
protected:

	OrderParameters _op;
	std::vector<parsed_condition> _conditions;
	std::string _order_parameters_file;
	std::string _ffs_file;
	char _state_str[2048];

	number _sqr_rcut;
	void _ffs_compute_forces(void);
	number pair_interaction_nonbonded_DNA_with_op(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces = false);

public:
	FFS_MD_CPUBackend();
	virtual ~FFS_MD_CPUBackend();

	void print_order_parameters(void);
	void init_ffs_from_file(const char *fname);
	bool check_stop_conditions(void);

	virtual void get_settings(input_file &inp);

	void init();

	void sim_step();
	char * get_op_state_str(void);
	virtual void print_observables();

};

#endif /* FFS_MD_CPUBACKEND_H_ */
