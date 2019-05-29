/**
 * @file    MCTras.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#ifndef MCTRAS_H_
#define MCTRAS_H_

#include "BaseMove.h"

template<typename number>
class MCTras : public BaseMove<number> {
	protected:
		number _delta;
		LR_vector<number> pos_old;

		number _verlet_skin;
	
	public:
		MCTras(ConfigInfo<number> *Info);
		virtual ~MCTras();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
};
#endif // MCTras.h
