/*
 * contact_map.h
 *
 * Created on Mar 25, 2019
 *       Author: Poppleton
 */

#ifndef CLOSEST_N_H_
#define CLOSEST_N_H_

#include "Observables/BaseObservable.h"

/**
 * @brief Outputs a contact map for the system
 */

class neighs {
	protected:
		int _max_size;
	public:
		neighs(int max_n);
		~neighs();
		void add_pair(int i, number d);
		std::vector<int> ns;
		std::vector<number> ds;
};


class ClosestN: public BaseObservable {

protected:
	number _rmax;
	int _nmax;

public:
	ClosestN();
	virtual ~ClosestN();
	virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable * make_ClosestN() { return new ClosestN(); }

#endif /*CLOSEST_N_H_ */
