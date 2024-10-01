/**
 * @file    RotateSite.h
 * @date    23/mar/2016
 * @author  flavio
 *
 */

#ifndef ROTATE_SITE_H_
#define ROTATE_SITE_H_

#include "BaseMove.h"


class RotateSite : public BaseMove {
	protected:
		number _delta;
		std::vector<int> allowed_sites;
		
	public:
		RotateSite();
		virtual ~RotateSite();

		virtual void get_settings(input_file &inp, input_file &sim_inp);
		void apply (llint curr_step);
};
#endif // ROTATE_SITE_H_
