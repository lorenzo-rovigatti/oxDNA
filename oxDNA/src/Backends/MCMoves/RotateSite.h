/**
 * @file    RotateSite.h
 * @date    23/mar/2016
 * @author  flavio
 *
 */

#ifndef ROTATE_SITE_H_
#define ROTATE_SITE_H_

#include "BaseMove.h"

template<typename number>
class RotateSite : public BaseMove<number> {
	protected:
		number _delta;
		std::vector<int> allowed_sites;
		
	public:
		RotateSite(ConfigInfo<number> *Info);
		virtual ~RotateSite();

		virtual void get_settings(input_file &inp, input_file &sim_inp);
		void apply (llint curr_step);
};
#endif // ROTATE_SITE_H_
