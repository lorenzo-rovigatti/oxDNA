/*
 * StressAutocorrelation.h
 *
 *  Created on: 09/jul/2022
 *      Author: lorenzo
 */

#ifndef STRESSAUTOCORRELATION_H_
#define STRESSAUTOCORRELATION_H_

#include "BaseObservable.h"

/**
 * @brief Outputs the instantaneous pressure, computed through the virial.
 *
 * The supported syntax is
 * @verbatim
 type = pressure (an observable that computes the osmotic pressure of the system)
 [stress_tensor = <bool> (if true, the output will contain 9 fields: the total pressure and the nine components of the stress tensor, xx, xy, xz, yx, yy, yz, zx, zy, zz)]
 @endverbatim
 */

class StressAutocorrelation: public BaseObservable {
	struct Level {
		uint m;
		uint p;
		const uint level_number;
		uint start_at;

		std::vector<double> data;
		std::vector<double> correlation;
		std::vector<uint> counter;
		double accumulator = 0.;
		uint accumulator_counter = 0;
		std::shared_ptr<Level> next = nullptr;

		Level(uint nm, uint np, uint l_number) :
			m(nm),
			p(np),
			level_number(l_number),
			data(p),
			correlation(p),
			counter(p) {
			start_at = (l_number == 0) ? 0 : p / m;
		}

		void add_value(double v) {
			data.insert(data.begin(), v);
			if(data.size() > p) {
				data.erase(data.begin() + p);
			}

			for(uint i = start_at; i < data.size(); i++) {
				correlation[i] += data[0] * data[i];
				counter[i]++;
			}

			accumulator += v;
			accumulator_counter++;

			if(accumulator_counter == m) {
				if(next == nullptr) {
					next = std::make_shared<Level>(m, p, level_number + 1);
				}

				next->add_value(accumulator / accumulator_counter);
				accumulator = 0.;
				accumulator_counter = 0;
			}
		}

		void get_times(double dt, std::vector<double>& times) {
			for(uint i = start_at; i < p; i++) {
				times.push_back(i * std::pow((double) m, (double) level_number) * dt);
			}

			if(next != nullptr) {
				next->get_times(dt, times);
			}
		}

		void get_acf(double dt, std::vector<double> &acf) {
			for(uint i = start_at; i < p; i++) {
				acf.push_back(correlation[i] / counter[i]);
			}

			if(next != nullptr) {
				next->get_acf(dt, acf);
			}
		}
	};

protected:
	std::vector<LR_vector> _old_forces, _old_torques;
	std::shared_ptr<Level> _sigma_xy, _sigma_yz, _sigma_zx, _N_xy, _N_yz, _N_xz;
	double _delta_t;

public:
	StressAutocorrelation();
	virtual ~StressAutocorrelation();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	void update_data(llint curr_step) override;

	std::string get_output_string(llint curr_step);
};

#endif /* STRESSAUTOCORRELATION_H_ */
