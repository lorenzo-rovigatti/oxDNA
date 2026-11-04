/**
 * @file    SkewTrap.h
 * @date    03/nov/2026
 * @author  Erik
 *
 */

#ifndef SKEWTRAP_H_
#define SKEWTRAP_H_

#include "BaseForce.h"

/**
 * @brief Asymmetric harmonic-like trap between two particles.
 *
 * With x = r - r0 the potential is minus the log of a skew-normal density with scale s and shape a:
 *
 *   V(x) = x^2 / (2 s^2) + ln(sqrt(2 pi) s) - ln erfc(w),   w = -a x / (s sqrt(2))
 *
 * On the steep side (w > 0) erfc(w) is evaluated through the scaled function erfcx(w) = exp(w^2) erfc(w),
 * which keeps both the potential and the force accurate for any x.
 */
class SkewTrap: public BaseForce {
private:
	int _particle;
	int _ref_id;

public:
	// these are public so that the CUDA backend can copy them
	BaseParticle *_p_ptr;
	number _r0;
	number _a; // shape parameter
	number _s; // scale parameter ("stdev")
	bool PBC;

	// precomputed constants
	number _inv_s2; // 1 / s^2
	number _w_coeff; // -a / (s sqrt(2)), so that w = _w_coeff * x
	number _g_coeff; // a s sqrt(2 / pi)
	number _V_const; // ln(sqrt(2 pi) s)

	SkewTrap();
	virtual ~SkewTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector force(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

	/// potential as a function of x = r - r0
	number potential_of_x(number x) const;
	/// -dV/dx as a function of x = r - r0 (positive values push the particles apart)
	number force_mag_of_x(number x) const;

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // SKEWTRAP_H
