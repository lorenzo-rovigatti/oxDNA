/*
 * MaxBackboneForce.h
 *
 * The cap on the force of FENE backbone bonds set by max_backbone_force (used to relax stressed structures).
 *
 * A FENE bond V(x) = -(eps / 2) log(1 - x^2 / delta2), where x is the deviation of the bond length from its
 * equilibrium value, exerts a force f(x) = eps x / (delta2 - x^2) that diverges at |x| = sqrt(delta2). With a maximum
 * force fmax the potential is replaced, for |x| > xmax (where f(xmax) = fmax), by
 *
 *   V(x) = finf |x| + (fmax - finf) xmax log|x| + c
 *
 * whose force decreases from fmax at xmax to finf (the "far" force) at large separations, and c makes the potential
 * continuous. The far force is either given in the input (max_backbone_force_far, oxDNA and oxRNA) or chosen so that
 * c = 0 (ANNaMo). With max_backbone_rate the maximum force grows with the simulation step:
 * fmax(step) = max_backbone_force + max_backbone_rate * step.
 *
 * An interaction with FENE bonds owns one object per kind of bond (e.g. the DNA and the RNA bonds of oxNA), passes
 * the input to get_settings() and the FENE parameters to set_fene(), calls update() before computing the bonded
 * interactions of a step and uses energy() and force() where applies() is true. Interactions whose FENE parameters
 * depend on the particles (oxDNA3) can use the static functions with their own parameters.
 */

#ifndef MAXBACKBONEFORCE_H_
#define MAXBACKBONEFORCE_H_

#include "../defs.h"

struct input_file;

class MaxBackboneForce {
public:
	/// how the far force is chosen
	enum class FarForce {
		INPUT,     ///< max_backbone_force_far (default 0.04), with c making the potential continuous
		CONTINUOUS ///< such that the potential is continuous with c = 0
	};

	explicit MaxBackboneForce(FarForce far_force = FarForce::INPUT);

	/// reads max_backbone_force, max_backbone_force_far and max_backbone_rate
	void get_settings(input_file &inp);

	/// sets the parameters of the FENE potential V(x) = -(eps / 2) log(1 - x^2 / delta2) that the cap applies to
	void set_fene(number eps, number delta2);

	/// logs the settings (to be called once, after set_fene())
	void log_settings(const char *interaction_name) const;

	bool enabled() const {
		return _enabled;
	}

	/// with max_backbone_rate, sets the maximum force of the current step; returns true if it changed
	bool update();

	number fmax() const {
		return _fmax;
	}

	/// the deviation from the FENE minimum beyond which the force is capped
	number xmax() const {
		return _xmax;
	}

	/// the force at large separations (finf)
	number far_force() const {
		return _finf;
	}

	/// the coefficient of log|x| in the capped potential, (fmax - finf) xmax
	number log_coefficient() const {
		return (_fmax - _finf) * _xmax;
	}

	/// the constant c of the capped potential
	number offset() const;

	/// true if the force of a bond with deviation x is capped
	bool applies(number x) const {
		return _enabled && fabs(x) > _xmax;
	}

	/// the capped potential at deviation x (|x| > xmax)
	number energy(number x) const {
		return energy(x, _fmax, _finf, _xmax, _eps, _delta2);
	}

	/// -dV/dx of the capped potential at deviation x (|x| > xmax)
	number force(number x) const {
		return force(x, _fmax, _finf, _xmax);
	}

	/// the deviation at which the force of the FENE potential with parameters eps and delta2 is fmax
	static number xmax(number fmax, number eps, number delta2);
	/// the capped potential at deviation x for the given parameters
	static number energy(number x, number fmax, number finf, number xmax, number eps, number delta2);
	/// -dV/dx of the capped potential at deviation x for the given parameters
	static number force(number x, number fmax, number finf, number xmax);

private:
	FarForce _far_force_mode;
	bool _enabled = false;
	bool _fene_set = false;
	number _fmax_initial = 0.;
	number _rate = 0.;
	llint _step = -1;

	number _fmax = 0.;
	number _finf = 0.;
	number _xmax = 0.;
	number _eps = 0.;
	number _delta2 = 0.;

	/// computes xmax (and, for FarForce::CONTINUOUS, finf) from fmax and the FENE parameters
	void _compute();
};

#endif /* MAXBACKBONEFORCE_H_ */
