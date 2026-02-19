/**
 * @file    RepulsiveKeplerPoinsot.h
 *
 * Repulsive "Kepler–Poinsot style" star repulsor (no file imports).
 * Implemented as the UNION of 12 pentagonal pyramidal spikes aligned with
 * icosahedral directions (same symmetry family as the small stellated dodecahedron).
 *
 * This is an analytic, equation-defined repulsor suitable for oxDNA external forces.
 *
 * Input (external forces file):
 * {
 *   particle = -1
 *   type = repulsive_kepler_poinsot
 *   stiff = 10.0
 *   rate = 0.0
 *   center = 0,0,0
 *
 *   # Geometry at step 0 (all lengths in simulation units):
 *   apex = 1.20          # spike tip distance along each normal (from center)
 *   base = 0.70          # spike base plane distance along each normal
 *   base_radius = 0.45   # pentagon circumradius at the base plane
 *
 *   # Optional:
 *   kappa = 25.0         # smoothness for soft-max of pentagon support (bigger = sharper)
 * }
 */

#ifndef REPULSIVEKEPLERPOINSOT_H_
#define REPULSIVEKEPLERPOINSOT_H_

#include "BaseForce.h"

class RepulsiveKeplerPoinsot : public BaseForce {
public:
    /// Center of the star repulsor
    LR_vector _centre;

    /// Linear growth rate per step: growth = 1 + rate * step (applied to apex/base/base_radius)
    number _rate;

    /// Geometry at step 0
    number _apex;        // > _base
    number _base;        // >= 0
    number _base_radius; // > 0

    /// Soft-max sharpness for pentagon support
    number _kappa;

    RepulsiveKeplerPoinsot();
    virtual ~RepulsiveKeplerPoinsot() {}

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

    LR_vector value(llint step, LR_vector &pos) override;
    number potential(llint step, LR_vector &pos) override;
};

#endif // REPULSIVEKEPLERPOINSOT_H_