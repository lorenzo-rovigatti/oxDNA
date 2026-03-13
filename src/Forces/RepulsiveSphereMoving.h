/**
 * @file    custom_force.h
 * @date    22/oct/2025
 * @author  Victor
 */

#ifndef REPULSIVESPHEREMOVING_H_
#define REPULSIVESPHEREMOVING_H_

#include <tuple>
#include <vector>
#include <string>

#include "BaseForce.h"

class RepulsiveSphereMoving : public BaseForce {
public:
    RepulsiveSphereMoving();
    ~RepulsiveSphereMoving() override = default;

    // oxDNA hooks
    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;
    LR_vector value(llint step, LR_vector &pos) override;
    number    potential(llint step, LR_vector &pos) override;
    
    // --- CUDA/host-side accessors (needed by CUDAForces.h) ---
    number    stiff()  const { return _stiff; }
    number    r0()     const { return _r0; }
    number    rate()   const { return _rate; }
    number    r_ext()  const { return _r_ext; }

    LR_vector origin() const { return _origin; }
    LR_vector target() const { return _target; }
    llint     steps()  const { return _steps; }


private:
    // Interpolated center at a given MD step
    LR_vector center_for_step(llint step) const;

    // Parameters
    number    _stiff{};          // harmonic stiffness
    number    _r0{};             // initial radius
    number    _rate{};           // growth rate (per step)
    number    _r_ext{1e10};      // outer cutoff

    // Positions
    LR_vector _center{0., 0., 0.};  // optional static center if used
    LR_vector _origin{0., 0., 0.};  // start of motion
    LR_vector _target{0., 0., 0.};  // end of motion
    llint     _steps{0};            // steps to go origin -> target
};

#endif // REPULSIVESPHEREMOVING_H_

