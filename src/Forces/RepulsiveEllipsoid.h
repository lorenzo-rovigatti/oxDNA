/**
 * @file    RepulsiveEllipsoid.h
 * @date    26/mag/2021
 * @author  Andrea
 *
 */

#ifndef REPULSIVEELLIPSOID_H_
#define REPULSIVEELLIPSOID_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into an ellipsoid.
 *
 * example section in the external forces file:

 \n
 {\n
 particle = -1    # acts on all particles\n
 type = ellipse \n
 r0 = 7.,7.,7.    # semi-axis    in simulation unit lengths\n
 stiff = 10.      # quite stiff to confine\n
 rate = 0.        # constant radius\n
 center = 0,0,0   # center position\n
 }\n\n

 * @verbatim
 stiff = <float> (stiffness of the repulsion.)
 [r0 = <float>,<float>,<float> (semi-axis of the ellipse, in simulation units.)]
 rate = <float> (rate of growth of the radius. Note that the growth is linear in timesteps/MC steps, not reduced time units.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles)
 [center = <float>,<float>,<float> (centre of the sphere, defaults to 0,0,0)]
 @endverbatim
 */
class RepulsiveEllipsoid: public BaseForce {
public:
	/// center of the ellipsoid
	LR_vector _centre;

	/// outer axes of the ellipsoid: for distances smaller than these values particles will not feel any force
	LR_vector _r_2;

	/// inner axes: for distances greater than these values particles will not feel any force
	LR_vector _r_1;

	RepulsiveEllipsoid();
	virtual ~RepulsiveEllipsoid() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // REPULSIVEELLIPSOID
