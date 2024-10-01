/**
 * @file    RepulsiveSphereSmooth.h
 * @date    28/nov/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIVESPHERESMOOTH_H_
#define REPULSIVESPHERESMOOTH_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a sphere.
 *
 * example section in the external forces file:

 \n
 {\n
 particle = -1    # acts on all particles\n
 type = sphere \n
 r0 = 7.          # radius is 7 simulation unit lengths\n
 stiff = 50.      # quite stiff. Good for MC, not for MD\n
 rate = 0.        # constant radius\n
 center = 0,0,0
 }\n\n

 * @verbatim
 stiff = <float> (stiffness of the repulsion.)
 r0 = <float> (radius of the sphere, in simulation units.)
 rate = <float> (rate of growth of the radius. Note that the growth is linear in timesteps/MC steps, not reduced time units.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles)
 [center = <float>,<float>,<float> (centre of the sphere, defaults to 0,0,0)]
 @endverbatim
 */

class RepulsiveSphereSmooth: public BaseForce {
public:
	/// center of the sphere
	LR_vector _center;

	/// initial radius of the sphere 
	number _r0;
	/// external radius: for distances greater than r_ext, particles will not feel any force
	number _r_ext;
	/// smoothness of the force
	number _smooth;
	/// inflection point of the force
	number _alpha;

	RepulsiveSphereSmooth();
	virtual ~RepulsiveSphereSmooth() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // REPULSIVESPHERESMOOTH
