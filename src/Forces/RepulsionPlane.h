/**
 * @file    RepulsionPlane.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIONPLANE_H_
#define REPULSIONPLANE_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a semispace.
 *
 * The definition of the plane is given by dir * (x,y,z) + position = 0.
 * Example section in the external forces file:

 \n
 {\n
 particle = -1    # acts on all particles\n
 type = repulsion_plane\n
 position = 7.    # in simulation unit lengths\n
 stiff = 50.      # quite stiff. Good for MC, not for MD\n
 dir = 0.,0.,1.   # points towards positive z\n
 }\n\n

 * @verbatim
 stiff = <float> (stiffness of the repulsion.)
 dir = <float>,<float>,<float> (the vector normal to the plane: it should point towards the half-plane where the repulsion is not acting.)
 position = <float> (defines the position of the plane along the direction identified by the plane normal.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles.)
 @endverbatim
 */

class RepulsionPlane: public BaseForce {
public:
	number _position;

	RepulsionPlane();
	virtual ~RepulsionPlane() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // MOVINGTRAP_H
