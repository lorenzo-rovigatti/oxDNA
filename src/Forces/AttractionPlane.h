/**
 * @file    AttractionPlane.h
 * @date    01/aug/2024
 * @author  Matthies, Tilibit 
 *
 */

#ifndef ATTRACTIONPLANE_H_
#define ATTRACTIONPLANE_H_

#include "BaseForce.h"

/**
 * @brief External force field that attracts particles towards a plane.
 *
 * The definition of the plane is given by dir * (x,y,z) + position = 0.
 * Example section in the external forces file:

\n
{\n
  particle = -1    # acts on all particles\n
  type = attraction_plane\n
  position = 7.    # in simulation unit lengths\n
  stiff = 50.      # quite stiff. Good for MC, not for MD\n
  dir = 0.,0.,1.   # points towards positive z\n
}\n\n

 * @verbatim
stiff = <float> (stiffness of the attraction.)
dir = <float>,<float>,<float> (the vector normal to the plane: it should point towards the half-plane where the repulsion is not acting.)
position = <float> (defines the position of the plane along the direction identified by the plane normal.)
particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles.)
@endverbatim
*/

class AttractionPlane : public BaseForce {
public:
	number _position;

	AttractionPlane ();
	virtual ~AttractionPlane() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // ATTRACTIONPLANE_H
