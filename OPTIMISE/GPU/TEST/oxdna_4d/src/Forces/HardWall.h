/**
 * @file    HardWall.h
 * @date    5/oct/2016
 * @author  Lorenzo
 *
 */

#ifndef HARDWALL_H_
#define HARDWALL_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a semispace. Given its "hard" nature, it may be used only iin Monte Carlo simulations.
 *
 * The definition of the plane is given by dir * (x,y,z) + position = 0.
 * Example section in the external forces file:

 \n
 {\n
 particle = -1          # acts on all particles\n
 type = hard_wall\n
 position = 7.          # in simulation unit lengths\n
 dir = 0.,0.,1.         # points towards positive z\n
 sigma = 1.5            # diameter of the wall
 }\n\n

 * @verbatim
 dir = <float>,<float>,<float> (the vector normal to the plane: it should point towards the half-plane that is allowed to the particles.)
 position = <float> (defines the position of the plane along the direction identified by the plane normal.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles.)
 [sigma = <float> ("Diameter" of the wall. It effectively rescales the distance between particle and wall. Defaults to 1.)]
 @endverbatim
 */
class HardWall: public BaseForce {
public:
	number _position;
	number _sigma;

	HardWall();
	virtual ~HardWall() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // HARDWALL_H
