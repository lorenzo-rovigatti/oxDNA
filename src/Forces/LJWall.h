/**
 * @file    LJWall.h
 * @date    04/dec/2015
 * @author  Lorenzo
 *
 */

#ifndef LJWALL_H_
#define LJWALL_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a semispace. It can feature a short-range attraction.
 *
 * This external force is similar in spirit to RepulsionPlane but it is much more tunable. The definition of the plane is given by dir * (x,y,z) + position = 0.
 * Example section in the external forces file:

 \n
 {\n
 particle = -1          # acts on all particles\n
 type = LJ_plane\n
 position = 7.          # in simulation unit lengths\n
 stiff = 50.            # quite stiff. Good for MC, not for MD\n
 dir = 0.,0.,1.         # points towards positive z\n
 sigma = 1.5            # diameter of the wall
 n = 10                 # the generalised Lennard-Jones exponent, its default value is 6
 generate_inside = true # useful when generating the initial configuration
 }\n\n

 * @verbatim
 dir = <float>,<float>,<float> (the vector normal to the plane: it should point towards the half-plane where the repulsion is not acting.)
 position = <float> (defines the position of the plane along the direction identified by the plane normal.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles.)
 [stiff = <float> (stiffness of the repulsion. Defaults to 1.)]
 [sigma = <float> ("Diameter" of the wall. It effectively rescales the distance between particle and wall. Defaults to 1.)]
 [n = <int> (Exponent of the 2n-n generalized Lennard-Jones expression. Defaults to 6.)]
 [only_repulsive = <bool> (If true, the interactio between particle and wall gets cut-off at the minimum, resulting in a purely-repulsive potential. Defaults to false.)]
 [generate_inside = <bool> (If true the wall-particle interaction may not diverge, even for negative distances. Useful when generating the starting configuration. Defaults to false)]
 @endverbatim
 */
class LJWall: public BaseForce {
private:
	bool _only_repulsive;
	bool _generate_inside;

public:
	int _n;
	number _position;
	number _sigma;
	number _cutoff;

	LJWall();
	virtual ~LJWall() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // LJWALL_H
