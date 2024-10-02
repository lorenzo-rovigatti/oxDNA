/**
 * @file    RepulsionPlaneMoving.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIONPLANEMOVING_H_
#define REPULSIONPLANEMOVING_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a semispace defined by a normal and the position of one or more reference particles.
 *
 * The definition of the plane is given by dir * (x,y,z) + position = 0.
 * Example section in the external forces file:

 \n
 {\n
 particle = -1                      # acts on all particles\n
 ref_particle = 10                  # the position of the plane will be defined with respect to particle 10's position
 type = repulsion_plane_moving\n
 stiff = 50.                        # quite stiff. Good for MC, not for MD\n
 dir = 0.,0.,1.                     # points towards positive z\n
 }\n\n

 * @verbatim
 stiff = <float> (stiffness of the repulsion.)
 dir = <float>,<float>,<float> (the vector normal to the plane: it should point towards the half-plane where the repulsion is not acting.)
 particle = <int> (index(es) of the particle(s) on which the force shall be applied. Can be a list of comma-separated indexes. If -1, the force will be exerted on all the particles.)
 ref_particle = <int> (index(es) of the particle(s) whose position along the normal will be used to define the repulsive plane(s). Can be a list of comma-separated indexes.)
 @endverbatim
 */

class RepulsionPlaneMoving: public BaseForce {
private:
	std::string _particles_string;
	std::string _ref_particles_string;

public:
	std::vector<BaseParticle *> ref_p_ptr;
	int low_idx, high_idx;

	RepulsionPlaneMoving();
	virtual ~RepulsionPlaneMoving() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif // REPULSIONPLANEMOVING_H_
