/**
 * @file    RepulsiveSphere.h
 * @date    28/nov/2014
 * @author  Flavio 
 *
 */

#ifndef REPULSIVESPHERE_H_
#define REPULSIVESPHERE_H_

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
}\n\n

 * @verbatim
stiff = <float> (stiffness of the repulsion.)
r0 = <float> (radius of the sphere, in simulation units.)
rate = <float> (rate of growth of the radius. Note that the growth is linear in timesteps/MC steps, not reduced time units.)
@endverbatim
*/
template<typename number>
class RepulsiveSphere : public BaseForce<number> {
private:
	int _particle;

	/// center of the sphere
	LR_vector<number> _center;

	/// initial radius of the sphere and rate of growth (linear in timesteps/MC steps, not reduced time units)
	number _r0, _rate;

	/// pointer to the box side
	number * _box_side_ptr;

public:
	RepulsiveSphere();
	virtual ~RepulsiveSphere() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);
};

#endif // REPULSIVESPHERE
