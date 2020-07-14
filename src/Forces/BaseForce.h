/**
 * @file    BaseForce.h
 * @date    18/oct/2011
 * @author  Flavio 
 *
 */

#ifndef BASEFORCE_H_
#define BASEFORCE_H_

#include <string>
#include "../defs.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/Utils.h"
#include "../Utilities/ConfigInfo.h"

// forward declarations of BaseParticle and BaseBox; needed to compile
class BaseParticle;
class BaseBox;

/**
 * @brief Base class for external forces. All external forces inherit from here.
 *
 * Note: This class contains public members with names starting with underscores.
 * We have to change scope policy due to the fact that GPU classes
 * require access to these members.
 *
 */

class BaseForce {
private:
	/**
	 * @brief A name of the group this force belongs to.
	 *
	 * Different forces can be grouped under the same name. This can be used to
	 * act separately on different forces. For example, it can be used together
	 * with the observable ForceEnergy to print only the energy due to specific
	 * groups of forces.
	 */
	std::string _group_name;

protected:
	/**
	 * @brief Adds the current force to the particle(s) listed in particle_string
	 *
	 * This method internally uses Utils::get_particles_from_string to extract a list of particles from the particle_string parameter
	 * and then use it to initialise all the particles contained therein.
	 *
	 * @param particles particle array
	 * @param particle_string a list of particles
	 * @param force_description an optional description (defaults to "generic force") that will be used in the logging messages
	 */
	void _add_self_to_particles(std::vector<BaseParticle *> &particles, std::string particle_string, std::string force_description = std::string("force"));

public:
	/**
	 * @brief standard members for forces
	 *
	 * we need these members to be public because
	 * we need access in order to copy these numbers
	 * to the GPU memory
	 */
	number _rate;
	number _F0;
	LR_vector _direction;
	LR_vector _pos0;
	number _stiff;
	int _site;
	BaseParticle * _p_ptr;

	BaseForce();
	virtual ~BaseForce();

	/** 
	 * @brief init function
	 *
	 * @return the list of particles it will act on and the force description
	 *
	 * This function initialises the force object and returns the list of particles it will act on
	 */
	virtual std::tuple<std::vector<int>, std::string> init(input_file &inp, BaseBox *box) = 0;

	virtual void set_group_name(std::string &name) {
		_group_name = name;
	}

	virtual std::string get_group_name() {
		return _group_name;
	}

	/**
	 * @brief returns value of the force (a vector)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual LR_vector value(llint step, LR_vector &pos) = 0;

	/**
	 * @brief returns value of the potential associated to the force (a number)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual number potential(llint step, LR_vector &pos) = 0;
};

using ForcePtr = std::shared_ptr<BaseForce>;

#endif /* BASEFORCE_H_ */
