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
	std::string _group_name = "default";
	std::string _id = "";
	std::string _type = "type_unread";

public:
	/**
	 * @brief standard members for forces
	 *
	 * we need these members to be public because
	 * we need access in order to copy these numbers
	 * to the GPU memory
	 */
	number _rate = 0.;
	number _F0 = 0.;
	LR_vector _direction = LR_vector(1., 0., 0.);
	LR_vector _pos0 = LR_vector(0., 0., 0.);
	number _stiff;
	BaseParticle *_p_ptr;

	BaseForce();
	virtual ~BaseForce();

	/** 
	 * @brief init function
	 *
	 * @return the list of particles it will act on and the force description
	 *
	 * This function initialises the force object and returns the list of particles it will act on
	 */
	virtual std::tuple<std::vector<int>, std::string> init(input_file &inp);

	void set_group_name(std::string &name) {
		_group_name = name;
	}

	std::string get_group_name() {
		return _group_name;
	}

	void set_id(std::string id) {
		_id = id;
	}

	std::string get_id() {
		return _id;
	}

	std::string get_type() {
		return _type;
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
