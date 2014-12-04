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

/**
 * @brief Base class for external forces. All external forces inherit from here.
 *
 * Note: This class contains public members with names starting with underscores.
 * We have to change scope policy due to the fact that GPU classes
 * require access to these members.
 *
 */

template <typename number> class BaseParticle;

template<typename number>
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
	LR_vector<number> _direction;
	LR_vector<number> _pos0;
	number _stiff;
	int _site;
	BaseParticle<number> * _p_ptr;

	BaseForce();
	virtual ~BaseForce();

	/**
	 * @brief get_settings function
	 *
	 * this function parses the input file
	 */
	virtual void get_settings(input_file &inp) = 0;

	/** 
	 * @brief init function
	 *
	 * This function initialises the force object and assignes 
	 * it to the relevant particles.
	 */
	virtual void init(BaseParticle<number> **particles, int N, number *box_side) = 0; 

	virtual void set_group_name(std::string &name) { _group_name = name; }
	virtual std::string get_group_name() { return _group_name; }

	/**
	 * @brief returns value of the force (a vector)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual LR_vector<number> value(llint step, LR_vector<number> &pos) = 0;
	
	/**
	 * @brief returns value of the potential associated to the force (a number)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual number potential (llint step, LR_vector<number> &pos) = 0;
};

template <typename number>
BaseForce<number>::BaseForce () {
	_F0 = -1.;
	_rate = -1.;
	_direction = LR_vector<number>(1., 0., 0.);
	_pos0 = LR_vector<number>(0., 0., 0.);
	_site = -1;
	_stiff = 0.;
	_p_ptr = P_VIRTUAL;
}

template <typename number>
BaseForce<number>::~BaseForce () {

}

#endif /* BASEFORCE_H_ */
