/**
 * @file    AttractionPlane.cpp
 * @date    01/aug/2024
 * @author  Matthies, Tilibit 
 *
 */

#include "AttractionPlane.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

AttractionPlane::AttractionPlane() :
				BaseForce() {
	_position = -1.;
}


std::tuple<std::vector<int>, std::string> AttractionPlane::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "position", &_position, 1);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	}
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "AttractionPlane");
	std::string description = Utils::sformat("AttractionPlane (stiff=%g, position=%g, dir=%g,%g,%g", _stiff, _position, _direction.x, _direction.y, _direction.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector AttractionPlane::value(llint step, LR_vector &pos) {
	number distance_from_plane = this->_direction*pos + this->_position;

	if(distance_from_plane >=  0.) {
		number force = -this->_stiff * 1.0;  // constant times * unit of length
		return force * this->_direction;
	}
	else
		return -(distance_from_plane*this->_stiff)*this->_direction;
}

number AttractionPlane::potential(llint step, LR_vector &pos) {
	number distance_from_plane = this->_direction*pos + this->_position;

	if(distance_from_plane >= 0.) 
		return  (number) (this->_stiff*distance_from_plane);
	else 
		return (number) (0.5*this->_stiff*SQR(distance_from_plane));
}

