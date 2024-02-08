/**
 * @file AttractionPlane.cpp
 * @brief Implementation file for the AttractionPlane class copied form Michael's file.
 * @author CopySubho
*/

#include "AttractionPlane.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

AttractionPlane::AttractionPlane() :
                BaseForce() {
    _particle = -1;
    _position = -1.f;
}

number AttractionPlane::potential(llint step, LR_vector &pos) {
    number distance = this->_direction*pos + this->_position; // distance from the plane
    if(distance >=0.) return (number) this->_stiff*distance;
    else return (number) (0.5*this->_stiff*SQR(distance));
}

LR_vector AttractionPlane::value(llint step, LR_vector &pos) {
    number distance = this->_direction*pos + this->_position; // distance from the plane
    if(distance >=0.) return this->_stiff*this->_direction;
    else return -1.f*this->_stiff*distance*this->_direction;
}

std::tuple<std::vector<int>, std::string> AttractionPlane::init(input_file &inp) {
    BaseForce::init(inp);
    getInputInt(&inp, "particle", &_particle, 1);
    getInputNumber(&inp, "position", &_position, 1);
    getInputNumber(&inp, "stiff", &_stiff, 1);

    int tmpi;
    double tmpf[3];
    std::string strdir;
    getInputString(&inp, "dir", strdir, 1);
    tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
    if(tmpi != 3)
        throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
    this->_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
    this->_direction.normalize();

    std::vector<BaseParticle *> &particles = CONFIG_INFO->particles();

    int N = particles.size();
    if(_particle >=N || N < -1) throw oxDNAException("Trying to add a AttractionPlane on non-existent particle %d. Aborting", _particle);
    if(_particle !=-1){
        OX_LOG(Logger::LOG_INFO, "Adding AttractionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on particle %d", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z, _particle);
        particles[_particle]->add_ext_force(this);
    }else{ //add force to all the particles 
        OX_LOG (Logger::LOG_INFO, "Adding AttractionPlane (stiff=%g, position=%g, dir=%g,%g,%g, on ALL particles", this->_stiff, this->_position, this->_direction.x, this->_direction.y, this->_direction.z);
		for(int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
    }
    std::string description = "Attraction Plane has been created successfully!!!!";
    return std::make_tuple(std::vector<int>{_particle}, description);
}