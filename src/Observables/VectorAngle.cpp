/*
 * VectorAngle.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Ferdinando Randisi
 */

#include "VectorAngle.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"
template<typename number>
VectorAngle<number>::VectorAngle() {
	
	_angle_index = 1;
	_first_particle_index = 0;
	_last_particle_index = -1;
}

template<typename number>
VectorAngle<number>::~VectorAngle() {

}

template<typename number>
void VectorAngle<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
	
	BaseParticle<number> **p = this->_config_info.particles;

	// check that _first_particle_index is less than _last_particle_index
	if (_first_particle_index > _last_particle_index && _last_particle_index > -1)
		throw oxDNAException("VectorAngle: first_particle_position must be equal or less than last_particle_position, but they are respectively %d and %d.",_first_particle_index,_last_particle_index);	

	if (_last_particle_index == -1){
		_last_particle_index = _first_particle_index;
		do{ // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself 
			_last_particle_index = p[_last_particle_index]->n5->index;

		}	while (p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	printf("first_particle = %d, last_particle = %d\n",_first_particle_index,_last_particle_index);
}

template<typename number>
void VectorAngle<number>::get_settings(input_file &my_inp, input_file &sim_inp) {

	// if find_particle_position is not found, then set it to 0.
	if (getInputInt(&my_inp,"first_particle_index",&_first_particle_index,0) == KEY_FOUND) {
		if (_first_particle_index < 0)
			throw oxDNAException("(VectorAngle.cpp) first_particle_index must be 0 or positive, but it was %d.",_first_particle_index);
	}
	// if find_particle_position is not found, then set it to N-2.
	if (getInputInt(&my_inp,"last_particle_index",&_last_particle_index,0) == KEY_FOUND) {
		if (_last_particle_index < 1)
			throw oxDNAException("(VectorAngle.cpp) last_particle_index must be 0 or positive, but it was %d.",_last_particle_index);
	}
	
	// get the angle to consider 
	if (getInputInt(&my_inp,"angle_index",&_angle_index,0) == KEY_FOUND){
		if (_angle_index < 0 || _angle_index > 3)
			throw oxDNAException("(VectorAngle.cpp) angle_index must be 0,1,2, or 3, but was set to %d.",_angle_index);
	}

}

template<typename number>
std::string VectorAngle<number>::get_output_string(llint curr_step) {
	int sign;
	string result; 
	BaseParticle<number> **p = this->_config_info.particles;
	double delta_omega=0;
	LR_vector<number> u, up, v, vp, f, fp;


	printf("sciacalla\n");
	int i = _first_particle_index;
	do{
		if ( _angle_index == -1){
			// omega + gamma = twisting angle
			u =	p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			f = p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;
			
			v = p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;
			
			if (up * (f.cross(fp)) > 0)
				sign = 1;
			else
				sign = -1;
			delta_omega = LRACOS(( f*fp + v*vp )/(1 + u*up)); 
			
			result += Utils::sformat("%14.14lf\t",delta_omega * sign)+' ';

		}
		else if ( _angle_index == 1){
			// u * u
			u =	p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			result += Utils::sformat("%14.14lf\t",u*up)+' ';
		}
		else if ( _angle_index == 2){
			// f * f
			f =	p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;

			result += Utils::sformat("%14.14lf\t",f*fp)+' ';
		}
		else if ( _angle_index == 3){
			// v * v
			v =	p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;

			result += Utils::sformat("%14.14lf\t",v*vp)+' ';
		}
		// some of the things that MeanVectorCosine chose to output instead.
		// kept here just in case they are needed.
		/*
		else if ( _angle_index == 2){
			// u * t
			u = p[i]->orientationT.v1;
			LR_vector<number> t = p[i]->n5->pos - p[i]->pos;
			average += u*t/t.module();
		}else if ( _angle_index == 3){
			if( p[i]->n5->n5 != P_VIRTUAL){
				// t * t
				LR_vector<number> t = p[i]->n5->pos - p[i]->pos;
				LR_vector<number> tp = p[i]->n5->n5->pos - p[i]->n5->pos;
				average += t*tp/(t.module()*tp.module());
			}
			

		} else {
			u =	p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			f = p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;
			
			v = p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;

			average +=( f*fp + v*vp )/(1 + u*up);
		}
		*/
	i = p[i]->n5->index;
	}while (i!= _last_particle_index);
	printf("fuori dal loop\n");
	return result;
}
/*
template<typename number>
std::string VectorAngle<number>::OLD_get_output_string(llint curr_step) {

	string result; number average = 0.;
	BaseParticle<number> **p = this->_config_info.particles;
	for( int i = _first_particle_id; i <= _last_particle_id; i++){
		if ( _vector_to_average == 1){
			_u =	p[i]->orientationT.v1;
			_up = p[i]->n5->orientationT.v1;

			average += _u*_up;
		} else if (_vector_to_average == 2){
			_f =	p[i]->orientationT.v2;
			_fp = p[i]->n5->orientationT.v2;

			average += _f*_fp;
		} else if (_vector_to_average == 3){
			_v =	p[i]->orientationT.v3;
			_vp = p[i]->n5->orientationT.v3;

			average += _v*_vp;
		} else if (_vector_to_average == 0){
			_u =	p[i]->orientationT.v1;
			_up = p[i]->n5->orientationT.v1;

			_f = p[i]->orientationT.v2;
			_fp = p[i]->n5->orientationT.v2;
			
			_v = p[i]->orientationT.v3;
			_vp = p[i]->n5->orientationT.v3;

			average +=( _f*_fp + _v*_vp )/(1 + _u*_up);
		}
//TODO: remove the -1 thing after we are sure it's equivalent to the 0 one
		 else if (_vector_to_average == -1){
			_u =	p[i]->orientationT.v2;
			_up = p[i]->n5->orientationT.v2;
			_v =	p[i]->orientationT.v3;
			_vp = p[i]->n5->orientationT.v3;
			average += (_u*_up)*(_v*_vp) - _u.cross(_up).module()*_v.cross(_vp).module();
		}
*/ /*
		printf("VectorAngle\ti = %d p[i]->pos.x %lf u.x %lf u'.x %lf u*u' %lf \n",i,p[i]->pos.x, u.x, up.x, u * up);
		printf("u.x %lf u.y %lf u.z %lf\n",u.x,u.y,u.z);
		printf("up.x %lf up.y %lf up.z %lf\n",up.x,up.y,up.z);
		printf("thex %lf they %lf thez %lf total %lf\n",u.x*up.x,u.y*up.y,u.z*up.z,u.x*up.x+u.y*up.y+u.z*up.z);
*/ /*
//		printf("average %lf\n",average);
	}
//		printf("average %lf\n",average*_one_over_number_of_values);
	result  =  Utils::sformat("%14.4lf", average*_one_over_number_of_values);
	return result;
}
*/

template class VectorAngle<float>;
template class VectorAngle<double>;
