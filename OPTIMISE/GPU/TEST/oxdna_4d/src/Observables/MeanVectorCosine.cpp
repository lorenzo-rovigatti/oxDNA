/*
 * MeanVectorCosine.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Ferdinando Randisi
 */

#include "MeanVectorCosine.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"

MeanVectorCosine::MeanVectorCosine() {
	_chain_id = -1;
}

MeanVectorCosine::~MeanVectorCosine() {

}

void MeanVectorCosine::init() {
	BaseObservable::init();

	// check that _vector_to_average is either 1,2, or 3.
	if(_vector_to_average > 3 || _vector_to_average < -1) {
		throw oxDNAException("MeanVectorCosine: vector_to_average must be 0,1,2 or 3, but was set to %d.", _vector_to_average);
	}
	// check that _first_particle_position is less than _last_particle_position
	if(_first_particle_position > _last_particle_position && _last_particle_position > -1) throw oxDNAException("MeanVectorCosine: first_particle_position must be equal or less than last_particle_position, but they are respectively %d and %d.", _first_particle_position, _last_particle_position);
	//set_first_last_particle_id(); //TODO: nel caso poi rimettere la chiamata a questa funzione
	//printf("first_last = %d %d\n",_first_particle_position,_last_particle_position);

}

void MeanVectorCosine::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	int tmp = 0;

	getInputInt(&my_inp, "chain_id", &tmp, 1);
	_chain_id = tmp;	//TODO: implement these checks also in the ParticlePosition observable.
	if(_chain_id < 0) throw oxDNAException("MeanVectorCosine, _chain_id should be >= 0 (_chain_id = %d).", _chain_id);
	//printf("_chain_id = %d\n",_chain_id);
	//if (_chain_id >= N_strands) throw oxDNAException("MeanVectorCosine, it should be _chain_id < N_strands (_chain_id = %d, N_strands = %d)",_chain_id,N_strands); //TODO: to implement this check, BaseObservable.h must be changed so that both N and N_strands are read into the config info.

	// if find_particle_position is not found, then set it to 0.
	if(getInputInt(&my_inp, "first_particle_position", &tmp, 0) == KEY_FOUND) {
		_first_particle_position = tmp;
		if(_first_particle_position < 0) throw oxDNAException("MeanVectorCosine, first_particle_position must be 0 or positive.");
	}
	else _first_particle_position = -1;
	// if find_particle_position is not found, then set it to N-2.
	if(getInputInt(&my_inp, "last_particle_position", &tmp, 0) == KEY_FOUND) {
		_last_particle_position = tmp;
		if(_last_particle_position < 0) throw oxDNAException("MeanVectorCosine, last_particle_position must be 0 or positive.");
	}
	else _last_particle_position = -1;

	// get the vector to average
	if(getInputInt(&my_inp, "vector_to_average", &tmp, 0) == KEY_FOUND) {
		_vector_to_average = tmp;
	}
	else _vector_to_average = 1;

}

std::string MeanVectorCosine::get_output_string(llint curr_step) {
	int sign;
	std::string result;
	double average = 0.;
	std::vector<BaseParticle *> &p = _config_info->particles();
	double delta_omega = 0;

	//FILE *fp = fopen("myv-1.txt","w");//TODO: comment this line when the problem is solved

	for(int i = 0; i < _config_info->N() - 2; i++) {
		if(_vector_to_average == -1) {
			_u = p[i]->orientationT.v1;
			_up = p[i]->n5->orientationT.v1;

			_f = p[i]->orientationT.v2;
			_fp = p[i]->n5->orientationT.v2;

			_v = p[i]->orientationT.v3;
			_vp = p[i]->n5->orientationT.v3;

			//sign = _up * (_f.cross(_fp)) > 0 ? 1 : -1;
			if(_up * (_f.cross(_fp)) > 0) sign = 1;
			else sign = -1;
			delta_omega = LRACOS((_f * _fp + _v * _vp) / (1 + _u * _up));

			average += delta_omega * sign;
			//fprintf(fp,"%d %d %.7lf %lf\n",i,sign,delta_omega,average);//TODO: comment this line when the problem is solved

		}
		else if(_vector_to_average == 1) {
			// u * u
			_u = p[i]->orientationT.v1;
			_up = p[i]->n5->orientationT.v1;

			average += _u * _up;
		}
		else if(_vector_to_average == 2) {
			// u * t
			_u = p[i]->orientationT.v1;
			LR_vector t = p[i]->n5->pos - p[i]->pos;
			average += _u * t / t.module();
		}
		else if(_vector_to_average == 3) {
			if(p[i]->n5->n5 != P_VIRTUAL) {
				// t * t
				LR_vector t = p[i]->n5->pos - p[i]->pos;
				LR_vector tp = p[i]->n5->n5->pos - p[i]->n5->pos;
				average += t * tp / (t.module() * tp.module());
			}

		}
		else {
			_u = p[i]->orientationT.v1;
			_up = p[i]->n5->orientationT.v1;

			_f = p[i]->orientationT.v2;
			_fp = p[i]->n5->orientationT.v2;

			_v = p[i]->orientationT.v3;
			_vp = p[i]->n5->orientationT.v3;

			average += (_f * _fp + _v * _vp) / (1 + _u * _up);
		}
	}

	if(_vector_to_average == -1) {
		if(average > M_PI || average < -M_PI) result = Utils::sformat("%14.4lf", -50.0);
		else result = Utils::sformat("%14.14lf", cos(average));
	}
	else if(_vector_to_average == 3) {
		result = Utils::sformat("%14.14lf", average / (_config_info->N() - 2));
	}
	else {
		result = Utils::sformat("%14.14lf", average / (_config_info->N() - 2));
	}
	//fclose(fp);abort();//TODO: comment this line when the problem is solved
	return result;
}

void MeanVectorCosine::set_first_last_particle_id() {
	int N = _config_info->N();
	std::vector<BaseParticle *> &p = _config_info->particles();
	int chain_counter = 0;
	int particle_counter = 0;
	int number_of_values = 0;
	//TODO: discuss these checks with Flavio and see if he agrees.
	char explanation[] = { "MeanVectorCosine: currently this observable only works if the order of strands is well defined in the topology file (e.g., the first particle is the first particle of a strand, the last particle is the last particle of a strand, and neighbouring particles are listed next to each other). Implementing the mean vector's cosine for other cases, in which the topology file is more messy, can easily be done, but design choices have to be taken to decide how this should be done." };	// An option could be to create a data structure which stores the chain-type, with particle index of first and last. This structure is probably necessary once we implement the mixed model.
	if(_config_info->particles()[0]->n3 != P_VIRTUAL) throw oxDNAException("%s\n --- The first particle in the particles array is not the first particle of a strand", explanation);
	if(_config_info->particles()[N - 1]->n5 != P_VIRTUAL) throw oxDNAException("%s\n --- The last particle in the particles array is not the last particle of a strand", explanation);

	_first_particle_id = -1;
	_last_particle_id = -1;
	for(int i = 0; i < N; i++) {
		if(p[i]->n3 != P_VIRTUAL && p[i]->n3->index != i - 1) throw oxDNAException("%s\n --- The particle %d is preceded by particle %d instead that by particle %d", i, p[i]->n3->index, i - 1);
		if(chain_counter == _chain_id) {	//if it gets to the chain of interest, find the index of the particles
			if(p[i]->n3 == P_VIRTUAL && _first_particle_position == -1) _first_particle_id = i;
			else if(particle_counter == _first_particle_position) _first_particle_id = i;
			if(p[i]->n5 == P_VIRTUAL && _last_particle_position == -1) _last_particle_id = i - 1;
			else if(particle_counter == _last_particle_position) _last_particle_id = i;

			//printf("_first_particle_id = %d _last_particle_id = %d\n",_first_particle_id,_last_particle_id);
			if(_first_particle_id != -1) {
				number_of_values++;
				printf("number_of_values %d\n", number_of_values);
			}
			// exit the loop if _last_particle_id is set or if the chain is over.
			if(_last_particle_id != -1 || chain_counter > _chain_id) break;
			particle_counter++;
		}
		if(p[i]->n5 == P_VIRTUAL) {
			chain_counter++;
			printf("Chain_counter %d\n", chain_counter);
		}	//when it finds a terminal particle, increment the chain counter.
	}
	// throw an exception if either the first or last particle were not found, or if the last particle is a terminal one.
	if(_first_particle_id == -1) throw oxDNAException("MeanVectorCosine: could not find particle %d of strand %d", _first_particle_position, _chain_id);
	if(_last_particle_id == -1) throw oxDNAException("MeanVectorCosine: could not find particle %d strand %d", _last_particle_position, _chain_id);
	if(p[_last_particle_id]->n5 == P_VIRTUAL) throw oxDNAException("MeanVectorCosine: particle specified by the argument last_particle_position must not be a terminal particle of the chain. If you want to set it to the end of the chain, just do not specify the argument last_particle_position.");
	_one_over_number_of_values = 1. / double(number_of_values);
	printf("number_of_values %d\n", number_of_values);
	return;
}
