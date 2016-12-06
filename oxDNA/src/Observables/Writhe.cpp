/*
 * Writhe.cpp
 *
 *  Created on: 16 Sep, 2015
 *      Author: Ferdinando Randisi
 */

#include "Writhe.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"
template<typename number>
Writhe<number>::Writhe() {

	_first_particle_index = 0;
	_last_particle_index = -1;

	_subdomain_size = -1;
	_N = -1;
	_use_default_go_around = true;
	_locate_plectonemes = false;
	_writhe_threshold = 0.28;
}

template<typename number>
Writhe<number>::~Writhe() {
	//deallocate them only if they have been allocated
	if (_N != -1){
		for (int i = 0; i < _N; i ++){
			delete []_writhe_integrand_values[i];
		}
		delete []_writhe_integrand_values;
	}

}

template<typename number>
void Writhe<number>::init(ConfigInfo<number> &config_info) {
	BaseObservable<number>::init(config_info);

	BaseParticle<number> **p = this->_config_info.particles;
	_N = *this->_config_info.N;

	//allocate the arrays - these take up unnecessary memory when the subdomain is smaller than N, but I don't think I care.
	_writhe_integrand_values = new number*[_N];
	for (int i = 0; i < _N; i ++){
		_writhe_integrand_values[i]= new number[ _N ];
	}
	for ( int i = 0; i < _N; i ++){
		for ( int j = 0; j < _N ;j++){
			_writhe_integrand_values[i][j] = -1e9;
		}
	}
	// check that _first_particle_index is in [0,N) and not the terminal particle
	if (_first_particle_index < 0 || _first_particle_index > _N-1){
		throw oxDNAException("Writhe: first_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",_N,_first_particle_index);	
	}
	if (p[_first_particle_index]->n5 == P_VIRTUAL){
		throw oxDNAException("Writhe: first_particle_index must not be the index of the last particle of a strand, otherwise which particle should be last_particle_index referring to?");
	}


	if (_last_particle_index == -1){
		_last_particle_index = _first_particle_index;
		do{ // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself 
			_last_particle_index = p[_last_particle_index]->n5->index;

		}	while (p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	if(_first_particle_index == _last_particle_index){
		_last_particle_index = p[_first_particle_index]->n3->index;
	}
	// check that _last_particle_index is in [0,N)
	if (_last_particle_index < 0 || _last_particle_index > _N-1){
		throw oxDNAException("Writhe: last_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",_N,_last_particle_index);	
	}
	// check that _first_particle_index is different than _last_particle_index
	if (_first_particle_index == _last_particle_index){
		throw oxDNAException("Writhe: first_particle_index can't be equal to last_particle_index. If you're talking about a circular molecule then just don't set last_particle_index.");
	}
	// if the molecule is not circular, the last bead is not considered since it just mirrors the behaviour of the first-but-last bead
	if (p[_last_particle_index]->n5 == P_VIRTUAL){
		_last_particle_index = p[_last_particle_index]->n3->index;
	}
	// check that first and last particle are on the same strand.
	if (p[_first_particle_index]->strand_id != p[_last_particle_index]->strand_id){
		throw oxDNAException("In observable writhe, the first particle (index %d) and the last particle (index %d) are not on the same strand. They're supposed to define a domain of topologically adjacent particles, but obviously they don't.",_first_particle_index,_last_particle_index);
	}
	// check that you can actually start from the first particle and get to the last one going forward along the chain.
	int test_index = p[_first_particle_index]->n5->index;
	int N_strand = 2;
	while(p[test_index]->n5 != P_VIRTUAL && test_index != _first_particle_index && test_index != _last_particle_index)
	{	
		test_index = p[test_index]->n5->index;
		N_strand++;
	} 
	if (p[test_index]->n5 == P_VIRTUAL){
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward.",_first_particle_index,_last_particle_index);
	}
	if (p[test_index]->n5->index == _first_particle_index && test_index != _last_particle_index){
		//throw oxDNAException("Bye!");
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward. This is very strange since they are both on the same strand as far as I know, so one of the developers (probably Ferdinando) messed something up. Please report the occurrence of this error to the developers.",_first_particle_index,_last_particle_index);
	}

	int default_subdomain_size = (_last_particle_index - _first_particle_index) -1;
	if ( _subdomain_size == -1){
		if (_locate_plectonemes) _subdomain_size = 35;
		else _subdomain_size = default_subdomain_size;
	}
	if ( _subdomain_size >= _last_particle_index - _first_particle_index){
		throw oxDNAException("In observable Writhe, subdomain_size %d should be strictly less than the difference between last_particle_index and first_particle_index.");
	}
	if (_use_default_go_around){
		// set the _go_around variable
		if (p[_last_particle_index]->n5->index == _first_particle_index && _subdomain_size != default_subdomain_size) _go_around = true;
		else _go_around = false;
	}

}

template<typename number>
void Writhe<number>::get_settings(input_file &my_inp, input_file &sim_inp) {

	getInputInt(&my_inp,"first_particle_index",&_first_particle_index,0);

	getInputInt(&my_inp,"last_particle_index",&_last_particle_index,0);

	getInputInt(&my_inp,"subdomain_size",&_subdomain_size,0);

	//if the user has set the value of go_around, keep track of it so that we don't overwrite it.
	if( getInputBool(&my_inp,"go_around",&_go_around,0) == KEY_FOUND){
		_use_default_go_around = false;
	}

	getInputBool(&my_inp,"locate_plectonemes",&_locate_plectonemes,0);
	getInputNumber(&my_inp,"writhe_threshold",&_writhe_threshold,0);

}

template<typename number>
std::string Writhe<number>::get_output_string(llint curr_step) {
	string result; 
	BaseParticle<number> **p = this->_config_info.particles;
	LR_vector<number> r, rp, t, tp;

	number writhe = 0;
	//number writhetemp = 0;
	// first compute the value of the writhe integrand
	for ( int i = _first_particle_index; i <= _last_particle_index; i++){
		BaseParticle<number> * i_n5_particle = p[i]->n5;

		t = (i_n5_particle->pos - p[i]->pos);
		t.normalize();
		r = p[i]->pos;

		// this loop starts from i+2 to disregard adjacent segments. 
		// if two segments are adjacent then r-r'=-t, so we have - t x t' . t = 0
		for ( int j = (int)(i-_subdomain_size -1); j < i; j++) {
			// if we should go around, then we have to compute the ''negative'' values of j as well.
			int jj;
			if (j < _first_particle_index){
				if (_go_around){
					jj = j + _last_particle_index +1 - _first_particle_index;
				}
				// otherwise we don't need to.
				else {
					continue;
				}
			}
			else{
				jj = j;
			}
			BaseParticle<number> * jj_n5_particle = p[jj]->n5;
			//the following condition should never be met - if it is, an edge case has not been treated properly.
			if (jj_n5_particle == P_VIRTUAL){
				throw oxDNAException("In observable writhe, Skipping i %d jj %d. This wasn't supposed to happen!",i,jj);
			}

			tp = (jj_n5_particle->pos - p[jj]->pos);
			tp.normalize();
			rp = p[jj]->pos;
			_writhe_integrand_values[i][jj] = _writheIntegrand(t,tp,r,rp);
		}
	} 
	//then perform the addition in order to compute the (local) writhe
	bool on_peak = false;
	double peak_value = -1;
	int  peak_position = -1;
	bool stop = false;
	for(int k = _first_particle_index; k <= _last_particle_index; k++){
		if ( k >= _last_particle_index - _subdomain_size && ! _go_around) break;
		writhe = 0;
		for(int i = k + 1; i < k + _subdomain_size; i++){
			for(int j = k; j < i; j++){
				int actual_i = i > _last_particle_index ? i - _last_particle_index + _first_particle_index + 1: i;
				int actual_j = j > _last_particle_index ? j - _last_particle_index + _first_particle_index + 1: j;
				if( _writhe_integrand_values[actual_i][actual_j] < -1e8){ 
					OX_LOG(Logger::LOG_INFO,"Obsevrable writhe: problem with %d %d %lf\n",actual_i,actual_j,_writhe_integrand_values[actual_i][actual_j]);
					stop = true;
				}

				writhe += _writhe_integrand_values[actual_i][actual_j];
			}
		}
		char temp[512]={0};
		if (_locate_plectonemes){
	  	// If the writhe is higher than the threshold, there's a peak.
	  	if (!on_peak){
	    	if (abs(writhe) > _writhe_threshold){
					
	      	on_peak = true;
	      	peak_value = abs(writhe);
	      	peak_position = k; 
				}
			}
	  	else{
	  		// When on a peak, look for the maximum.
	    	if (abs(writhe) > peak_value){
	      	peak_value = abs(writhe);
	      	peak_position = (int)(k+_subdomain_size/2.);
				}
	  		// When the writhe goes below the threshold, the peak is over.
	    	else if (abs(writhe) < _writhe_threshold){ 
	      	on_peak = false;
					if (peak_position > _last_particle_index){
						peak_position -= ( _last_particle_index - _first_particle_index + 1);
					}
					sprintf(temp,"%d\t",peak_position);
					result += std::string(temp);
				}
			}

		}
		else{
			sprintf(temp,"%14.14lf\n",writhe);
			result += std::string(temp);
		}
	}
	if (stop) throw oxDNAException("Dying badly because of problems with the observable writhe. That's not supposed to happen!");


	return result;
}


/* 
//I'm labelling segments, so the i index goes from 0 to the label of the first_but_last particle, N-2
// There are two for loops: the first computes the values where j goes from i+1 to i+_subdom_size-1
// The other computes the ones in which j goes from i+1 to N-2.
//	Example for N = 7, _subdomain_size = 4
//
//i	
//
//6|_|_|_|_|_|_|/|
//5|_|_|_|_|_|/|o|
//4|_|_|_|_|/|o|o|
//3|_|_|_|/|x|x|x|
//2|_|_|/|x|x|x|_|
//1|_|/|x|x|x|_|_|
//0|/|x|x|x|_|_|_|
//  0 1 2 3 4 5 6  j
//
//  The first for loop computes the x values, the second computes the o values.
//  Notice that since j goes from i+1 to something, j is labeled starting from 0 and ending to _subdomain_size-1
*/

template class Writhe<float>;
template class Writhe<double>;
