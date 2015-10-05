/*
 * Contacts.cpp
 *
 *  Created on: Aug 31, 2015
 *      Author: Ferdinando Randisi
 */

#include "Contacts.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"
template<typename number>
Contacts<number>::Contacts() {
	
	_first_particle_index = 0;
	_last_particle_index = -1;
	_neighbours_to_ignore = 2;

	_contact_distance = 1.5;
	_only_outermost_contacts = false;
}

template<typename number>
Contacts<number>::~Contacts() {

}

template<typename number>
void Contacts<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
	
	BaseParticle<number> **p = this->_config_info.particles;
	const int N = *this->_config_info.N;

	// check that _first_particle_index is in [0,N) and not the terminal particle
	if (_first_particle_index < 0 || _first_particle_index > N-1)
		throw oxDNAException("Contacts: first_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",N,_first_particle_index);	
	if (p[_first_particle_index]->n5 == P_VIRTUAL){
		throw oxDNAException("Contacts: first_particle_index must not be the index of the last particle of a strand, otherwise which particle should be last_particle_index referring to?");
	}


	if (_last_particle_index == -1){
		_last_particle_index = _first_particle_index;
		do{ // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself 
			_last_particle_index = p[_last_particle_index]->n5->index;

		}	while (p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	// check that _last_particle_index is in [0,N)
	if (_last_particle_index < 0 || _last_particle_index > N-1)
		throw oxDNAException("Contacts: last_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",N,_last_particle_index);	
	// check that _first_particle_index is different than _last_particle_index
	if (_first_particle_index == _last_particle_index){
		throw oxDNAException("Contacts: first_particle_index can't be equal to last_particle_index. If you're talking about a circular molecule then just don't set last_particle_index.");
	}

	// if the molecule is not circular, the last bead is not considered since it just mirrors the behaviour of the first-but-last bead
	if (p[_last_particle_index]->n5 == P_VIRTUAL){
		_last_particle_index = p[_last_particle_index]->n3->index;
	}
	// check that first and last particle are on the same strand.
	if (p[_first_particle_index]->strand_id != p[_last_particle_index]->strand_id){
		throw oxDNAException("In observable Contacts, the first particle (index %d) and the last particle (index %d) are not on the same strand. They're supposed to define a domain of topologically adjacent particles, but obviously they don't.");
	}
	// check that you can actually start from the first particle and get to the last one going forward along the chain.
	int test_index = p[_first_particle_index]->n5->index;
	while(p[test_index]->n5 != P_VIRTUAL && test_index != _first_particle_index && test_index != _last_particle_index)
	{	
		test_index = p[test_index]->n5->index;
	} 
	if (p[test_index]->n5 == P_VIRTUAL){
		throw oxDNAException("In observable Contacts, could not get from particle %d to particle %d by going forward.\n Currently the observable Contacts works only for particles on the same strand. It can easily be modified for particles on different strands/free, but then instead of looping around particles using a list we must use the indices directly.");
	
	}
	if (p[test_index]->n5->index == _first_particle_index){
		throw oxDNAException("In observable Contacts, could not get from particle %d to particle %d by going forward. This is very strange since they are both on the same strand as far as I know, so one of the developers (probably Ferdinando) messed something up. Please report the occurrence of this error to the developers.");
	}
	
	//check that neighbours_to_ignore is 0 or greater.
	if (_neighbours_to_ignore < 0){
		throw oxDNAException("In observable Contacts, neighbours_to_ignore has to be 0 or higher, but was set to %d",_neighbours_to_ignore);
	}
	
	//check that contact_distance is actually 0 or higher
	if (_contact_distance < 0){
		throw oxDNAException("In observable Contacts, contact_distance has to be 0 or higher, but was set to %g",_neighbours_to_ignore);
	}
}

template<typename number>
void Contacts<number>::get_settings(input_file &my_inp, input_file &sim_inp) {

	getInputInt(&my_inp,"first_particle_index",&_first_particle_index,0);

	getInputInt(&my_inp,"last_particle_index",&_last_particle_index,0);
	
	getInputInt(&my_inp,"neighbours_to_ignore",&_neighbours_to_ignore,0);
	getInputNumber(&my_inp,"contact_distance",&_contact_distance,0);
	getInputBool(&my_inp,"only_outermost_contacts",&_only_outermost_contacts,0);


}

template<typename number>
std::string Contacts<number>::get_output_string(llint curr_step) {
	

	BaseParticle<number> **p = this->_config_info.particles;
	string result; 
	
	int first_but_last_particle = p[_last_particle_index]->n3->index;
	int innermost_particle = _first_particle_index;
	
	// loop over all the particles, starting from the first
	int i_particle = _first_particle_index;
	while (i_particle != first_but_last_particle){
		int j_particle = _last_particle_index;
		// j_particle goes from the last particle to the chain to the end
		// of the domain in which to look contacts (i.e. innermost_particle).
		// The end of the domain is given either by nth_i_neighbour or by innermost_particle,
		// whichever appears before. If only the outermost contacts are reported, innermost
		// particle will be set to the outermost particle to make a contact.
		bool refresh_innermost_particle = false;
		int nth_i_neighbour = i_particle;
		for (int i = 0; i <= _neighbours_to_ignore; i++){
			if ( innermost_particle == nth_i_neighbour ){
				refresh_innermost_particle = true;
			}
			nth_i_neighbour = p[nth_i_neighbour]->n5->index;
		}
		if ( refresh_innermost_particle ){
			innermost_particle = nth_i_neighbour;
		}
		
		// if i_particle and j_particle are within contact distance, print their indices.
		while (j_particle >= innermost_particle){
			LR_vector<number> r=p[i_particle]->pos - p[j_particle]->pos;
			//printf("%d %d %lf\n",i_particle,j_particle,r.module());

			if( r.module() < _contact_distance ){
				result += Utils::sformat("%d %d\t",i_particle,j_particle);	
				// if only the outermost contacts should be reported, then skip to the next i-particle
				// and don't loop through the j_particles before it.
				if( _only_outermost_contacts){
					innermost_particle = p[j_particle]->n5->index;
					break;
				}
			}
			
		
			j_particle = p[j_particle]->n3->index;
		}
	
		i_particle =p[i_particle]->n5->index; 
	}
	/*
	// ____________________________ Previous function
	int sign, number_of_values=0;
	string result; 
	BaseParticle<number> **p = this->_config_info.particles;
	double delta_omega=0;
	number result_number=0;
	LR_vector<number> u, up, v, vp, f, fp;


	int i = _first_particle_index;
	do{
		if ( _angle_index == 0){
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
			delta_omega = LRACOS(( f*fp + v*vp )/(1 + u*up))/(2*PI); 
			if(_print_local_details){
				result += Utils::sformat("%14.14lf\t",delta_omega * sign)+' ';
			}
			else{
				result_number += delta_omega * sign;
			}

		}
		else if ( _angle_index == 1){
			// u * u
			u =	p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			if(_print_local_details){
				result += Utils::sformat("%14.14lf\t",u*up)+' ';
			}
			else{
				result_number += u*up;
			}
		}
		else if ( _angle_index == 2){
			// f * f
			f =	p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;

			if(_print_local_details){
				result += Utils::sformat("%14.14lf\t",f*fp)+' ';
			}
			else{
				result_number += f*fp;
			}
		}
		else if ( _angle_index == 3){
			// v * v
			v =	p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;

			if(_print_local_details){
				result += Utils::sformat("%14.14lf\t",v*vp)+' ';
			}
			else{
				result_number += v*vp;
			}
		}
		// some of the things that MeanVectorCosine chose to output instead.
		// kept here just in case they are needed.
		*/
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
		/*
		i = p[i]->n5->index;
		number_of_values++;
	}while (i!= _last_particle_index);
	if (! _print_local_details){
		if (_angle_index == 0){
			result = Utils::sformat("%14.14lf",result_number);
		}
		else{
			result = Utils::sformat("%14.14lf",result_number/number_of_values);
		}
	}
	*/
	
	return result;
}
/*
template<typename number>
std::string Contacts<number>::OLD_get_output_string(llint curr_step) {

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
		printf("Contacts\ti = %d p[i]->pos.x %lf u.x %lf u'.x %lf u*u' %lf \n",i,p[i]->pos.x, u.x, up.x, u * up);
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

template class Contacts<float>;
template class Contacts<double>;
