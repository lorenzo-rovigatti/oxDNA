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
}

template<typename number>
Writhe<number>::~Writhe() {

}

template<typename number>
void Writhe<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
	
	BaseParticle<number> **p = this->_config_info.particles;
	const int N = *this->_config_info.N;

	// check that _first_particle_index is in [0,N) and not the terminal particle
	if (_first_particle_index < 0 || _first_particle_index > N-1)
		throw oxDNAException("Writhe: first_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",N,_first_particle_index);	
	if (p[_first_particle_index]->n5 == P_VIRTUAL){
		throw oxDNAException("Writhe: first_particle_index must not be the index of the last particle of a strand, otherwise which particle should be last_particle_index referring to?");
	}


	if (_last_particle_index == -1){
		_last_particle_index = _first_particle_index;
		do{ // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself 
			_last_particle_index = p[_last_particle_index]->n5->index;

		}	while (p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	// check that _last_particle_index is in [0,N)
	if (_last_particle_index < 0 || _last_particle_index > N-1)
		throw oxDNAException("Writhe: last_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.",N,_last_particle_index);	
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
		throw oxDNAException("In observable writhe, the first particle (index %d) and the last particle (index %d) are not on the same strand. They're supposed to define a domain of topologically adjacent particles, but obviously they don't.");
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
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward.");
	}
	if (p[test_index]->n5->index == _first_particle_index){
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward. This is very strange since they are both on the same strand as far as I know, so one of the developers (probably Ferdinando) messed something up. Please report the occurrence of this error to the developers.");
	}
	
	if ( _subdomain_size == -1){
		_subdomain_size = N;
	}
	//declare the array that does the thing. TODO: this is actually bigger than it should be, but I can
	// initialise everything to a very negative value to make sure I don't go out of bounds.
	_writhe_integrand_values = new number*[N-1];
	for (int i = 0; i < N-1; i ++){
		_writhe_integrand_values[i]= new number[_subdomain_size - 1];
	}
	for ( int i = 0; i < N-1; i ++){
		for ( int j = 0; j < _subdomain_size - 1;j++){
			_writhe_integrand_values[i][j] = -1e9;
		}
	}
}

template<typename number>
void Writhe<number>::get_settings(input_file &my_inp, input_file &sim_inp) {

	getInputInt(&my_inp,"first_particle_index",&_first_particle_index,0);

	getInputInt(&my_inp,"last_particle_index",&_last_particle_index,0);

	// get the subdomain size.
	getInputInt(&my_inp,"subdomain_size",&_subdomain_size,0);

}

template<typename number>
std::string Writhe<number>::get_output_string(llint curr_step) {
	string result; 
	BaseParticle<number> **p = this->_config_info.particles;
	const int N = *this->_config_info.N;
	LR_vector<number> r, rp, t, tp;
	
	number writhe = 0;
	number writhetemp = 0;
	int i_p_index = _first_particle_index;
	for ( int i = 0; i < N - 2; i++){
		BaseParticle<number> * i_n5_particle = p[i_p_index]->n5;
		
		t = (i_n5_particle->pos - p[i_p_index]->pos);
		t.normalize();
		r = p[i_p_index]->pos;
		
		int j_p_index = i_p_index;
		// this loop starts from i+2 to disregard adjacent segments. 
		// if two segments are adjacent then r-r'=-t, so we have - t x t' . t = 0
		for ( int j = i + 2; j < N - 1; j++) {
			j_p_index = p[j_p_index]->n5->index;
			BaseParticle<number> * j_n5_particle = p[j_p_index]->n5;

			//printf("i=%dj=%d,ii=%d,jj=%d ",i,j,i_p_index,j_p_index);
			tp = (j_n5_particle->pos - p[j_p_index]->pos);
			tp.normalize();
			rp = p[j_p_index]->pos;
			//printf("%lf %lf %lf, %lf %lf %lf ",t.x,t.y,t.z,tp.x,tp.y,tp.z);
			//printf("%lf %lf %lf, %lf %lf %lf ",r.x,r.y,r.z,rp.x,rp.y,rp.z);
			
			writhetemp = _writheIntegrand(t,tp,r,rp);
			writhe += writhetemp;
			//printf("Writhe = %lf\t",writhe);
			//printf("Writhetemp = %14.14lf\n",writhetemp);

		}
		i_p_index = p[i_p_index]->n5->index;

	} 
			//printf("Writhe = %lf\n",writhe);
	return Utils::sformat("%14.14lf",writhe);

	
	/* I'm fairly sure that I'm doing something meaningful, but before doing this I should implement a simple writhe calculation that doesn't drive me insane and at the same time works. I should therefore be able to compute at least the global writhe with it. This is what the above stuff does.
	//TODO
	//summary of what's left to do:
	// populate the array with all the values of _writheIntegrand that are needed, and the others are neglected.
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
	int i_p_index = _first_particle_index;
	for ( int i = 0; i <= N - _subdomain_size; i++){
		BaseParticle<number> * i_n5_particle = p[i_p_index]->n5;

		t = (i_n5_particle->pos - p[i_p_index]->pos).normalize();
		r = p[i_p_index]->pos;
		
		int j_p_index = i_p_index;
		for ( int j = 0; j < _subdomain_size-1; j++){
			j_p_index = p[j_p_index]->n5->index;
			BaseParticle<number> * j_n5_particle = p[j_p_index]->n5;

			tp = (j_n5_particle->pos - p[j_p_index]->pos).normalize();
			rp = p[i_p_index]->pos;

			_writhe_integrand_values[i][j] = _writheIntegrand(t,tp,r,rp);
		}
		i_p_index = p[i_p_index]->n5->index;
	}
	for ( int i = N - _subdomain_size; i < N-2; i++){
		BaseParticle<number> * i_n5_particle = p[i_p_index]->n5;

		t = (i_n5_particle->pos - p[i_p_index]->pos).normalize();
		r = p[i_p_index]->pos;
		
		int j_p_index = i_p_index;
		// according to the table above, j should go from i + 1 to N - 2, so when we use an index that goes
		// from 0 to something it has to go until N - 1 - i
		for ( int j = 0; j < N - 1 - i ; j++){
			j_p_index = p[j_p_index]->n5->index;
			BaseParticle<number> * j_n5_particle = p[j_p_index]->n5;

			tp = (j_n5_particle->pos - p[j_p_index]->pos).normalize();
			rp = p[i_p_index]->pos;

			_writhe_integrand_values[i][j] = _writheIntegrand(t,tp,r,rp);
		}
		i_p_index = p[i_p_index]->n5->index;

	}
	
	// Now compute the integral on the first subdomain.
	number writhe = 0;
	int i_p_index = _first_particle_index;
	for ( int i = 0; i < _subdomain size - 1; i++){
		for (int j = 0; j < _subdomain_size - 1 - i; j++){
			writhe += _writhe_integrand_values[i][j];

		}
	}
		*/
	
	// if _subdomain_size==N_strand, just sum them all and call it a day.
	// figure out the first, and print it. If _subdomain_size == N_strand, then that's all we need, and we can
	// call it a day. Otherwise,
	// subtract the bottom row and add the rightmost column, and print it. That's the second.
	// do this until the end, every time adding to the string to output.
	//
	// All the following lines are kept here for inspiration.
	// ____________________________________________________________________________
/*

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
		*/
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
	
	return result;
		*/
}
/*
template<typename number>
std::string Writhe<number>::OLD_get_output_string(llint curr_step) {

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
		printf("Writhe\ti = %d p[i]->pos.x %lf u.x %lf u'.x %lf u*u' %lf \n",i,p[i]->pos.x, u.x, up.x, u * up);
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

template class Writhe<float>;
template class Writhe<double>;
