/*
 * Writhe.cpp
 *
 *  Created on: 16 Sep, 2015
 *      Author: Ferdinando Randisi
 */

#include "Writhe.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"

Writhe::Writhe() {

	_first_particle_index = 0;
	_last_particle_index = -1;

	_subdomain_size = -1;
	_use_default_go_around = true;
	_locate_plectonemes = false;
	_writhe_threshold = 0.28;
	_print_space_pos = false;
	_print_size = false;
	_print_left_right = false;

	_contact_threshold = 5.;
	_size_outer_threshold = 30;
	_minimum_plectoneme_size = 1;
	_bending_angle_number_segments = 0;
	_particles_are_bases = true;

	_writhe_integrand_values = nullptr;
}

Writhe::~Writhe() {
	//deallocate them only if they have been allocated
	if(_writhe_integrand_values != nullptr) {
		for(int i = 0; i < _config_info->N(); i++) {
			delete[] _writhe_integrand_values[i];
		}
		delete[] _writhe_integrand_values;
	}

}

void Writhe::init() {
	BaseObservable::init();

	std::vector<BaseParticle *> &p = _config_info->particles();
	int N = _config_info->N();

	//allocate the arrays - these take up unnecessary memory when the subdomain is smaller than N, but I don't think I care.
	_writhe_integrand_values = new number*[N];
	for(int i = 0; i < N; i++) {
		_writhe_integrand_values[i] = new number[N];
	}
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			_writhe_integrand_values[i][j] = -1e9;
		}
	}
	// check that _first_particle_index is in [0,N) and not the terminal particle
	if(_first_particle_index < 0 || _first_particle_index > N - 1) {
		throw oxDNAException("Writhe: first_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.", N, _first_particle_index);
	}
	if(p[_first_particle_index]->n5 == P_VIRTUAL) {
		throw oxDNAException("Writhe: first_particle_index must not be the index of the last particle of a strand, otherwise which particle should be last_particle_index referring to?");
	}

	if(_last_particle_index == -1) {
		_last_particle_index = _first_particle_index;
		do { // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself
			_last_particle_index = p[_last_particle_index]->n5->index;

		} while(p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	if(_first_particle_index == _last_particle_index) {
		_last_particle_index = p[_first_particle_index]->n3->index;
	}
	// check that _last_particle_index is in [0,N)
	if(_last_particle_index < 0 || _last_particle_index > N - 1) {
		throw oxDNAException("Writhe: last_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.", N, _last_particle_index);
	}
	// check that _first_particle_index is different than _last_particle_index
	if(_first_particle_index == _last_particle_index) {
		throw oxDNAException("Writhe: first_particle_index can't be equal to last_particle_index. If you're talking about a circular molecule then just don't set last_particle_index.");
	}
	// if the molecule is not circular, the last bead is not considered since it just mirrors the behaviour of the first-but-last bead
	if(p[_last_particle_index]->n5 == P_VIRTUAL) {
		_last_particle_index = p[_last_particle_index]->n3->index;
	}
	//else if(_print_size or _print_left_right) throw oxDNAException("In observable writhe, print_size = true and print_left_right = true are incompatible with circular strands.");
	// check that first and last particle are on the same strand.
	if(p[_first_particle_index]->strand_id != p[_last_particle_index]->strand_id) {
		throw oxDNAException("In observable writhe, the first particle (index %d) and the last particle (index %d) are not on the same strand. They're supposed to define a domain of topologically adjacent particles, but obviously they don't.", _first_particle_index, _last_particle_index);
	}
	// check that you can actually start from the first particle and get to the last one going forward along the chain.
	int test_index = p[_first_particle_index]->n5->index;
	int N_strand = 2;
	while(p[test_index]->n5 != P_VIRTUAL && test_index != _first_particle_index && test_index != _last_particle_index) {
		test_index = p[test_index]->n5->index;
		N_strand++;
	}
	if(p[test_index]->n5 == P_VIRTUAL) {
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward.", _first_particle_index, _last_particle_index);
	}
	if(p[test_index]->n5->index == _first_particle_index && test_index != _last_particle_index) {
		//throw oxDNAException("Bye!");
		throw oxDNAException("In observable writhe, could not get from particle %d to particle %d by going forward. This is very strange since they are both on the same strand as far as I know, so one of the developers (probably Ferdinando) messed something up. Please report the occurrence of this error to the developers.", _first_particle_index, _last_particle_index);
	}

	int default_subdomain_size = (_last_particle_index - _first_particle_index) - 1;
	if(_subdomain_size == -1) {
		if(_locate_plectonemes) _subdomain_size = 35;
		else _subdomain_size = default_subdomain_size;
	}
	if(_subdomain_size >= _last_particle_index - _first_particle_index) {
		throw oxDNAException("In observable Writhe, subdomain_size %d should be strictly less than the difference between last_particle_index %d and first_particle_index %d.", _subdomain_size, _last_particle_index, _first_particle_index);
	}
	if(_use_default_go_around) {
		// set the _go_around variable
		if(p[_last_particle_index]->n5->index == _first_particle_index && _subdomain_size != default_subdomain_size) _go_around = true;
		else _go_around = false;
	}

}

void Writhe::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "first_particle_index", &_first_particle_index, 0);
	getInputInt(&my_inp, "last_particle_index", &_last_particle_index, 0);
	getInputInt(&my_inp, "subdomain_size", &_subdomain_size, 0);

	getInputBool(&my_inp, "print_space_position", &_print_space_pos, 0);
	getInputBool(&my_inp, "print_size", &_print_size, 0);
	getInputBool(&my_inp, "print_left_right", &_print_left_right, 0);
	if(_print_size and _print_left_right) throw oxDNAException("Writhe.cpp: print_size and print_left_right can't be both true.");
	getInputNumber(&my_inp, "contact_threshold", &_contact_threshold, 0);
	getInputInt(&my_inp, "size_outer_threshold", &_size_outer_threshold, 0);
	getInputInt(&my_inp, "minimum_plectoneme_size", &_minimum_plectoneme_size, 0);
	if(_minimum_plectoneme_size < 1) throw oxDNAException("In observable writhe, minimum_plectoneme_size must be 1 or greater.");
	getInputInt(&my_inp, "bending_angle_number_segments", &_bending_angle_number_segments, 0);

	//if the user has set the value of go_around, keep track of it so that we don't overwrite it.
	if(getInputBool(&my_inp, "go_around", &_go_around, 0) == KEY_FOUND) {
		_use_default_go_around = false;
	}

	getInputBool(&my_inp, "locate_plectonemes", &_locate_plectonemes, 0);
	getInputNumber(&my_inp, "writhe_threshold", &_writhe_threshold, 0);
	if(_writhe_threshold < 0) throw oxDNAException("In observable writhe, writhe_threshold is used to compare the absolute value of the local writhe, so it can't be set to a negative value.");

	// check whether we're using the oxDNA/oxRNA models, or kTEP
	std::string inter_type("DNA");
	if(getInputString(&sim_inp, "interaction_type", inter_type, 0) == KEY_FOUND) {
		if(inter_type.substr(0, 3) == "DNA" or inter_type.substr(0, 3) == "RNA") {
			_particles_are_bases = true;
		}
		else _particles_are_bases = false;
	}	// the default interaction is DNA
	else _particles_are_bases = true;
	if(_particles_are_bases) throw oxDNAException("Writhe observable NOT IMPLEMENTED for DNA or RNA.");

}

std::string Writhe::get_output_string(llint curr_step) {
	std::string result;
	std::vector<BaseParticle *> &p = _config_info->particles();
	int time = _config_info->curr_step;
	LR_vector r, rp, t, tp;

	number writhe = 0;
	//number writhetemp = 0;
	// first compute the value of the writhe integrand
	for(int i = _first_particle_index; i <= _last_particle_index; i++) {
		BaseParticle * i_n5_particle = p[i]->n5;

		t = (i_n5_particle->pos - p[i]->pos);
		t.normalize();
		r = p[i]->pos;

		// this loop starts from i+2 to disregard adjacent segments. 
		// if two segments are adjacent then r-r'=-t, so we have - t x t' . t = 0
		for(int j = (int) (i - _subdomain_size - 1); j < i; j++) {
			// if we should go around, then we have to compute the ''negative'' values of j as well.
			int jj;
			if(j < _first_particle_index) {
				if(_go_around) {
					jj = j + _last_particle_index + 1 - _first_particle_index;
				}
				// otherwise we don't need to.
				else {
					continue;
				}
			}
			else {
				jj = j;
			}
			BaseParticle * jj_n5_particle = p[jj]->n5;
			//the following condition should never be met - if it is, an edge case has not been treated properly.
			if(jj_n5_particle == P_VIRTUAL) {
				throw oxDNAException("In observable writhe, Skipping i %d jj %d. This wasn't supposed to happen!", i, jj);
			}

			tp = (jj_n5_particle->pos - p[jj]->pos);
			tp.normalize();
			rp = p[jj]->pos;
			_writhe_integrand_values[i][jj] = _writheIntegrand(t, tp, r, rp);
		}
	}
	//then perform the addition in order to compute the (local) writhe
	bool on_peak = false;
	double peak_value = -1;
	int peak_position = -1;
	bool stop = false;
	int final_particle_index = 0;
	if(_go_around) final_particle_index = _last_particle_index;
	else final_particle_index = _last_particle_index - _subdomain_size - 1;
	for(int k = _first_particle_index; k <= final_particle_index; k++) {
		/*
		 if(k >= _last_particle_index - _subdomain_size && !_go_around){
		 printf("break: k = %d\n",k);
		 break;
		 }
		 */
		writhe = 0;
		for(int i = k + 1; i < k + _subdomain_size; i++) {
			for(int j = k; j < i; j++) {
				int actual_i = i > _last_particle_index ? i - _last_particle_index + _first_particle_index + 1 : i;
				int actual_j = j > _last_particle_index ? j - _last_particle_index + _first_particle_index + 1 : j;
				if(_writhe_integrand_values[actual_i][actual_j] < -1e8) {
					OX_LOG(Logger::LOG_INFO,"Obsevrable writhe: problem with %d %d %lf\n",actual_i,actual_j,_writhe_integrand_values[actual_i][actual_j]);
					stop = true;
				}

				writhe += _writhe_integrand_values[actual_i][actual_j];
			}
		}
		char temp[512] = { 0 };
		if(_locate_plectonemes) {
			// If the writhe is higher than the threshold, there's a peak.
			if(!on_peak) {
				// don't get on a peak on the last particle index, because you will already have considered that at the beginning
				if(fabs(writhe) > _writhe_threshold and k != final_particle_index) {
					on_peak = true;
					peak_value = fabs(writhe);
					peak_position = k;
					//printf("Getting on peak at position %d\n",k);
				}
			}
			else {
				// When on a peak, look for the maximum.
				if(fabs(writhe) > peak_value) {
					peak_value = fabs(writhe);
					peak_position = (int) (k + _subdomain_size / 2.);
				}
				// When the writhe goes below the threshold, or when we have reached the end of the chain, the peak is over.
				if(fabs(writhe) < _writhe_threshold or k == final_particle_index) {
					on_peak = false;
					//printf("Getting off peak at position %d\n",k);
					if(peak_position > _last_particle_index) {
						peak_position -= (_last_particle_index - _first_particle_index + 1);
					}
					result += std::string(temp);
					// compute the size of the plectoneme
					if(_print_size or _print_left_right) {
						if((peak_position - _minimum_plectoneme_size) >= 0 && (peak_position + _minimum_plectoneme_size) < _config_info->N()) {		//TODO: probably remove this if statement - the pointers should keep track of things by checking for P_VIRTUAL
						//--old way
						//BaseParticle *left = p[peak_position - _minimum_plectoneme_size];
						//BaseParticle *right = p[peak_position + _minimum_plectoneme_size];
						//--
							BaseParticle *left = p[peak_position];
							BaseParticle *right = p[peak_position];
							for(int l = 0; l < _minimum_plectoneme_size; l++) {
								if(left->n3 != P_VIRTUAL) left = left->n3;
								if(right->n5 != P_VIRTUAL) right = right->n5;
							}
							bool done = false;

							while(!done) {
								BaseParticle *new_left = left;
								BaseParticle *new_right = right;
								number min_distance = 1e6;
								//--old if statement
								//for(int pleft = 1; pleft < _size_outer_threshold && (left->index - pleft) >= 0; pleft++) {
								//	BaseParticle *curr_left = p[left->index - pleft];
								BaseParticle *curr_left = left->n3;
								for(int pleft = 1; pleft < _size_outer_threshold && curr_left != P_VIRTUAL; pleft++) {
									//--old if statement
									//for(int pright = 1; pright < _size_outer_threshold && (right->index + pright) < *_config_info->N; pright++) {
									//BaseParticle *curr_right = p[right->index + pright];
									BaseParticle *curr_right = right->n5;
									for(int pright = 1; pright < _size_outer_threshold && curr_right != P_VIRTUAL; pright++) {
										number curr_dist = (curr_right->pos - curr_left->pos).module();
										if(curr_dist < min_distance) {
											min_distance = curr_dist;
											new_left = curr_left;
											new_right = curr_right;
										}
										curr_right = curr_right->n5;
									}
									curr_left = curr_left->n3;
								}
								if(min_distance < _contact_threshold) {
									left = new_left;
									right = new_right;
									// this will be triggered for circular molecules
									if(left == right or left->n3 == right) done = true;
								}
								else done = true;
							}

							//--old size
							//int size = right->index - left->index;
							int size = 0, max_size = _config_info->N() + 10;
							BaseParticle * step_counter = left;
							for(int l = 1; (l <= max_size and step_counter != right); l++) {
								size = l;
								step_counter = step_counter->n5;
							}
							if(size == 0) size = _config_info->N();
							//printf("Just computed: peak %d size %d left %d right %d\n",peak_position, size, left->index, right->index);
							if(size == max_size) {
								OX_LOG(Logger::LOG_INFO,"Obsevrable writhe: problem with plectoneme on position %d: size is equal to max_size (%d). Left = %d Right = %d.",peak_position, size, left->index, right->index);

							}

							// we don't print info about plectonemes that are smaller than _minimum_plectoneme_size
							if(size > _minimum_plectoneme_size * 2) {
								if(_print_size) result += Utils::sformat(" %d %d", peak_position, size);
								else if(_print_left_right) result += Utils::sformat(" %d %d %d", peak_position, left->index, right->index);

								// print the spatial position of the particle on the tip of the plectoneme
								if(_print_space_pos) {
									BaseParticle *p_tip = p[peak_position];
									result += Utils::sformat(" %lf %lf %lf", p_tip->pos.x, p_tip->pos.y, p_tip->pos.z);
								}
								// print the bending angle averaged over a few segments around the tip bead
								if(_bending_angle_number_segments > 0) {
									int offset = 0;
									number average_angle = 0;
									for(int i = 0; i < _bending_angle_number_segments; i++) {
										int i_1 = peak_position + offset;
										int i_2 = i_1 + 1;
										number angle = LRACOS(p[i_1]->orientationT.v1 * p[i_2]->orientationT.v1);
										printf("%d %g %d %g %g\n", i_1, p[i_1]->orientationT.v1.x, i_2, p[i_2]->orientationT.v1.x, angle);
										average_angle += angle;
										// this offset starts at 0, and then goes -1, 1, -2, 2, -3, 3, etc.
										offset = ((offset >= 0) * 1 + offset) * (-1);
									}
									average_angle /= _bending_angle_number_segments;
									result += Utils::sformat(" %lf", average_angle);
								}
								result += "\t";
							}
						}
					}
					else result += Utils::sformat("%d\t", peak_position);
				}
			}
		}
		else {
			sprintf(temp, "%14.14lf\n", writhe);
			result += std::string(temp);
		}
	}
	if(on_peak) throw oxDNAException("Writhe.cpp never got off a peak on time %d. Most likely a bug!!!", time);
	if(stop) throw oxDNAException("Dying badly because of problems with the observable writhe. That's not supposed to happen!");

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
