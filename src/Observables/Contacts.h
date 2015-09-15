/*
 * Contacts.h
 *
 *  Created on: 31 Aug,  2015
 *      Author: Ferdinando Randisi
 */

#ifndef CONTACTS_H_
#define CONTACTS_H_

#include "BaseObservable.h"
/**
 * @brief Outputs a list of all the particle pairs that are within a given distance of each other.
 *        This observable has been thought to be used for a singly-bound polymer model (e.g. ssDNA, or a double strand of DNA modelled with the TEP model) but should in principle being usable for other systems as well.

 * To use this observable, use type = contacts.

 * This observable takes 5 optional arguments
 * @verbatim

first_particle_index = <int> (defaults to 0. index of the first particle to consider. All the particles coming before this one will be ignored.)

last_particle_index = <int> (defaults to the index of the first-but last bead in the same strand as the first particle. Therefore, if I have a strand of N beads, the last one will be the one with index N-2. This is because the last bead is atypical in the TEP model (e.g. it's aligned with the vector before it rather than the one in front of it.). index of the last particle to consider. All the particles coming before this one will be ignored.)

neighbours_to_ignore = <int> (defalts to 1. Number of neighbours to ignore before-after each particle. E.g., if equals to 1, contacts between given first-neighbours will never be reported, if equals to 2, contacts between second neighbours will never be reported, etc).

contact_distance = <number> (defaults to 1. A contact is defined if the centers of mass of the particles is lower than this value).

only_outermost_contacts = <bool> (defaults to false. if true, contacts nested within other contacts will not be reported. E.g. if the i-th monomer is linked to both the i-1-th and the i+1-th monomer, and the contacts are 10-40, 10-25, 13-32, 12-48 and 45-60, only 10-40, 12-48 and 45-60 will be reported, since 10-25 and 13-32 are both nested inside 10-40. This is only get one result per plectoneme. Whatch out though, since this will report clashes between a plectoneme and the chain/other plectonemes. Telling a plectoneme and a plectoneme contact just by using the contact map might be non-banal ).


TODO: this can be made more efficient implementing adjacency lists. I'm not doing it right now because I don't think I'll need it.

@endverbatim
 */

template<typename number>
class Contacts : public BaseObservable<number> {
private:
	// arguments
	int _first_particle_index;
	int _last_particle_index;
	int _neighbours_to_ignore;
	bool _only_outermost_contacts;
	number _contact_distance;

public:
	Contacts();
	virtual ~Contacts();

	virtual void init(ConfigInfo<number> &config_info);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* CONTACTS_H_ */
