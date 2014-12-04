/**
 * @file    BaseMove.h
 * @date    30/apr/2014
 * @author  flavio
 *
 */

#ifndef MCMOVES_H_
#define MCMOVES_H_

#include "../../Utilities/Logger.h"
#include "../../Particles/BaseParticle.h"
#include "../../Observables/BaseObservable.h"

using namespace std;

/**
 * @brief Abstract class defining the MC moves. All the other moves inherit from this one.
 *
 */
template<typename number>
class BaseMove {
	protected:
		/// number of attempted moves
		unsigned int _attempted;

		/// number of accepted moves
		unsigned int _accepted;

		/// Temperature at which the simulation is run
		number _T;

		/// Information on the system.
		ConfigInfo<number> * _Info;

		/// submoves
		vector<BaseMove *> _submoves;

	
	public:
		BaseMove(ConfigInfo<number> * Info);
		virtual ~BaseMove();	

		/// relative probability with which the move is attempted
		number prob;

		/// array of moved particles
		vector<int> moved;

		/// method to get the move settings from the input file
		virtual void get_settings(input_file &inp, input_file &sim_inp);

		/// method to do any initialization
		virtual void init ();

		/// helper function to compute the energy of each particle; internally, it calls the interaction;
		number particle_energy (BaseParticle<number> *p);

		/// helper function to compute the system energy; internally, it calls the interaction
		number system_energy ();

		/// method that applies the move to the system. Each child class must have it.
		virtual void apply (llint curr_step) = 0;
};

template<typename number>
BaseMove<number>::BaseMove (ConfigInfo<number> * Info) {
	_attempted = 0;
	_accepted = 0;
	_T = -1.;
	_Info = Info;
	prob = (number) -1.f;
}

template<typename number>
BaseMove<number>::~BaseMove () {

}

template<typename number>
void BaseMove<number>::get_settings (input_file &inp, input_file &sim_inp) {
	getInputNumber (&sim_inp, "T", &_T, 1);
	//OX_LOG(Logger::LOG_INFO, "(BaseMoves.cpp) Here is T %g", _T);
}

template<typename number>
void BaseMove<number>::init() {
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h) BaseMove generic init...");
}

template <typename number>
number BaseMove<number>::particle_energy (BaseParticle<number> * p) {
	number res = (number) 0.f;
	if(p->n3 != P_VIRTUAL) res += _Info->interaction->pair_interaction_bonded(p, p->n3);
	if(p->n5 != P_VIRTUAL) res += _Info->interaction->pair_interaction_bonded(p, p->n5);

	if (_Info->interaction->get_is_infinite() == true) return (number) 1.e12;

	std::vector<BaseParticle<number> *> neighs = _Info->lists->get_neigh_list(p);
	for(unsigned int n = 0; n < neighs.size(); n++) {
		BaseParticle<number> *q = neighs[n];
		res += _Info->interaction->pair_interaction_nonbonded(p, q);
		if(_Info->interaction->get_is_infinite() == true) {
			return (number) 1.e12;
		}
	}
	return res;
}

template <typename number>
number BaseMove<number>::system_energy() {
	number res = (number) 0.f;

	for (int i = 0; i < *_Info->N; i ++) {
		BaseParticle<number> *p = _Info->particles[i];
		if(p->n3 != P_VIRTUAL) res += _Info->interaction->pair_interaction_bonded(p, p->n3);
		// we omit E(p,p->n5) because it gets counted as E(q, q->n3);
		std::vector<BaseParticle<number> *> neighs = _Info->lists->get_neigh_list(p);
		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle<number> *q = neighs[n];
			if (p->index < q->index) {
				res += _Info->interaction->pair_interaction_nonbonded(p, q);
				if(_Info->interaction->get_is_infinite() == true) return (number) 1.e12;
			}
		}
	}

	return res;
}
#endif
