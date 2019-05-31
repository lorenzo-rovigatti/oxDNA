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
		llint _attempted;

		/// number of accepted moves
		llint _accepted;

		/// Temperature at which the simulation is run
		number _T;

		/// Information on the system.
		ConfigInfo<number> * _Info;

		/// should the move compute the energy before acting (default=true, use false only for hard potentials with no attraction)
		bool _compute_energy_before;

		/// submoves
		vector<BaseMove *> _submoves;

		/// whether to adjust moves during equilibration
		bool _adjust_moves;

		/// restrict the move to this particle type
		int _restrict_to_type;

		/// number of equilibration steps
		llint _equilibration_steps;

		/// target acceptance rate for the move (used only when adjust_moves=true)
		number _target_acc_rate;

		/// auxiliary factors to adjust acceptance, used internally
		number _rej_fact, _acc_fact;

		/// type of the move
		std::string _name;

	public:
		BaseMove();

		virtual ~BaseMove();

		/// relative probability with which the move is attempted
		number prob;

		/// array of moved particles
		vector<int> moved;

		/// method to get the move settings from the input file
		virtual void get_settings(input_file &inp, input_file &sim_inp);

		/// method to do any initialization
		virtual void init ();

		/// fills a string with the info
		virtual void log_parameters();

		/// helper function to compute the energy of each particle; internally, it calls the interaction;
		number particle_energy (BaseParticle<number> *p);

		/// helper function to compute the system energy; internally, it calls the interaction
		number system_energy ();

		/// method that applies the move to the system. Each child class must have it.
		virtual void apply (llint curr_step) = 0;

		/// method that gets the ratio of accepted moves
		virtual double get_acceptance() {
			if (_attempted > 0) return _accepted / (double) _attempted;
			else return 0.;
		}
};

template<typename number>
BaseMove<number>::BaseMove ():  _Info(&ConfigInfo<number>::ref_instance()) {
	_attempted = 0;
	_accepted = 0;
	_T = -1.;
	//_Info = Info;
	prob = (number) 1.f;
	_target_acc_rate = 0.25;
	_equilibration_steps = 0;
	_adjust_moves = false;
	_compute_energy_before = true;
	_restrict_to_type = -1;
}

template<typename number>
BaseMove<number>::~BaseMove () {

}

template<typename number>
void BaseMove<number>::get_settings (input_file &inp, input_file &sim_inp) {
	char raw_T[1024];
	getInputString(&sim_inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);

	getInputString (&inp, "type", _name, 1);
	getInputNumber (&inp, "prob", &prob, 0);
	getInputLLInt(&sim_inp, "equilibration_steps", &_equilibration_steps, 0);
	getInputBool(&inp, "adjust_moves", &_adjust_moves, 0);
	getInputNumber(&inp, "target_acc_rate", &_target_acc_rate, 0);
	getInputBool(&inp, "compute_energy_before", &_compute_energy_before, 0);
	getInputInt(&inp, "restrict_to_type", &_restrict_to_type, 0);
}

template<typename number>
void BaseMove<number>::init() {
	_rej_fact = 1.001;
	_acc_fact = 1. + 0.001 * (_target_acc_rate - 1.) / (0.001 - _target_acc_rate);
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h) BaseMove init for move type %s;", _name.c_str());
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h)               (b, a) = (%g, %g), target_acc_rate=%g, adjust_moves=%d,", (_rej_fact - 1.), (_acc_fact - 1.), _target_acc_rate, _adjust_moves);
	OX_LOG(Logger::LOG_INFO, "(BaseMove.h)               compute_energy_before=%d, restrict_to_type=%d", _compute_energy_before, _restrict_to_type);
}

template<typename number>
void BaseMove<number>::log_parameters() {
	OX_LOG(Logger::LOG_INFO, "Logging parameters for MC2 Move with type = %s", _name.c_str());
	OX_LOG(Logger::LOG_INFO, "\tprob = %g", prob);
	OX_LOG(Logger::LOG_INFO, "\tadjust_moves = %d", int(_adjust_moves));
	OX_LOG(Logger::LOG_INFO, "\ttarget_acc_rate = %g", _target_acc_rate);
	OX_LOG(Logger::LOG_INFO, "\tcompute_energy_before = %d", int(_compute_energy_before));
	OX_LOG(Logger::LOG_INFO, "\trestrict_to_type = %d", int(_restrict_to_type));
}

template <typename number>
number BaseMove<number>::particle_energy(BaseParticle<number> * p) {
	number res = (number) 0.f;
	typename vector<ParticlePair<number> >::iterator it = p->affected.begin();
	for(; it != p->affected.end(); it++) {
		number de = _Info->interaction->pair_interaction_bonded(it->first, it->second);
		res += de;
	}
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
