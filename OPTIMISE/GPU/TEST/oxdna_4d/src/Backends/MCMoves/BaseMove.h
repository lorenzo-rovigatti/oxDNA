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
class BaseMove {
	protected:
		/// number of attempted moves
		llint _attempted;

		/// number of accepted moves
		llint _accepted;

		/// Temperature at which the simulation is run
		number _T;

		/// Information on the system.
		std::shared_ptr<ConfigInfo> _Info;

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

		virtual void _on_T_update();

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
		virtual void init();

		/// fills a string with the info
		virtual void log_parameters();

		/// helper function to compute the energy of each particle; internally, it calls the interaction;
		number particle_energy (BaseParticle *p);

		/// helper function to compute the system energy; internally, it calls the interaction
		number system_energy();

		/// method that applies the move to the system. Each child class must have it.
		virtual void apply(llint curr_step) = 0;

		/// method that gets the ratio of accepted moves
		virtual double get_acceptance() {
			if (_attempted > 0) {
				return _accepted / (double) _attempted;
			}
			else {
				return 0.;
			}
		}
};

using MovePtr = std::shared_ptr<BaseMove>;

#endif
