/**
 * @file    MoleculeVolumeMove.h
 * @date    21/apr/2020
 * @author  lorenzo
 */

#ifndef MOLECULE_VOLUME_MOVE_H_
#define MOLECULE_VOLUME_MOVE_H_

#include "BaseMove.h"

class MoleculeVolumeMove : public BaseMove {
	protected:
		number _delta = 0.;

		number _verlet_skin = 0.;
		number _P = 0.;
		bool _isotropic = true;

	public:
		MoleculeVolumeMove();
		virtual ~MoleculeVolumeMove();

		void apply (llint curr_step);
		virtual void init ();
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void log_parameters();
};
#endif // MOLECULE_VOLUME_MOVE_H_
