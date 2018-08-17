/**
 * @file    MCMovePatchyShape.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#ifndef MCMOVEPS_H_
#define MCMOVEPS_H_
#define MCMOVE_CUSTOM

#include "../../../../../src/Backends/MCMoves/BaseMove.h"
#include "../../Interactions/PatchyShapeInteraction.h"
#include "../../Particles/PatchyShapeParticle.h"

template<typename number>
class MCMovePatchyShape : public BaseMove<number> {
	protected:
		number _delta;
		number _delta_rotation;

		LR_matrix<number> _orientation_old, _orientationT_old;
		LR_vector<number> pos_old;

		number _verlet_skin;
	
	public:
		MCMovePatchyShape();
		virtual ~MCMovePatchyShape();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init(void);

};


extern "C" BaseMove<float> * make_MCMovePatchyShape_float()   { return new MCMovePatchyShape<float> () ; }
extern "C" BaseMove<double> * make_MCMovePatchyShape_double() {return new MCMovePatchyShape<double> () ; }

#endif // MCMovePatchyShape.h



