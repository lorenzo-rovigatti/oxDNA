/**
 * @file    VMMC.h
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */
#ifndef VMMC_H_
#define VMMC_H_

#include "BaseMove.h"

class VMMC : public BaseMove {
	struct movestr {
	    int seed, type, seed_strand_id;
	    LR_vector t;
	    LR_matrix R;
	    LR_matrix Rt;
	};

	protected:
		number _delta_tras;
		number _delta_rot;
		LR_vector pos_old;

		number _max_move_size, _max_move_size_sqr;
		int _max_cluster_size;
		
		number _build_cluster();

		std::vector<bool> _inclust();

		std::vector<BaseParticle *>  _particles_old;
		std::vector<int> _clust;
		
		inline void _store_particle(BaseParticle *src);
		inline void _restore_particle(BaseParticle *src);
		inline void _re_restore_particle(BaseParticle *src);

		inline void _move_particle(movestr *moveptr, BaseParticle *p);

		number VMMC_link(double E_new, double E_old) { return (1. - exp((1. / this->_T) * (E_old - E_new)));}
		inline number _next_rand () {return drand48();}

		number build_cluster (movestr *moveptr, int maxsize);
		number build_cluster_old (movestr *moveptr, int maxsize);
	
	public:
		VMMC();
		virtual ~VMMC();

		void apply (llint curr_step);
		virtual void get_settings(input_file &inp, input_file &sim_inp);
		virtual void init();
};
#endif // VMMC_H_
