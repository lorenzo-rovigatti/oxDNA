#ifndef CUDACCGINTERACTION_H_
#define CUDACCGINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/CCGInteraction.h"

class CUDACCGInteraction: public CUDABaseInteraction, public CCGInteraction{
public:
    CUDACCGInteraction();
    virtual ~CUDACCGInteraction();
    void get_settings(input_file &inp);
    void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	};
    virtual void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
    
};

#endif