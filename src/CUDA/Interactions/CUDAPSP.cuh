// CUDA Patchy Spring Particles

#ifndef CUDA_PSP_CUH
#define CUDA_PSP_CUH

#include "CUDABaseInteraction.h"
#include "../../Interactions/PSPInteraction.h" //cupu interactions
#include "../cuda_utils/CUDA_lr_common.cuh"
#include "../Lists/CUDASimpleVerletList.h"
#include "../Lists/CUDANoList.h"

#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/transform.h>

class CUDAPSPInteraction: public CUDABaseInteraction, public PSPInteraction {
public:
    

    CUDAPSPInteraction();
    virtual ~CUDAPSPInteraction();
    void get_settings(input_file &inp) override;
    void cuda_init(int N) override;
    c_number get_cuda_rcut() {
        return this->get_rcut();
    }
    void compute_forces(CUDABaseList *lists,c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torque, LR_bonds *d_bonds, CUDABox *d_box);

};

#endif // CUDA_PSP_CUH