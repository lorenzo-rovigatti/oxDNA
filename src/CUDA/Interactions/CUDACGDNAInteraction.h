/*
 * CUDACGDNAInteraction.h
 *
 *  Created on: 7/27/22
 *      Author: jonah
 */

#ifndef CUDACGDNAINTERACTION_H_
#define CUDACGDNAINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/CGDNAInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model with ANM protein model as implemented in DNANMInteraction.
 */
class CUDACGDNAInteraction: public CUDABaseInteraction, public CGDNAInteraction {
public:

    bool _use_debye_huckel;
    bool _use_oxDNA2_coaxial_stacking;
    bool _use_oxDNA2_FENE;
    // copied from DNA2Interaction.h (CPU) (and change number -> float), the least bad way of doing things
    float _salt_concentration;
    bool _debye_huckel_half_charged_ends;
    float _debye_huckel_prefactor;
    float _debye_huckel_lambdafactor;

    float _debye_huckel_RC; // this is the maximum interaction distance between backbones to interact with DH
    float _debye_huckel_RHIGH; // distance after which the potential is replaced by a quadratic cut-off
    float _debye_huckel_B; // prefactor of the quadratic cut-off
    float _minus_kappa;

    float _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_sqr_rcut;
    float _pro_base_sigma, _pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_sqr_rcut;
    float _pro_rstar, _pro_b, _pro_rcut, _pro_sqr_rcut, _pro_sigma;

    c_number *_pro_spring_eqdist, *_pro_spring_potential; //Temp arrays for parameter storage
    c_number *_gs_spring_eqdist, *_gs_spring_potential; //Temp arrays for parameter storage

    c_number *_h_gs_other_exc_vol_params;
    c_number *_d_gs_other_exc_vol_params;

    c_number *_h_gs_gs_exc_vol_params;
    c_number *_d_gs_gs_exc_vol_params;

    int _pro_offset;
    int _gs_offset;

    size_t _pro_spring_param_size_number, _gs_spring_param_size_number;

    bool _read_par;

    //compressed parameter arrays
    c_number *_h_pro_aff_gamma, *_d_pro_aff_gamma;
    c_number *_h_pro_aff_eqdist, *_d_pro_aff_eqdist;

    c_number *_h_gs_aff_gamma, *_d_gs_aff_gamma;
    c_number *_h_gs_aff_eqdist, *_d_gs_aff_eqdist;

    int *_h_pro_affected, *_d_pro_affected;
    int *_h_gs_affected, *_d_gs_affected;

    int *_pro_affected_len, *_h_pro_affected_indx, *_d_pro_affected_indx;
    int *_gs_affected_len, *_h_gs_affected_indx, *_d_gs_affected_indx;

    explicit CUDACGDNAInteraction(bool btp);
    virtual ~CUDACGDNAInteraction();

    void get_settings(input_file &inp);
    void cuda_init(int N);
    c_number get_cuda_rcut() {
        return this->get_rcut();
    }
    virtual void _on_T_update();

    virtual void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

#endif /* CUDACGDNAINTERACTION_H_ */
