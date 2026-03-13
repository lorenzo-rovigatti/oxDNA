/*
 * CUDADNANMInteraction.h
 *
 *  Created on: 6/11/19
 *      Author: jonah
 */

#ifndef CUDADNANMINTERACTION_H_
#define CUDADNANMINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNANMInteraction.h"

/**
 * @brief CUDA implementation of the oxDNA model with ANM protein model as implemented in DNANMInteraction.
 */
class CUDADNANMInteraction: public CUDABaseInteraction, public DNANMInteraction {
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

    int *_d_is_strand_end = nullptr;

    float _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_sqr_rcut;
    float _pro_base_sigma, _pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_sqr_rcut;
    float _pro_rstar, _pro_b, _pro_rcut, _pro_sqr_rcut, _pro_sigma;
    float _ktor, _kbend;

    c_number *_spring_eqdist, *_spring_potential; //Temp arrays for parameter storage

    c_number *_d_ang_params, *_h_ang_params;
    c_number *_d_ang_kbkt, *_h_ang_kbkt;
    size_t _spring_param_size_number, _ang_param_size;

    bool _read_par;

    //compressed parameter arrays
    c_number *_h_aff_gamma, *_d_aff_gamma;
    c_number *_h_aff_eqdist, *_d_aff_eqdist;
    int *_h_affected, *_d_affected;
    int *_affected_len, *_h_affected_indx, *_d_affected_indx;

    int offset; //Will only come into play if proteins are after dna in topology file (particle id wise). Adjusts proteins index for the spring parameter arrays
    explicit CUDADNANMInteraction(bool btp);
    virtual ~CUDADNANMInteraction();

    void get_settings(input_file &inp);
    void cuda_init(int N);
    c_number get_cuda_rcut() {
        return this->get_rcut();
    }
    virtual void _on_T_update();
    virtual void _init_strand_ends(LR_bonds *d_bonds);

    virtual void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box);
};

#endif /* CUDADNANMINTERACTION_H_ */
