/**
 * @brief Hybrid DNA/ANM Model
 *
 *
 *
 * Jonah Feb. 2021
 * 
 * To use the oxDNA2 model with ANM Protein, set
 *
 * interaction_type = CGDNA
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (set parameter file for protein component, set to none for DNA only sims)
 massfile = string (set massfile for simulations w/ different massed items)
 */

#ifndef CGDNA_INTERACTION_H
#define CGDNA_INTERACTION_H

#include "DNANMInteraction.h"


class CGDNAInteraction: public DNANMInteraction {

protected:
    int ngs; // number of generic sphere particles
    int npep; // number of protein strands
    int ngstrands; // number of gs strands

    void load_extra_file(std::string &massfile); // always the mass/ radii file
    std::map<int, number> radii;

    number _gs_exc_vol_stiffness;

    // indexed by type
    number (CGDNAInteraction::*_interaction_matrix_bonded[3])(BaseParticle *, BaseParticle *, bool, bool){
        &CGDNAInteraction::_dna_bonded,
        &CGDNAInteraction::_pro_bonded,
        &CGDNAInteraction::_gs_bonded
    };
    // indexed by type
    number (CGDNAInteraction::*_interaction_matrix_nonbonded[9]) (BaseParticle *, BaseParticle *, bool, bool){
            &CGDNAInteraction::_dna_nonbonded,  // 0 0
            &CGDNAInteraction::_dna_pro_nonbonded,  // 0 1
            &CGDNAInteraction::_dna_gs_nonbonded,  // 0 2
            &CGDNAInteraction::_dna_pro_nonbonded,  // 1 0
            &CGDNAInteraction::_pro_nonbonded,  // 1 1
            &CGDNAInteraction::_pro_gs_nonbonded,  // 1 2
            &CGDNAInteraction::_dna_gs_nonbonded,  // 2 0
            &CGDNAInteraction::_pro_gs_nonbonded,  // 2 1
            &CGDNAInteraction::_gs_nonbonded,  // 2 2
    };


    bool _fill_in_constants(number interaction_radius, number& sigma, number& rstar, number &b, number &rc)
    {
        rstar = interaction_radius;
        sigma = rstar + 0.05f;
        b = _calc_b(sigma,rstar);
        rc = _calc_rc(sigma,rstar);
        return _check_repulsion_smoothness(sigma,rstar,b,rc);
    }

    bool _fill_in_constants_quart(number interaction_radius, number& sigma, number& rstar, number &b, number &rc)
    {
        rstar = interaction_radius;
        sigma = rstar + 0.001f;
        b = _calc_b_quart(sigma,rstar);
        rc = _calc_rc_quart(sigma,rstar);
        return _check_repulsion_smoothness_quart(sigma,rstar,b,rc);
    }

    static number pow6(number x){
        return SQR(x)*SQR(x)*SQR(x);
    }

    static number pow14(number x){
        return SQR(x)*SQR(x)*SQR(x)*SQR(x)*SQR(x)*SQR(x)*SQR(x);
    }

    static number _calc_b(number sigma,number R)
    {
        return 36.0f * pow6(sigma)  * SQR((- pow6(R) + 2.0f * pow6(sigma)))   / (pow14(R) *(-pow6(R)+ pow6(sigma)));
    }

    static number _calc_b_quart(number sigma,number R)
    {
        return 81.0f * pow6(sigma)  * SQR(SQR((- pow6(R) + 2.0f * pow6(sigma))))   / (4.f * pow14(R) * CUB(-pow6(R)+ pow6(sigma)));
    }

    static number _calc_rc_quart(number sigma,number R)
    {
        return R *(5.0f *  pow6(R) -8.0f * pow6(sigma)  ) /(3.0f * (pow6(R) -2.0f *  pow6(sigma) ));
    }

    static number _calc_rc(number sigma,number R)
    {
        return R *(4.0f *  pow6(R) -7.0f * pow6(sigma)  ) /(3.0f * (pow6(R) -2.0f *  pow6(sigma) ));
    }

    bool _check_repulsion_smoothness(number sigma,number rstar, number b, number rc)
    {
        /*TODO IMPLEMENT ME NICER!*/
//        return true;
        // printf
        // printf("##### STARTING checking for %f\n",rstar);
        //double tolerance = 0.2;
        LR_vector rr(0,0,rstar-0.2);
        LR_vector forcer(0,0,0);
        double stiff = 1.f;
        number old_energy =  _repulsive_lj(rr,forcer, true, sigma, b, rstar, rc, stiff);
        number oldforce = forcer.norm();
        for (double x = rstar-0.2; x < rc+0.01; x += 0.0001)
        {
            LR_vector r(0,0,x);
            LR_vector force(0,0,0);
            number energy_new = _repulsive_lj(r,force, true, sigma, b, rstar, rc, stiff);
            number force_new = force.norm();
            //printf("#### %f %f %f \n",x,energy_new,force_new);
            //if( fabs(old_energy - energy_new) > tolerance || fabs(force_new - oldforce) > tolerance || old_energy < energy_new)
            if(  old_energy < energy_new || oldforce < force_new)
            {
                throw oxDNAException("Non-monotonous change in the repulsion potential for distance  r = %f",x);
                return false;
            }
            else
            {
                old_energy = energy_new;
                oldforce = force_new;
            }
        }
        return true;
    };

    bool _check_repulsion_smoothness_quart(number sigma,number rstar, number b, number rc)
    {
        /*TODO IMPLEMENT ME NICER!*/
//        return true;
        // printf
        // printf("##### STARTING checking for %f\n",rstar);
        //double tolerance = 0.2;
        LR_vector rr(0,0,rstar-0.2);
        LR_vector forcer(0,0,0);
        double stiff = 1.f;
        number old_energy =  _repulsive_lj_quart(rr,forcer, true, sigma, b, rstar, rc, stiff);
        number oldforce = forcer.norm();
        for (double x = rstar-0.2; x < rc+0.01; x += 0.0001)
        {
            LR_vector r(0,0,x);
            LR_vector force(0,0,0);
            number energy_new = _repulsive_lj_quart(r,force, true, sigma, b, rstar, rc, stiff);
            number force_new = force.norm();
            //printf("#### %f %f %f \n",x,energy_new,force_new);
            //if( fabs(old_energy - energy_new) > tolerance || fabs(force_new - oldforce) > tolerance || old_energy < energy_new)
            if(  old_energy < energy_new || oldforce < force_new)
            {
                throw oxDNAException("Non-monotonous change in the repulsion potential for distance  r = %f",x);
                return false;
            }
            else
            {
                old_energy = energy_new;
                oldforce = force_new;
            }
        }
        return true;
    };


public:
    enum {
        GS_EXC_VOL = 12,
        GS_DNA_EXC_VOL = 13,
        GS_PRO_EXC_VOL = 14
        //1-11 are assigned in parent classes DNANM->DNA2->DNA
    };

    int gs_subtype_num;
    std::array<int, 3> topology_order;
    number *_gs_other_exc_vol_params;
    number *_gs_gs_exc_vol_params;

    explicit CGDNAInteraction(bool btp); //btn controls whether bending/torsion potential is applied
    ~CGDNAInteraction() override;
    void get_settings(input_file &inp) override;
    void allocate_particles(std::vector<BaseParticle *> &particles) override;
    void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
    number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    void check_input_sanity(std::vector<BaseParticle *> &particles) override;
    void init() override;
    number _repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
    number _repulsive_lj_quart(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
    number _gs_dna_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _gs_pro_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _gs_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _dna_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _pro_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _pro_gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _dna_gs_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _dna_pro_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _gs_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _pro_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    number _dna_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    int particle_id(int type);
};

#endif /* CGDNA_INTERACTION_H */
