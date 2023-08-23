/**
 * @brief Hybrid RNA/ANM Model
 *
 *
 *
 * Jonah Feb. 2021
 * 
 * To use the oxRNA2 model with ANM Protein, set
 *
 * interaction_type = RNANM
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (set parameter file for protein component, set to none for RNA only sims)
 massfile = string (set massfile for simulations w/ different massed items)
 */

#ifndef RNANM_INTERACTION_H
#define RNANM_INTERACTION_H

#include "RNAInteraction2.h"


class RNANMInteraction: public RNA2Interaction {

protected:
    int nrna;//How many particles of RNA type: Used in allocate_particles
    int npro;//How many bonds b/t different particle types
    int nrnas;//Number of Strands that are rna, rest assumed to be protein
    int _firststrand; //+ for rna, - for protein

    // excluded volume parameters for quartic LJ on backbone/protein + base/protein
    number _pro_backbone_sigma, _pro_backbone_rstar, _pro_backbone_b, _pro_backbone_rcut, _pro_backbone_stiffness;
    number _pro_base_sigma,_pro_base_rstar, _pro_base_b, _pro_base_rcut, _pro_base_stiffness;
    std::map<std::pair<int, int>, double> _rknot; //Both maps used just as they are in ACInteraction
    std::map<std::pair<int, int>, std::pair<char, double> > _potential;
    std::map<int, number> masses;
    bool _parameter_kbkt; //Controls whether kb/kt values are global or read from parameter file
    std::map<int, std::pair<double, double> > _ang_stiff; // Stores per particle pair, kb kt values
    number _k_bend, _k_tor;
    std::map<int, std::vector <double> > _ang_vals;
    number _pro_sigma, _pro_rstar, _pro_b, _pro_rcut; // Protein-protein quartic LJ params same as in ANM
    number _pro_base_sqr_rcut, _pro_backbone_sqr_rcut, _pro_sqr_rcut, _pro_rna_sqr_rcut;

public:
    enum {
        SPRING = 8,
        PRO_EXC_VOL = 9,
        PRO_RNA_EXC_VOL = 10,
        PRO_ANG_POT = 10
        //Assigned 8 9 and 10 so it won't overwrite the already existing RNA function pointers in the _int_map
    };

    char _parameterfile[500];
    std::string _massfile;
    bool _angular;

    explicit RNANMInteraction(bool btp);
    virtual void load_massfile(std::string &filename);
    virtual ~RNANMInteraction();
    void get_settings(input_file &inp) override;
    void allocate_particles(std::vector<BaseParticle *> &particles) override;
    void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
//    void load_masses(std::string &massfile);
    number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number _protein_rna_exc_volume(BaseParticle *p,BaseParticle *q, bool compute_r, bool update_forces);
    number _protein_rna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
    void check_input_sanity(std::vector<BaseParticle *> &particles) override;
    void init() override;
    number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
    number _protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_ang_pot(BaseParticle *p, BaseParticle*q, bool compute_r, bool update_forces);
};

#endif /* RNANM_INTERACTION_H */
