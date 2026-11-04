/**
 * @brief Hybrid DNA/ANM Model
 *
 *
 *
 * Jonah,Subhajit Feb. 2026
 * 
 * To use the oxDNA2 model with ANM Protein, set
 *
 * interaction_type = DNANM
 *
 * in the input file
 *
 * Input options:
 *
 * @verbatim
 parfile = string (set parameter file for protein component, set to none for DNA only sims)
 massfile = string (set massfile for simulations w/ different massed items)
 topology = mixed legacy topology format:
            <N> <N_strands> <ndna> <npro> <ndnas>
 */

#ifndef DNANM_INTERACTION_H
#define DNANM_INTERACTION_H

#include "DNA2Interaction.h"
#include <vector>


class DNANMInteraction: public DNA2Interaction {

protected:
    int ndna;//How many particles of DNA type: Used in allocate_particles
    int npro;//How many bonds b/t different particle types
    int ndnas;//Number of Strands that are dna, rest assumed to be protein
    int _firststrand; //+ for dna, - for protein

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
    number _pro_base_sqr_rcut, _pro_backbone_sqr_rcut, _pro_sqr_rcut, _pro_dna_sqr_rcut;

    number _pro_pro_sqr_rcut;
    bool _protein_dna_electrostatics;
    bool _protein_protein_electrostatics;
    bool _exclude_intragroup_electrostatics;

    number _protein_charge_scale;
    number _protein_dna_dh_scale;
    number _protein_protein_dh_scale;
    number _dna_effective_charge_e;

    std::string _protein_charge_file;

    // Indexed in global particle space.
    std::vector<number> _protein_charge_reduced;
    std::vector<int> _protein_group;

    // Explicit mapping from protein-local charge-file indices
    // to the mixed DNANM/DNACT topology.
    std::vector<int> _protein_local_to_global;
    std::vector<int> _protein_topology_strand;
    std::vector<char> _protein_residue_code;

    int _protein_charge_format_version;
    std::string _protein_charge_model;
    std::string _protein_charge_source;

    void _prepare_protein_charge_topology_metadata(int total_particles);
    void _load_protein_charge_file(int total_particles);
    void _update_electrostatic_cutoffs();

    number _scaled_debye_huckel(
        const LR_vector &r,
        number charge_product,
        LR_vector &force,
        bool update_forces
    );

    number _protein_dna_debye_huckel(
        BaseParticle *p,
        BaseParticle *q,
        bool compute_r,
        bool update_forces
    );

    number _protein_protein_debye_huckel(
        BaseParticle *p,
        BaseParticle *q,
        bool compute_r,
        bool update_forces
    );

public:
    enum {
        SPRING = 8,
        PRO_EXC_VOL = 9,
        PRO_DNA_EXC_VOL = 10,
        PRO_ANG_POT = 11,
        PRO_DNA_DH = 12,
        PRO_PRO_DH = 13
        //Assigned 8 9 and 10 so it won't overwrite the already existing DNA function pointers in the _int_map
    };

    char _parameterfile[512];
    std::string _massfile;
    bool _angular;

    int get_id(int btype);
    explicit DNANMInteraction(bool btp); //btn controls whether bending/torsion potential is applied
    virtual void load_massfile(std::string &filename);
    ~DNANMInteraction() override;
    void get_settings(input_file &inp) override;
    void allocate_particles(std::vector<BaseParticle *> &particles) override;
    void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
    void read_parameter_file(std::vector<BaseParticle *> &particles);

    number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q,
                                 bool compute_r = true, bool update_forces = false) override;
    number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) override;
    number _protein_dna_exc_volume(BaseParticle *p,BaseParticle *q, bool compute_r, bool update_forces);
    number _protein_dna_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces, number &sigma, number &b, number &rstar, number &rcut, number &stiffness);
    void check_input_sanity(std::vector<BaseParticle *> &particles) override;
    void init() override;
    virtual number _protein_repulsive_lj(const LR_vector &r, LR_vector &force, bool update_forces);
    virtual number _protein_exc_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
    virtual number _protein_ang_pot(BaseParticle *p, BaseParticle*q, bool compute_r, bool update_forces);
};

#endif /* DNANM_INTERACTION_H */
