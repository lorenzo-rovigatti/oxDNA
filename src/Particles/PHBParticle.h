// Subhajit

//btype = color and 100 is colorless

#ifndef PHBParticle_H
#define PHBParticle_H

#include "CCGParticle.h"
// #include "PatchyShapeParticle.h" //tried this but failed pathetically

struct Patch {
    LR_vector position; //the position of the patch with respect to the CM of the particle
    LR_vector a1;  //vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
    LR_vector a2; // vector which needs to be parallel with a2 on the complementary patch
    LR_vector a1static={1,0,0},a2static={1,0,0};
    int id=0; //the id of the patch; it is used during initialization to assign patches to particles according to input file; sets the type of patch
    int index ; //this is the unique index of the patch in the simulation
    bool active = true; //is the patch on or not
    int locked_to_particle=-1; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
    int locked_to_patch=-1; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
    number locked_energy=0;

    int color=-1; //this is the color of the patch
    number strength=1;  //sets the strength of the interaction
    // number a1_x=1, a1_y=0, a1_z=0;
    // number a2_x=1, a2_y=0, a2_z=0;


    // Patch() {};

    // Patch(LR_vector _a1_xyz, LR_vector _a2_xyz, LR_vector _position, int _id,int _color, number _strength=1.0,  bool _active = true) :
    //     position(_position), id(_id), active(_active),   color(_color), strength(_strength),a1static(_a1_xyz),a2static(_a2_xyz){};

    void set_lock(int particle=-1,int patch=-1,number energy=0)
    {
        locked_to_particle = particle;
        locked_to_patch = patch;
        locked_energy = energy;
    }

    bool is_locked(void) {return locked_to_particle >= 0;}

    int get_color(void) {return color;}

    bool locked_to(int particle_id,int patch_id) {return is_locked() && (particle_id == locked_to_particle && locked_to_patch == patch_id);}
    bool locked_to_particle_id(int particle_id)  {return is_locked() && particle_id == locked_to_particle;}
    void get_lock(int& particle_id, int& patch_id) {particle_id = locked_to_particle; patch_id =  locked_to_patch;}
    number get_lock_energy() {return locked_energy;}
    void unlock() {locked_to_particle = -1;}

};

class PHBParticle: public CCGParticle {
public:

    //strange variables
    number th_b_0=0,beta_b_0=0,xu_bending=0.952319757,xk_bending=1.14813301,kb1=0,kb2=0.80,kt_pref=1;
    PHBParticle();
    PHBParticle(int totPatches, int totVertexes=0);
    virtual ~PHBParticle();

    //Patchy settings
    std::vector<Patch> patches; //store patches
    std::vector<LR_vector> vertexes; //store vertexes
    virtual void copy_from(PHBParticle p);
    void add_patch(Patch &patch);
    void set_base_patches();
    void set_positions();
    void add_patch(Patch patch,int position);

    int get_patch_color(int patchid)  {return this->patches[patchid].get_color();}
    // PHBParticle& operator = (const PHBParticle& b) {this->copy_from(b);  return *this;}
    virtual bool is_rigid_body() {
		return true;
	}

    void set_vertexes(void);
	void set_icosahedron_vertexes(void);

	bool locked_to_particle_id(int particle_id); // {return is_locked() && particle_id == locked_to_particle;}
	void unlock_patches(void);
};

#endif