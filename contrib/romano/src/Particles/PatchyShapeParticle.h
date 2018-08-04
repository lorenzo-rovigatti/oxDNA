/*
 * PatchyShapeParticle.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef PATCHYSPARTICLE_H_
#define PATCHYSPARTICLE_H_

#include "../../../../src/Particles/BaseParticle.h"

/// A structure describing the patch; Each particle can have multiple patches, positioned at different places; The patches are directional and each
/// patch interacts only with its specific complementary patch;
template <typename number>
struct Patch {
 LR_vector<number> position; //the position of the patch with respect to the CM of the particle
 LR_vector<number> a1;  //vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
 LR_vector<number> a2; // vector which needs to be parallel with a2 on the complementary patch
 int id; //the id of the patch; it is used during initialization to assign patches to particles according to input file; sets the type of patch
 int index ; //this is the unique index of the patch in the simulation
 bool active; //is the patch on or not
 int locked_to_particle; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
 int locked_to_patch; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
 number locked_energy;

 int color; //this is the color of the patch
 number strength;  //sets the strength of the interaction
 number a1_x, a1_y, a1_z;
 number a2_x, a2_y, a2_z;


 Patch() {active = false; id = 0; color = -1; strength = 1; a1_x = a1_y = a1_z = a2_x = a2_y = a2_z = 0; set_lock(-1,-1,0);}

 Patch(LR_vector<number> _a1_xyz, LR_vector<number> _a2_xyz, LR_vector<number> _position, int _id,int _color, number _strength=1.0,  bool _active = true) :
	 position(_position), id(_id), active(_active),   color(_color), strength(_strength)
 {
	 a1_x = _a1_xyz.x;
	 a1_y = _a1_xyz.y;
	 a1_z = _a1_xyz.z;
	 a2_x = _a2_xyz.x;
	 a2_y = _a2_xyz.y;
	 a2_z = _a2_xyz.z;

	 set_lock(-1,-1,0);
 }

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
 number get_lock_energy(void) {return locked_energy;}

 void unlock() {locked_to_particle = -1;}

};

///Excluded volume center
template <typename number>
struct ExcVolCenter {
 LR_vector<number> position; //the position of the excvol with respect to the CM of the particle
 number radius;  //radius of the excvol
 //LR_vector<number> a2; // vector which needs to be parallel with a2 on the complementary patch

};

/**
 * @brief Incapsulates a patchy particle with 2, 3, or 4 spherical patches. Used by PatchyInteraction.
 */
template<typename number>
class PatchyShapeParticle : public BaseParticle<number> {
public:
	//LR_vector<number> *_base_patches;
	int N_patches; //number of patches  = number of patches
	int N_vertexes; //number of vertices of the shape; 0 = sphere

    Patch<number> *patches;
    LR_vector<number> *_vertexes;


	void _set_base_patches();

public:
	PatchyShapeParticle(int N_patches=1 , int type = 0,int N_vertexes=0);
	PatchyShapeParticle(const PatchyShapeParticle<number> &b)
	{patches = 0; this->_vertexes =  0; this->copy_from(b);}

	virtual ~PatchyShapeParticle();

	void set_positions();

	virtual void copy_from(const BaseParticle<number> &);

	PatchyShapeParticle<number>& operator = (const PatchyShapeParticle<number>& b) {this->copy_from(b);  return *this;}
	void add_patch(Patch<number> &patch,int position);

	int get_patch_color(int patchid)  {return this->patches[patchid].get_color();}
	virtual bool is_rigid_body() {
		return true;
	}

	void _set_vertexes(void);
	void _set_icosahedron_vertexes(void);

	bool locked_to_particle_id(int particle_id); // {return is_locked() && particle_id == locked_to_particle;}
	void unlock_patches(void);


};

#endif /* PLPATCHYPARTICLE_H_ */
