#include "PatchyShapeParticle.h"
#include "../Utilities/oxDNAException.h"

#define HALF_ISQRT3 0.28867513459481292f


PatchyShapeParticle::PatchyShapeParticle(int _N_patches, int _type, int _N_vertexes) :  BaseParticle() {
	this->type = _type;
	N_patches =  _N_patches;
	N_vertexes = _N_vertexes; // this is unnecessary

	// this->N_int_centers = N_patches+N_vertexes; //This is now a vector
	
	if(N_patches + N_vertexes> 0)
	{
		// this->int_centers = new LR_vector[N_patches+N_vertexes]; //Don't think one need to initilize this.
		this->int_centers.resize(N_patches+N_vertexes);
		this->patches.resize(N_patches);
		this->vertexes.resize(N_vertexes);

	}
	else
	{
		this->int_centers.clear();
		this->patches.clear();
		this->vertexes.clear();
	}
	//_set_base_patches();
}


PatchyShapeParticle::~PatchyShapeParticle() {
}


void PatchyShapeParticle::copy_from(const BaseParticle &b)
{
  const PatchyShapeParticle *bb = dynamic_cast<const PatchyShapeParticle *>(&b);
  if (bb == 0)
	  throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");

  if(!(this->int_centers.size() == b.int_centers.size() && this->N_patches == bb->N_patches && this->N_vertexes == bb->N_vertexes)){
	this->int_centers.resize(bb->int_centers.size());
	this->N_patches = bb->N_patches;
	this->N_vertexes = bb->N_vertexes;

	// this->int_centers = new LR_vector[bb->N_int_centers];
	patches.resize(bb->N_patches);
	vertexes.resize(bb->N_vertexes);
  }

//   BaseParticle::copy_from(b);

  for(int i =0 ; i < this->N_patches; i++)	this->patches[i] = bb->patches[i];
  for(int i =0 ; i < this->N_vertexes; i++)	this->vertexes[i] = bb->vertexes[i];

}

 void PatchyShapeParticle::add_patch(Patch &patch,int position) {

	if(position < 0 || position >= static_cast<int>(this->int_centers.size()))
		throw oxDNAException ("Could process patch id, please check that the patches of id %d are correct. Aborting",position);
	patches[position] = patch;
}



void PatchyShapeParticle::_set_base_patches() {
	for(uint i = 0; i < this->int_centers.size(); i++) {
		patches[i].a1.normalize();
		patches[i].a2.normalize();
	}
}


void PatchyShapeParticle::set_positions() {
	/*
	 if(this->index == 0)
        printf("I am in set_positions for particle %d, N_patches=%d, N_vertices=%d, N_int_centers=%d \n",this->index,this->N_patches,this->N_vertexes,this->N_int_centers);

    */
    //printf("Setting position with %f %f %f\n",this->orientationT.v1.x,this->orientationT.v1.y,this->orientationT.v1.z);
	for(int i = 0; i < this->N_patches; i++)
    {
		this->int_centers[i] = this->orientation * this->patches[i].position;
        this->patches[i].a1 = (patches[i].a1_x * this->orientationT.v1) + (patches[i].a1_y * this->orientationT.v2) + (patches[i].a1_z * this->orientationT.v3); //possibly can be accelerated
        this->patches[i].a2 = (patches[i].a2_x * this->orientationT.v1) + (patches[i].a2_y * this->orientationT.v2) + (patches[i].a2_z * this->orientationT.v3); //possibly can be accelerated

    }
	for(int i = 0; i < this->N_vertexes; i++)
	{
			this->int_centers[this->N_patches+i] = this->orientation * this->vertexes[i];
	}

}


void PatchyShapeParticle::unlock_patches(void) {
	for(int i = 0; i < this->N_patches; i++)
	{
		this->patches[i].set_lock(); //cleans the lock
	}
}


bool PatchyShapeParticle::locked_to_particle_id(int particle_id)
{
	for(int i = 0; i < this->N_patches; i++)
	{
		if (this->patches[i].locked_to_particle_id(particle_id) )
			 return true;
	}

	return false;
}




void PatchyShapeParticle::_set_vertexes()
{
  if(N_vertexes == 12)
  {
	  _set_icosahedron_vertexes();
  }
  else {
	  throw oxDNAException("Unsupported number of vertexes: %d\n",N_vertexes);
  }
}


void PatchyShapeParticle::_set_icosahedron_vertexes() {
	double t = (1. + sqrt(5.))/2;      // golden radius
	double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
	number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
	number b = t / (2. * r);

	if (this->N_vertexes != 12)
	{
		throw oxDNAException("If you are using icosahedron, you need 12 vertices, but we have %d",this->N_vertexes);
	}

	// 12 vertexes of the icosahedron; circular
	// permutations of (0, +/-a, +/-b)
	vertexes[ 0] = LR_vector( 0.,  a,  b);
	vertexes[ 1] = LR_vector( 0.,  a, -b);
	vertexes[ 2] = LR_vector( 0., -a,  b);
	vertexes[ 3] = LR_vector( 0., -a, -b);

	vertexes[ 4] = LR_vector(  b, 0.,  a);
	vertexes[ 5] = LR_vector(  b, 0., -a);
	vertexes[ 6] = LR_vector( -b, 0.,  a);
	vertexes[ 7] = LR_vector( -b, 0., -a);

	vertexes[ 8] = LR_vector(  a,  b, 0.);
	vertexes[ 9] = LR_vector(  a, -b, 0.);
	vertexes[10] = LR_vector( -a,  b, 0.);
	vertexes[11] = LR_vector( -a, -b, 0.);

	// we now need to figure out which vertexes belong to which face
	// angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
	// the possible angles between all pairs are: ~63.4, ~116.6 e 180.
	// each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)


	number thres = 0.;  // threshold angle is 90 degrees
	for (int i = 0; i < 12; i ++) {
		for (int j = 0; j < i; j ++) {
			for (int k = 0; k < j; k ++) {
				if ((vertexes[i]*vertexes[j] > thres) &&
				    (vertexes[i]*vertexes[k] > thres) &&
				    (vertexes[j]*vertexes[k] > thres)) {
				}
			}
		}
	}





	this->set_positions();
}

// template class PatchyShapeParticle<float>;
// template class PatchyShapeParticle<double>;
