/*
 * PatchyShapeParticle.cpp
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#include "PatchyShapeParticle.h"
#include "../../../../src/Utilities/oxDNAException.h"

#define HALF_ISQRT3 0.28867513459481292f

PatchyShapeParticle::PatchyShapeParticle(int _N_patches, int _type, int _N_vertexes) :  BaseParticle() {
	this->type = _type;
	N_patches =  _N_patches;
	N_vertexes = _N_vertexes;

	if(N_patches + N_vertexes > 0)
	{
		this->int_centers = std::vector<LR_vector>(N_patches+N_vertexes);
		this->patches = new Patch[N_patches];
		this->_vertexes = new LR_vector[N_vertexes];

	}
	else
	{
        // constructor to std::vector<LR_vector> of length 0 is implicit
		this->patches = 0;
		this->_vertexes = 0;
	}
	//_set_base_patches();
}

PatchyShapeParticle::~PatchyShapeParticle() {
	delete[] patches;
	delete [] _vertexes;
}

void PatchyShapeParticle::copy_from(const BaseParticle &b)
{
    const PatchyShapeParticle *bb = dynamic_cast<const PatchyShapeParticle *>(&b);
    if (bb == 0)
    {
        throw oxDNAException("Can't convert particle to PatchyShapeParticle by dynamic cast'. Aborting");
    }

    delete [] this->patches;
    delete [] this->_vertexes;

    this->N_patches = bb->N_patches;
    this->N_vertexes = bb->N_vertexes;

    this->int_centers = std::vector<LR_vector>(bb->int_centers.size());
    patches = new Patch[bb->N_patches];
    _vertexes = new LR_vector[bb->N_vertexes];

    BaseParticle::copy_from(b);

    for(int i =0 ; i < this->N_patches; i++)
    {
        // this->int_centers[i] = bb->int_centers[i];
        this->patches[i] = bb->patches[i];
    }
    for(int i =0 ; i < this->N_vertexes; i++)
    {
        // this->int_centers[i] = bb->int_centers[i];
        this->_vertexes[i] = bb->_vertexes[i];
    }
}

void PatchyShapeParticle::add_patch(Patch &patch,int position) {

	if(position < 0 || position >= this->N_int_centers())
	{
		 throw oxDNAException ("Could process patch id, please check that the patches of id %d are correct. Aborting",position);
	}
	patches[position] = patch;
}


void PatchyShapeParticle::_set_base_patches() {
	for(int i = 0; i < this->N_int_centers(); i++) {
		patches[i].a1.normalize();
		patches[i].a2.normalize();
		//patches[i] *= 0.5;
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
			this->int_centers[this->N_patches+i] = this->orientation * this->_vertexes[i];
	}

	/*
	if(this->index == 0)
	{
		for(int i = 0; i < this->N_int_centers; i++)
			{
			printf("%d center: %f %f %f vertex: %f %f %f \n",i,this->int_centers[i].x,this->int_centers[i].y,this->int_centers[i].z,this->_vertexes[i].x,this->_vertexes[i].y,this->_vertexes[i].z);
			}
	}
	*/

}

void PatchyShapeParticle::unlock_patches() {
	for(int i = 0; i < this->N_patches; i++)
	{
		this->patches[i].set_lock(); //cleans the lock
	}
}

void PatchyShapeParticle::set_lock(int patch, int other_particle, int other_patch, number energy) {
    assert(patch < N_patches);
    assert(-1 < patch);
    patches[patch].set_lock(other_particle, other_patch, energy);
}

void PatchyShapeParticle::unlock(int patch) {
    assert(patch < N_patches);
    assert(-1 < patch);
    patches[patch].unlock();
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
	_vertexes[ 0] = LR_vector( 0.,  a,  b);
	_vertexes[ 1] = LR_vector( 0.,  a, -b);
	_vertexes[ 2] = LR_vector( 0., -a,  b);
	_vertexes[ 3] = LR_vector( 0., -a, -b);

	_vertexes[ 4] = LR_vector(  b, 0.,  a);
	_vertexes[ 5] = LR_vector(  b, 0., -a);
	_vertexes[ 6] = LR_vector( -b, 0.,  a);
	_vertexes[ 7] = LR_vector( -b, 0., -a);

	_vertexes[ 8] = LR_vector(  a,  b, 0.);
	_vertexes[ 9] = LR_vector(  a, -b, 0.);
	_vertexes[10] = LR_vector( -a,  b, 0.);
	_vertexes[11] = LR_vector( -a, -b, 0.);

	// we now need to figure out which vertexes belong to which face
	// angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
	// the possible angles between all pairs are: ~63.4, ~116.6 e 180.
	// each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)


	int nface = 0;
	number thres = 0.;  // threshold angle is 90 degrees
	for (int i = 0; i < 12; i ++) {
		for (int j = 0; j < i; j ++) {
			for (int k = 0; k < j; k ++) {
				if ((_vertexes[i]*_vertexes[j] > thres) &&
					(_vertexes[i]*_vertexes[k] > thres) &&
					(_vertexes[j]*_vertexes[k] > thres)) {
						//_faces[3*nface + 0] = i;
						//_faces[3*nface + 1] = j;
						//_faces[3*nface + 2] = k;
						/*printf ("\n\n%d %d %d @\n", i, j, k);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[i].x, _vertexes[i].y, _vertexes[i].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[j].x, _vertexes[j].y, _vertexes[j].z);
						printf ("%7.5g %7.5g %7.5g\n", _vertexes[k].x, _vertexes[k].y, _vertexes[k].z);
						printf ("  %g %g %g\n", 4.f*(_vertexes[i]*_vertexes[j]), 4.f*(_vertexes[i]*_vertexes[k]), 4.*(_vertexes[j]*_vertexes[k]));*/
						nface ++;
				}
			}
		}
	}

	this->set_positions();
}