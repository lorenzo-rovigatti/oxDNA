// Subhajit 

#include "PHBParticle.h"

PHBParticle::PHBParticle():CCGParticle(){
    this->int_centers.resize(patches.size()+vertexes.size());
}

PHBParticle::PHBParticle(int totPatches,int totVertexes){
    patches.resize(totPatches);
    vertexes.resize(totVertexes);
}

PHBParticle::~PHBParticle(){
}

void PHBParticle::copy_from(PHBParticle p){
    patches.resize(p.patches.size());
    vertexes.resize(p.vertexes.size());
    for(uint i=0;i<patches.size();i++) patches[i]=p.patches[i];
    for(uint i=0;i<vertexes.size();i++) vertexes[i]=p.vertexes[i];
}

void PHBParticle::add_patch(Patch &patch) {
	patches.push_back(patch);
	int_centers.resize(patches.size()+vertexes.size());
}

// void PHBParticle::set_base_patches() {
// 	for(uint i = 0; i < this->int_centers.size(); i++) {
// 		patches[i].a1.normalize();
// 		patches[i].a2.normalize();
// 	}
// }

void PHBParticle::set_positions(){
	#pragma omp parallel for
    for(uint i=0;i<patches.size();i++){
        int_centers[i]=orientation*patches[i].position;
        // patches[i].a1=patches[i].a1static.x*orientationT.v1+patches[i].a1static.y*orientationT.v2+patches[i].a1static.z*orientationT.v3;
        // patches[i].a2=patches[i].a2static.x*orientationT.v1+patches[i].a2static.y*orientationT.v2+patches[i].a2static.z*orientationT.v3;
    }
    // for(uint i=0;i<vertexes.size();i++) int_centers[patches.size()+i] = orientation*vertexes[i];
}
void PHBParticle::add_patch(Patch patch,int position) {

	if(position < 0 || position >= static_cast<int>(this->int_centers.size()))
		throw oxDNAException ("Could process patch id, please check that the patches of id %d are correct. Aborting",position);
	patches[position] = patch;
}

bool PHBParticle::locked_to_particle_id(int particle_id)
{
	for(uint i = 0; i < patches.size(); i++)
		if (this->patches[i].locked_to_particle_id(particle_id) ) 
            return true;
            
	return false;
}

/// old reminents

// void PHBParticle::set_vertexes()
// {
//   if(vertexes.size() == 12)
//   {
// 	  set_icosahedron_vertexes();
//   }
//   else {
// 	  throw oxDNAException("Unsupported number of vertexes: %d\n",vertexes.size());
//   }
// }

// void PHBParticle::set_icosahedron_vertexes() {
// 	double t = (1. + sqrt(5.))/2;      // golden radius
// 	double r = 2. * sin(2. * M_PI/5);  // radius of icosahedron of edge lenght two
// 	number a = 1./ (2. * r);           // will generate patches at a distance 0.5 from center
// 	number b = t / (2. * r);

// 	if (vertexes.size() != 12)
// 	{
// 		throw oxDNAException("If you are using icosahedron, you need 12 vertices, but we have %d",this->vertexes.size());
// 	}

// 	// 12 vertexes of the icosahedron; circular
// 	// permutations of (0, +/-a, +/-b)
// 	vertexes[ 0] = LR_vector( 0.,  a,  b);
// 	vertexes[ 1] = LR_vector( 0.,  a, -b);
// 	vertexes[ 2] = LR_vector( 0., -a,  b);
// 	vertexes[ 3] = LR_vector( 0., -a, -b);

// 	vertexes[ 4] = LR_vector(  b, 0.,  a);
// 	vertexes[ 5] = LR_vector(  b, 0., -a);
// 	vertexes[ 6] = LR_vector( -b, 0.,  a);
// 	vertexes[ 7] = LR_vector( -b, 0., -a);

// 	vertexes[ 8] = LR_vector(  a,  b, 0.);
// 	vertexes[ 9] = LR_vector(  a, -b, 0.);
// 	vertexes[10] = LR_vector( -a,  b, 0.);
// 	vertexes[11] = LR_vector( -a, -b, 0.);

// 	// we now need to figure out which vertexes belong to which face
// 	// angle between vertexes of the same face: 1.1071487177940436, ~63.435 degrees
// 	// the possible angles between all pairs are: ~63.4, ~116.6 e 180.
// 	// each vertex has 5 vertexes at 63, 5 more at 116 and 1 at 180 (opposite)


// 	number thres = 0.;  // threshold angle is 90 degrees
// 	for (int i = 0; i < 12; i ++) {
// 		for (int j = 0; j < i; j ++) {
// 			for (int k = 0; k < j; k ++) {
// 				if ((vertexes[i]*vertexes[j] > thres) &&
// 				    (vertexes[i]*vertexes[k] > thres) &&
// 				    (vertexes[j]*vertexes[k] > thres)) {
// 				}
// 			}
// 		}
// 	}





// 	this->set_positions();
// }