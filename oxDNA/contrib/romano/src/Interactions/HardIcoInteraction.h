/*
 * HardIcoInteraction.h
 *
 *  Created on: 16/Lug/2018
 *      Author: Flavio
 */

#ifndef HARDICOINTERACTION_H_
#define HARDICOINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Interactions/InteractionUtils.h"
#include "../Particles/Icosahedron.h"

/**
 * @brief Interaction class to simulate Hard Icosahedra.
 *
 * interaction_type = HardIco
 *
 * @verbatim
@endverbatim
 */
template <typename number>
class HardIcoInteraction: public BaseInteraction<number, HardIcoInteraction<number> > {
protected:
	inline number _hi_pot (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces);
	
	number _delta;
	number _tworinscribed;
	int * _close_vertexes;
	
public:
	enum {
		HardIco = 0
	};

	HardIcoInteraction();
	virtual ~HardIcoInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
};

template<typename number>
inline number HardIcoInteraction<number>::_hi_pot (BaseParticle<number> *ap, BaseParticle<number> *aq, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	/*
	for (int k = 0; k < 12; k ++) {
		LR_vector<number> tmp = *r - ap->int_centers[k];
		printf ("%3d % 8.5g % 8.5g % 8.5g %g\n", k, ap->int_centers[k].x, ap->int_centers[k].y, ap->int_centers[k].z, sqrt(tmp.norm()));
	}
	printf ("\n");
	for (int k = 0; k < 12; k ++) {
		LR_vector<number> tmp = *r + aq->int_centers[k];
		printf ("%3d % 8.5g % 8.5g % 8.5g %g\n", k, aq->int_centers[k].x, aq->int_centers[k].y, aq->int_centers[k].z, sqrt(tmp.norm()));
	}
	printf ("\n");
*/
	Icosahedron<number> * p = NULL, * q = NULL;
	LR_vector<number> my_r;
	if (ap->index < aq->index) {
		p = dynamic_cast< Icosahedron<number> *> (ap);
		q = dynamic_cast< Icosahedron<number> *> (aq);
		my_r = *r;
	}
	else {
		p = dynamic_cast< Icosahedron<number> *> (aq);
		q = dynamic_cast< Icosahedron<number> *> (ap);
		my_r = this->_box->min_image (p, q);
		//printf("switch!\n");
	}
	
	number rnorm = my_r.norm();
	if (rnorm > this->_sqr_rcut) return (number) 0.f;
	
	if (rnorm > (number) 1.) return (number) 0.f; // todo: modify for patches
	
	// radius of inscribed sphere: sqrt((1./12.) + (1./(6.*sqrt(5.)))) 
	if (rnorm < _tworinscribed * _tworinscribed) {
		this->set_is_infinite(true);
		return (number) 1.0e12;
	} 

	// now, for each dodecahedron, we must find the two closest vertexes
	/*
	number max1 = (number) 200.0f;
	number max2 = (number) 200.0f;
	int v1 = -1;
	int v2 = -1;

	for (int k = 0; k < 12; k ++) {
		//number a =  (p->int_centers[k] * (my_r));
		//number b = -(q->int_centers[k] * (my_r));
		LR_vector<number> va =   my_r + q->int_centers[k];
		LR_vector<number> vb =  -my_r + p->int_centers[k];
		number a = sqrt(va.norm());
		number b = sqrt(vb.norm());
		printf ("k, a, b: %d % 8.5f % 8.5f\n", k, a, b);
		if (a < max1) {
			max1 = a;
			v1 = k;
		}
		if (b < max2) {
			max2 = b;
			v2 = k;
		}
	}
	printf ("selecting %d and %d\n", v1, v2);
	*/
	
	/*
	number min = (number) 2. * rnorm; // this is actually min
	int v1 = -1;
	int v2 = -1;
	for (int k = 0; k < 12; k ++) {
		for (int l = 0; l < 12; l ++) {
			LR_vector<number> tmp = my_r + q->int_centers[k] - p->int_centers[l];
			if (tmp.norm() < min) {
				min = tmp.norm();
				v1 = l;
				v2 = k;
			}
		}
	}*/
	//printf ("Reselecting %d and %d\n", v1, v2);
	
	/*
	LR_vector<number> s1, s2;
	// check edges from p with faces from q
	s1 = p->int_centers[v1] - my_r;
	for (int i = 0; i < 5; i ++) {     // for each edge startging from p->int_centers[v1]
		s2 = p->int_centers[_close_vertexes[5*v1 + i]] - my_r;
		for (int j = 0; j < 5; j ++) { // for each face of q
			//printf ("checking %2d %2d %2d %2d %2d\n", v1, _close_vertexes[5*v1 + i], v2, _close_vertexes[5*v2+j], _close_vertexes[5*v2+((j+1)%5)]);
			bool check = InteractionUtils::edge_triangle_intersection(
					(s1),                                                    // starting point
					(s2),                                                    // end point
					(q->int_centers[v2]),                                // first vertex of triangle
					(q->int_centers[_close_vertexes[5 * v2 + j]]),           // second vertex of triangle
					(q->int_centers[_close_vertexes[5 * v2 + ((j + 1)%5)]])  // third vertex of triangle
				);
			if (check) {
				this->set_is_infinite(true);
				return (number)1.0e12;
			}
		}
	}
	
	s1 = q->int_centers[v2] + my_r;
	for (int i = 0; i < 5; i ++) {     // for each edge startging from q->int_centers[v2]
		s2 = q->int_centers[_close_vertexes[5*v2 + i]] + my_r;
		for (int j = 0; j < 5; j ++) { // for each face of p
			bool check = InteractionUtils::edge_triangle_intersection(
					(s1),                                                    // starting point
					(s2),                                                    // end point
					(p->int_centers[v1]),                                // first vertex of triangle
					(p->int_centers[_close_vertexes[5 * v1 + j]]),           // second vertex of triangle
					(p->int_centers[_close_vertexes[5 * v1 + ((j + 1)%5)]])  // third vertex of triangle
				);
			if (check) {
				this->set_is_infinite(true);
				return (number)1.0e12;
			}
		}
	}
	*/

	/*
	
	for (int i = 0; i < 12; i ++) {
		LR_vector<number> tmp = p->pos + p->int_centers[i];
		printf ("%+9.6f %+9.6f %+9.6f ", tmp.x, tmp.y, tmp.z);
		for (int j = 0; j < 5; j ++) 
			printf ("%2d ", _close_vertexes[5*i + j]);
		printf (" \n");
	}
	for (int i = 0; i < 12; i ++) {
		LR_vector<number> tmp = q->pos + q->int_centers[i];
		printf ("%+9.6f %+9.6f %+9.6f ", tmp.x, tmp.y, tmp.z);
		for (int j = 0; j < 5; j ++) 
			printf ("%2d ", _close_vertexes[5*i + j]);
		printf (" \n");
	}*/
	
	//throw oxDNAException ("you wish you were good, but you aren't");
	
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", p->int_centers[5][0], p->int_centers[5][1], p->int_centers[5][2]);
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", p->int_centers[3][0], p->int_centers[3][1], p->int_centers[3][2]);
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", q->int_centers[3][0], q->int_centers[3][1], q->int_centers[3][2]);
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", q->int_centers[11][0], q->int_centers[11][1], q->int_centers[11][2]);
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", q->int_centers[7][0], q->int_centers[7][1], q->int_centers[7][2]);
	
	// now double check
	int dcheck = 0;
	for (int i = 0; i < 12; i ++) {
		LR_vector<number> S1 = p->int_centers[i] - my_r;
		LR_vector<number> S2;
		LR_vector<number> P1;
		LR_vector<number> P2;
		LR_vector<number> P3;
		for (int j = 0; j < 5; j ++) {
			int other = _close_vertexes[5*i + j];
			if (other > i) {
				S2 = p->int_centers[other] - my_r;
				for (int k = 0; k < 12; k ++) {
					P1 = q->int_centers[k];
					for (int l = 0; l < 5; l ++) {
						P2 = q->int_centers[_close_vertexes[5*k+l]];
						P3 = q->int_centers[_close_vertexes[5*k+((l+1)%5)]];
						//printf("\n\nggg\n");
						//printf("%+8.5f %+8.5f %+8.5f\n", S1.x, S1.y, S1.z);
						//printf("%+8.5f %+8.5f %+8.5f\n", S2.x, S2.y, S2.z);
						//printf("%+8.5f %+8.5f %+8.5f\n", P1.x, P1.y, P1.z);
						//printf("%+8.5f %+8.5f %+8.5f\n", P2.x, P2.y, P2.z);
						//printf("%+8.5f %+8.5f %+8.5f\n", P3.x, P3.y, P3.z);
						//printf ("checking %2d %2d -- %2d %2d %2d\n", i, other, k, _close_vertexes[5*k + l], _close_vertexes[5*k + ((l+1)%5)]);
						//if (k < _close_vertexes[5*k+l] && _close_vertexes[5*k+l] < _close_vertexes[5*k+((l+1)%5)])
						//if (k < _close_vertexes[5*k+l] && _close_vertexes[5*k+l] < _close_vertexes[5*k+((l+1)%5)])
						if (InteractionUtils::edge_triangle_intersection(S1,S2,P1,P2,P3) == true) {
						//	printf ("found here p, q (line %d): i, j, k, l: %2d %2d %2d %2d -- %2d %2d %2d %2d %2d\n", __LINE__, i,j,k,l,
										//i, other, k, _close_vertexes[5*k+l], _close_vertexes[5*k+((l+1)%5)]);
							dcheck ++;
							//throw oxDNAException("prt");
							this->set_is_infinite(true);
							return (number) 1.0e12;
						}
					}
		}}}
	}
	for (int i = 0; i < 12; i ++) {
		LR_vector<number> S1 = q->int_centers[i] + my_r;
		LR_vector<number> S2;
		LR_vector<number> P1;
		LR_vector<number> P2;
		LR_vector<number> P3;
		for (int j = 0; j < 5; j ++) {
			int other = _close_vertexes[5*i + j];
			if (other > i) {
				S2 = q->int_centers[other] + my_r;
				for (int k = 0; k < 12; k ++) {
					P1 = p->int_centers[k];
					for (int l = 0; l < 5; l ++) {
						P2 = p->int_centers[_close_vertexes[5*k+l]];
						P3 = p->int_centers[_close_vertexes[5*k+((l+1)%5)]];
						//printf ("checking %2d %2d -- %2d %2d %2d\n", i, other, k, _close_vertexes[5*k + l], _close_vertexes[5*k + ((l+1)%5)]);
						//if (k < _close_vertexes[5*k+l] && _close_vertexes[5*k+l] < _close_vertexes[5*k+((l+1)%5)])
						if (InteractionUtils::edge_triangle_intersection(S1,S2,P1,P2,P3) == true) {
						//	printf ("found here q, p (line %d): i, j, k, l: %2d %2d %2d %2d -- %2d %2d %2d %2d %2d\n", __LINE__, i,j,k,l,
									//	i, other, k, _close_vertexes[5*k+l], _close_vertexes[5*k+((l+1)%5)]);
							dcheck ++;
							this->set_is_infinite(true);
							return (number) 1.0e12;
						}
					}	
		}}}
	}
	//if (dcheck != 0) throw oxDNAException ("should be zero...");
	//throw oxDNAException ("stop ");
	
	return (number) 0.;
}

extern "C" IBaseInteraction<float> * make_float() { return new HardIcoInteraction<float> (); }
extern "C" IBaseInteraction<double> * make_double() { return new HardIcoInteraction<double> (); }

#endif /* HARDICOINTERACTION_H_ */

