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
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
};

template<typename number>
inline number HardIcoInteraction<number>::_hi_pot (BaseParticle<number> *ap, BaseParticle<number> *aq, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");
	
	Icosahedron<number> * p = NULL, * q = NULL;
	LR_vector<number> my_r;
	if (ap->index < aq->index) {
		p = static_cast< Icosahedron<number> *> (ap);
		q = static_cast< Icosahedron<number> *> (aq);
		my_r = *r;
	}
	else {
		p = static_cast< Icosahedron<number> *> (aq);
		q = static_cast< Icosahedron<number> *> (ap);
		my_r = this->_box->min_image (p, q);
		//printf("switch!\n");
	}
	
	number rnorm = my_r.norm();
	if (rnorm > this->_sqr_rcut) return (number) 0.f;
	
	if (rnorm > (number) 1.0f) return (number) 0.f; // todo: modify for patches
	
	// radius of inscribed sphere: sqrt((1./12.) + (1./(6.*sqrt(5.)))) 
	if (rnorm < _tworinscribed * _tworinscribed) {
		this->set_is_infinite(true);
		return (number) 1.0e12;
	} 

	// now, for each dodecahedron, we must find the two closest vertexes
	number max1 = (number) -2.0f;
	number max2 = (number) -2.0f;
	number max11 = -3.;
	number max12 = -4.;
	number max21 = -3.;
	number max22 = -4.;
	int v1 = -1;
	int v2 = -1;
	int v11 = -1;
	int v12 = -1;
	int v21 = -1;
	int v22 = -1;
	
	// finds the three closest vertexes. Can possibly be optimized
	for (int k = 0; k < 6; k ++) {
		number a =  (p->int_centers[k] * (my_r));
		number b = -(q->int_centers[k] * (my_r));
		if (a > max1) {
			max12 = max11;
			max11 = max1;
			max1 = a;
			v12 = v11;
			v11 = v1;
			v1 = k;
		} else if (a > max11) {
			max12 = max11;
			max11 = a;
			v12 = v11;
			v11 = k;
		} else if (a > max12) {
			max12 = a;
			v12 = k;
		}
		//printf ("-        %2d % 8.5f % 8.5f % 8.5f %2d %2d %2d\n", k, max1, max11, max12, v1, v11, v12);
		if (b > max2) {
			max22 = max21;
			max21 = max2;
			max2 = b;
			v22 = v21;
			v21 = v2;
			v2 = k;
		} else if (b > max21) {
			max22 = max21;
			max21 = b;
			v22 = v21;
			v21 = k;
		} else if (b > max22) {
			max22 = b;
			v22 = k;
		}
		//printf ("-        %2d % 8.5f % 8.5f % 8.5f %2d %2d %2d\n", k, max2, max21, max22, v2, v21, v22);
		
		int kk = k + 6;
		a = -a;
		b = -b;
		if (a > max1) {
			max12 = max11;
			max11 = max1;
			max1 = a;
			v12 = v11;
			v11 = v1;
			v1 = kk;
		} else if (a > max11) {
			max12 = max11;
			max11 = a;
			v12 = v11;
			v11 = kk;
		} else if (a > max12) {
			max12 = a;
			v12 = kk;
		}
		//printf ("-        %2d % 8.5f % 8.5f % 8.5f %2d %2d %2d\n", k, max1, max11, max12, v1, v11, v12);
		if (b > max2) {
			max22 = max21;
			max21 = max2;
			max2 = b;
			v22 = v21;
			v21 = v2;
			v2 = kk;
		} else if (b > max21) {
			max22 = max21;
			max21 = b;
			v22 = v21;
			v21 = kk;
		} else if (b > max22) {
			max22 = b;
			v22 = kk;
		}
	}
	
	// Early exit with no overlap; this is the Separating Axis Theorem, using
	// the distance as the candedate axis. If we find light, it means
	// there cannot be an ovelrap
	//number dd = rmod - max2/(0.5*rmod) - max1/(0.5*rmod);
	number dd = rnorm - max2/0.5 - max1/0.5; // we avoid a sqrt()
	if (dd > (number) 0.)
		return 0.;
	
	// we try the two closest faces
	LR_vector<number> pf = ((number)(1./3.))*(p->int_centers[v1] + p->int_centers[v11] + p->int_centers[v12]);
	
	bool would_exit = false;
	number rin = 0.5 * _tworinscribed;
	
	LR_vector<number> d = pf / rin; // now normalized
	number qpf = (my_r + q->int_centers[v2] - pf) * d;
	number tmp = (my_r + q->int_centers[v21] - pf) * d;
	if (tmp < qpf) qpf = tmp;
	tmp = (my_r + q->int_centers[v22] - pf) * d;
	if (tmp < qpf) qpf = tmp;
	// we use a safe threshold of 0.015 to account for non-orthonormal rotation matrixes
	if (qpf > (number) 0.015) would_exit = true;
	if (would_exit) return (number) 0.;
	
	LR_vector<number> qf = ((number)(1./3.))*(q->int_centers[v2] + q->int_centers[v21] + q->int_centers[v22]);
	d = qf / rin;
	number ppf = (p->int_centers[v1] - my_r - qf) * d;
	tmp = (p->int_centers[v11] - my_r - qf) * d;
	if (tmp < ppf) ppf = tmp;
	tmp = (p->int_centers[v12] - my_r - qf) * d;
	if (tmp < ppf) ppf = tmp;
	// we use a safe threshold of 0.015 to account for non-orthonormal rotation matrixes
	if (ppf > (number) 0.015) would_exit = true;
	if (would_exit) return (number) 0.;
	

	// check edges from p with faces from q
	LR_vector<number> s1, s2;
	s1 = p->int_centers[v1] - my_r;
	for (int i = 0; i < 5; i ++) {     // for each edge startging from p->int_centers[v1]
		s2 = p->int_centers[_close_vertexes[5*v1 + i]] - my_r;
		for (int j = 0; j < 5; j ++) { // for each face of q
			bool check = InteractionUtils::edge_triangle_intersection(
					(s1),                                                    // starting point
					(s2),                                                    // end point
					(q->int_centers[v2]),                                    // first vertex of triangle
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
					(p->int_centers[v1]),	                                 // first vertex of triangle
					(p->int_centers[_close_vertexes[5 * v1 + j]]),           // second vertex of triangle
					(p->int_centers[_close_vertexes[5 * v1 + ((j + 1)%5)]])  // third vertex of triangle
				);
			if (check) {
				this->set_is_infinite(true);
				return (number)1.0e12;
			}
		}
	}
	
	return (number) 0.;
	
	// below, there is a thorough check that goes through all faces etc. to make
	// sure that everything is ok
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
	}
	
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", p->int_centers[5][0], p->int_centers[5][1], p->int_centers[5][2]);
	//printf ("%+8.5f %+8.5f %+8.5f @@\n", q->int_centers[3][0], q->int_centers[3][1], q->int_centers[3][2]);
	
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
						if (InteractionUtils::edge_triangle_intersection(S1,S2,P1,P2,P3) == true) {
						//	printf ("found here p, q (line %d): i, j, k, l: %2d %2d %2d %2d -- %2d %2d %2d %2d %2d\n", __LINE__, i,j,k,l,
										//i, other, k, _close_vertexes[5*k+l], _close_vertexes[5*k+((l+1)%5)]);
							dcheck ++;
							//throw oxDNAException("prt");
							//this->set_is_infinite(true);
							//return (number) 1.0e12;
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
						if (InteractionUtils::edge_triangle_intersection(S1,S2,P1,P2,P3) == true) {
						//	printf ("found here q, p (line %d): i, j, k, l: %2d %2d %2d %2d -- %2d %2d %2d %2d %2d\n", __LINE__, i,j,k,l,
									//	i, other, k, _close_vertexes[5*k+l], _close_vertexes[5*k+((l+1)%5)]);
							dcheck ++;
							//this->set_is_infinite(true);
							//return (number) 1.0e12;
						}
					}	
		}}}
	}
	if (dcheck != 0) throw oxDNAException ("should be zero %d %d", p->index, q->index);
	//throw oxDNAException ("stop ");
	
	return (number) 0.;
	*/
}

extern "C" BaseInteraction<float> * make_float() { return new HardIcoInteraction<float> (); }
extern "C" BaseInteraction<double> * make_double() { return new HardIcoInteraction<double> (); }

#endif /* HARDICOINTERACTION_H_ */

