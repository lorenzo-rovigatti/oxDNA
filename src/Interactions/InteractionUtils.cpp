/*
 * InteractionUtils.cpp
 *
 *  Created on: 15 ott 2019
 *      Author: lorenzo
 */

#include "InteractionUtils.h"

/// SAT to evaluate whether two boxes overlap
bool InteractionUtils::box_overlap(BaseParticle *p, BaseParticle *q, LR_vector dr, number lx, number ly, number lz) {
	// here we compute the 15 potential separating axes
	LR_vector sep[15];

	sep[0] = p->orientation.v1;
	sep[1] = p->orientation.v2;
	sep[2] = p->orientation.v3;

	sep[3] = q->orientation.v1;
	sep[4] = q->orientation.v2;
	sep[5] = q->orientation.v3;

	sep[6] = p->orientation.v1.cross(q->orientation.v1);
	sep[7] = p->orientation.v1.cross(q->orientation.v2);
	sep[8] = p->orientation.v1.cross(q->orientation.v3);
	sep[9] = p->orientation.v2.cross(q->orientation.v1);
	sep[10] = p->orientation.v2.cross(q->orientation.v2);
	sep[11] = p->orientation.v2.cross(q->orientation.v3);
	sep[12] = p->orientation.v3.cross(q->orientation.v1);
	sep[13] = p->orientation.v3.cross(q->orientation.v2);
	sep[14] = p->orientation.v3.cross(q->orientation.v3);

	for(int k = 6; k < 15; k++)
		sep[k].normalize();

	// now we have the separating vectors; we should look for Ra and Rb
	number Ra, Rb;
	for(int k = 0; k < 15; k++) {
		Ra = fabs(p->orientation.v1 * sep[k]) * lx / 2. + fabs(p->orientation.v2 * sep[k]) * ly / 2. + fabs(p->orientation.v3 * sep[k]) * lz / 2.;

		Rb = fabs(q->orientation.v1 * sep[k]) * lx / 2. + fabs(q->orientation.v2 * sep[k]) * ly / 2. + fabs(q->orientation.v3 * sep[k]) * lz / 2.;
		if(fabs(dr * sep[k]) > (Ra + Rb)) {
			return false;
		}
	}
	return true;
}

/// vega and Lago's function; assumes the lengths of the line segments are the same
/// and that the vectors u1 and u2 are normalized

bool InteractionUtils::spherocylinder_overlap(LR_vector dr, LR_vector u1, LR_vector u2, number length) {

	// this function computes a distance. Some early exit things might be implemented
	number hlength = length / 2.;
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;

	// the following if gets almost never executed, it's just here for checks
	if(cc < 1.e-9 && cc > 1.e-12) {
		lambda = (drdotu1 - u1dotu2 * drdotu2) / cc;
		mu = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if(!(fabs(lambda) <= hlength && fabs(mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if(aux1 > aux2) {
				lambda = copysign(hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if(fabs(mu) > hlength) mu = copysign(hlength, mu);
			}
			else {
				mu = copysign(hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
		number dr1 = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;
		number lambda1 = lambda;
		number mu1 = mu;

		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if(fabs(mu) > hlength) mu = copysign(hlength, mu);
		number dr2 = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;

		//printf ("%g %g ##\n", sqrt(dr1), sqrt(dr2));
		if((dr1 < 1.) != (dr2 < 1.)) {
			printf("INCONSISTENT (cc=%g - %g)\n", cc, u1dotu2 * u1dotu2);
			printf("sqrt(dr1) = %g, sqrt(dr2) = %g ##\n", sqrt(dr1), sqrt(dr2));
			printf("lambda, mu, drdotu1, drdotu2, u1dotu2 %g %g %g %g %g\n", lambda1, mu1, drdotu1, drdotu2, u1dotu2);
			printf("lambda, mu, drdotu1, drdotu2, u1dotu2 %g %g %g %g %g\n", lambda, mu, drdotu1, drdotu2, u1dotu2);
			printf("u1 = np.array([% 12.9g, % 12.9g, % 12.9g])\n", u1.x, u1.y, u1.z);
			printf("u2 = np.array([% 12.9g, % 12.9g, % 12.9g])\n", u2.x, u2.y, u2.z);
			//abort();
		}
	}

	if(cc < 1.e-12) {
		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if(fabs(mu) > hlength) mu = copysign(hlength, mu);
	}
	else {
		// line segments not parallel
		lambda = (drdotu1 - u1dotu2 * drdotu2) / cc;
		mu = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if(!(fabs(lambda) <= hlength && fabs(mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if(aux1 > aux2) {
				lambda = copysign(hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if(fabs(mu) > hlength) mu = copysign(hlength, mu);
			}
			else {
				mu = copysign(hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	if(dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1 <= (number) 1.f) return true;
	else return false;
}

// distance between two spherocylinders

LR_vector InteractionUtils::spherocylinder_vector(LR_vector dr, LR_vector u1, LR_vector u2, number length) {

	// this function computes a distance. Some early exit things might be implemented
	number hlength = length / 2.;
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;

	if(cc < 1.e-12) {
		// parallel line segments
		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if(fabs(mu) > hlength) mu = copysign(hlength, mu);
	}
	else {
		// line segments not parallel
		lambda = (drdotu1 - u1dotu2 * drdotu2) / cc;
		mu = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if(!(fabs(lambda) <= hlength && fabs(mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if(aux1 > aux2) {
				lambda = copysign(hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if(fabs(mu) > hlength) mu = copysign(hlength, mu);
			}
			else {
				mu = copysign(hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	LR_vector res = dr - lambda * u1 + mu * u2;

	return res;
}

// distance between two spherocylinders

number InteractionUtils::spherocylinder_distance(LR_vector dr, LR_vector u1, LR_vector u2, number length) {
	LR_vector r = InteractionUtils::spherocylinder_vector(dr, u1, u2, length);
	return sqrt(r.norm());
}

bool InteractionUtils::cylinder_sphere_overlap(LR_vector dr, LR_vector n, number height, number R) {
	// dr goes from cylinder to sphere; cylinder has radius 0.5
	number drpar = dr * n;
	number drpersq = dr.norm() - drpar * drpar;

	if(fabs(drpar) < (height / 2.)) {
		// sphere and cylindrical tube may overlap
		if(drpersq < (0.5 + R) * (0.5 + R)) return true;
		else return false;
	}
	else {
		// we check if the sphere is "on top"
		if(drpersq < 0.5 * 0.5) {
			if(fabs(drpar) < height / 2. + R) return true;
			else return false;
		}

		// here we check the disk/sphere overlap between
		// one face of the cylinder and the sphere;
		// we choose the face closest to the sphere
		LR_vector mydr;
		if(drpar > 0.) mydr = dr - n * (height / 2.);
		else mydr = dr + n * (height / 2.);

		// early ejection
		if(mydr.norm() > (0.5 + R) * (0.5 + R)) return false;

		// we get two vectors perpendicular to n
		LR_vector x = dr - n * drpar;
		x.normalize();
		LR_vector y = n.cross(x);
		return disk_sphere_overlap(mydr, x, y, n, R);
	}

	/* // different version of the above, where the calculation
	 * // of drpersq is within the if, and the check for the sphere
	 * // position "on top" of the cylinder is removed. The performaces
	 * // seem basically unaffected
	 if (fabs (drpar) < (height / 2.)) {
	 // sphere and cylindrical tube may overlap
	 number drpersq = dr.norm() - drpar * drpar;
	 if (drpersq < (0.5 + R) * (0.5 + R)) return true;
	 else return false;
	 }
	 else {
	 // here we check the disk/sphere overlap between
	 // one face of the cylinder and the sphere;
	 // we choose the face closest to the sphere
	 LR_vector mydr;
	 if (drpar > 0.) mydr = dr - n * (height / 2.);
	 else mydr = dr + n * (height / 2.);

	 // early ejection
	 if (mydr.norm() > (0.5 + R) * (0.5 + R)) return false;

	 // we get two vectors perpendicular to n
	 LR_vector x = dr - n * drpar;
	 x.normalize();
	 LR_vector y = n.cross(x);
	 return disk_sphere_overlap (mydr, x, y, n, R);
	 }
	 */
}

bool InteractionUtils::disk_sphere_overlap(LR_vector dr, LR_vector x, LR_vector y, LR_vector z, number R) {
	// assumes disk is in 0, 0, 0, with radius is 0.5;
	// dr goes from disk to sphere

	// sphere inside disk, but disk perimeter does not intersect sphere (like a ball in a glass)
	// Closest point from sphere's center to the plane identified by disk
	number drn = dr * z;
	LR_vector P = dr - z * drn;
	if(fabs(drn) < R && P.norm() < 0.25) return true;

	// we now check if the perimeter intersects the sphere
	number phi = atan2(dr * y, dr * x);

	// two possible values for phi, one corresponding to the maximum distance and one to the
	// minimum distance
	P = (x * 0.5 * cos(phi)) + (y * 0.5 * sin(phi));
	if((dr - P).norm() < R * R) return true;
	//P = (x * cos(phi + M_PI)) + (y  * sin(phi + M_PI));
	//if ((dr - P).norm() < R*R) return true;
	//equivalent to above two commented lines:
	//if ((dr + P).norm() < R * R) printf ("HHHH\n");
	//if ((dr + P).norm() < R * R) return true;
	// (the above condition is redundant, it is taken care by atan2)

	return false;
}

/// overlap between cylinders; combines box and spherocylinder overlap
/*

 bool InteractionUtils::cylinder_overlap (BaseParticle *p, BaseParticle * q, LR_vector dr, number length) {
 // one more thing that can be optimized is, instead of checking
 // (box && box) && (spherocylinder && spherocylinder)
 // we may want to check
 // (box && spherocylinder) && (box && spherocylinder)
 // depending on which is faster... box && spherocylinder has not been written yet...

 if (spherocylinder_overlap (dr, p->orientation.v3, q->orientation.v3, length)) return box_overlap (p, q, dr, (number)1.f, (number)1.f, length);
 else return false;
 }*/

bool InteractionUtils::cylinder_overlap(BaseParticle *p, BaseParticle * q, LR_vector dr, number length) {

	number hlength = length / 2.;
	// dr = pos_1 - pos_2 with PBC; thus, we invert it;
	// dr = -dr; NOT NEEDED SINCE LAST MAJOR UPDATE
	LR_vector u1 = p->orientation.v3;
	LR_vector u2 = q->orientation.v3;
	u1.normalize();
	u2.normalize();
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;
	if(cc < 1.e-6) {
		// special case of (almost) parallel line segments
		if(drdotu1 != (number) 0.f) {
			// parallel line segments, on different lines
			lambda = copysign(hlength, drdotu1);
			mu = lambda * u1dotu2 - drdotu2;
			if(fabs(mu) > hlength) mu = copysign(hlength, mu);
		}
		else {
			// parallel line segments, along the same line
			lambda = (number) 0.f;
			mu = (number) 0.f;
		}
	}
	else {
		// line segments not parallel
		lambda = (drdotu1 - u1dotu2 * drdotu2) / cc;
		mu = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if(!(fabs(lambda) <= hlength && fabs(mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if(aux1 > aux2) {
				lambda = copysign(hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if(fabs(mu) > hlength) mu = copysign(hlength, mu);
			}
			else {
				mu = copysign(hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if(fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	// now the distance is dr + mu * u2 - lambda * u1
	LR_vector dist = dr + mu * u2 - lambda * u1;
	number distnormsq = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;
	// if no overlap between spherocylinders we eject;
	if(!(distnormsq <= (number) 1.f)) return false;

	// if the distance goes through both cylindrical rims, we can eject now
	// and report with an overlap
	if(fabs(lambda + (dist * u1)) < hlength && fabs(mu + (-dist * u2)) < hlength) {
		//printf ("%s %d\n", __FILE__, __LINE__);
		return true;
	}

	//if (box_overlap(p, q, dr, (number)1.f, (number)1.f, length) == false) return false;

	// We now handle disk-disk overlaps.
	LR_vector v = u1.cross(u2);
	v.normalize();
	for(int i = -1; i < 2; i += 2) { // face on cylinder 1
		for(int j = -1; j < 2; j += 2) { // face on cylinder 2
			// we put p in (0, 0, 0) and q in (drx, dry, drz);
			LR_vector D1 = i * hlength * u1;
			LR_vector D2 = dr + j * hlength * u2;

			// we preempt the check by checking if spheres overlap first
			if((D1 - D2).norm() > (number) 1.f) continue;

			// line in between planes: g = g0 * eta * v
			LR_vector g0;

			number den = u1.y * u2.z - u2.y * u1.z;
			//if (fabs(den) > 1.e9 || fabs(den) < 1.e-9) fprintf (stderr, "WARNING: large den %g\n", den);
			number D1n1 = D1 * u1;
			number D2n2 = D2 * u2;
			g0.x = 0; // one of the coordinates is arbitrary
			g0.y = (D1n1 * u2.z - D2n2 * u1.z) / den;
			g0.z = (D2n2 * u1.y - D1n1 * u2.y) / den;

			// now, the line g = g0 + eta * v identifies the intersection
			// P1 is the point that belongs to g and is perpendicular to
			// the center of the first disk
			LR_vector tmpv = g0 - D1;
			number eta = -tmpv * v;
			LR_vector P1 = D1 + tmpv + eta * v;

			tmpv = g0 - D2;
			eta = -tmpv * v;
			LR_vector P2 = D2 + tmpv + eta * v;

			number P1D1norm = (P1 - D1).norm();
			number P2D2norm = (P2 - D2).norm();

			// no overlap in this case
			if(P1D1norm > 0.25 || P2D2norm > 0.25) continue;

			number S1 = sqrt(0.25 - P1D1norm);
			number S2 = sqrt(0.25 - P2D2norm);

			// overlap?
			if((P1 - P2).module() < S1 + S2) {
				//printf ("%s %d\n", __FILE__, __LINE__);
				return true;
			}
		}
	}

	for(int i = -1; i < 2; i += 2) { // face on cylinder 1
		LR_vector Dj = i * hlength * u1;
		LR_vector Ci = dr;
		LR_vector Ui = Ci + u2 * ((Dj - Ci) * u2);

		LR_vector tmpv = Dj - Ui;
		number tmpvnorm = tmpv.norm();
		// sphere on disk and cylinder don't overlap
		if(tmpvnorm > (number) 1.f) continue;

		// sphere on disk is in the continuation zone of the cylinder; this case
		// would have been detected in the disk-disk overlap
		if(tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u2) > hlength) continue;
		// center of disk within the cylinder
		if(tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u2) < hlength) return true;

		// early check: does the cylinder axis intersect disk within the cylinder?
		// if so, we can return immediately if the intersection point is inside the
		// cylinder.
		number d = (Dj - Ci) * u1 / u1dotu2;
		if(fabs(d) < hlength) {
			LR_vector Z = Ci + u2 * d;
			if((Z - Dj).norm() < 0.5 * 0.5) {
				return true;
			}
		}

		// we find the root of lambda + drn - sin(phi) - cos(phi)
		// with lambda = lamnda(phi) = -drn + R cos(phi) nx + R sin(phi) ny;
		number drx = (Ci - Dj) * p->orientation.v1;
		number dry = (Ci - Dj) * p->orientation.v2;
		number drn = (Ci - Dj) * u2;
		number nx = u2 * p->orientation.v1;
		number ny = u2 * p->orientation.v2;

		number cosphi, sinphi, func, df;
		number num, den;

		// first, we check wether we have the same sign:
		lambda = -hlength;
		number rho, drho, ddrho, f1, f2;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num * num + den * den);
		drho = (num * ny + den * nx) / rho;
		if(rho < 0.5) f1 = lambda + drn - rho * drho;
		else f1 = lambda + drn - 0.5 * drho;

		lambda = hlength;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num * num + den * den);
		drho = (num * ny + den * nx) / rho;
		if(rho < 0.5) f2 = lambda + drn - rho * drho;
		else f2 = lambda + drn - 0.5 * drho;

		if(f1 * f2 > 0.) {
			// f1 and f2 have the same sign... thus, the minimum distance
			// between axis and disk is not contained in the cylinder. Thus
			// again, if this were an overlap, it would mean that the disk
			// overlaps with one of the faces of the cylinder, which has been
			// checked already. Thus even more, we can safely continue;
			continue;
		}

		bool do_bisec = false;
		lambda = 0.;
		number w = 0.95;
		int cnt = 0;
		do {
			num = dry + lambda * ny;
			den = drx + lambda * nx;
			rho = sqrt(num * num + den * den);
			if(rho < 0.5) {
				func = lambda + drn - num * ny - den * nx;
				df = 1 - nx * nx - ny * ny;
				//ddf = 0.;
			}
			else {
				drho = (num * ny + den * nx) / rho;
				ddrho = (nx * nx + ny * ny - drho * drho) / rho;
				//dddrho = -3. * drho * ddrho / rho;
				func = lambda + drn - 0.5 * drho;
				df = 1. - 0.5 * ddrho;
				//ddf = - 0.5 * dddrho;
			}

			if(fabs(func) < 1.e-6) break;

			//printf ("%d %g %g %g %g HHH\n", cnt, lambda, func, df, ddf);
			if(cnt > 20) {
				//fprintf (stderr, "Too slow convergence. Aborting at line %d of file %s\n", __LINE__, __FILE__);
				do_bisec = true;
				break;
			}
			//if (cnt % 100 == 0 && cnt > 10) fprintf (stderr, "File %s, line %d: slow convergence, %d %g %g\n", __FILE__, __LINE__, cnt, func, df);

			if(cnt > 0 && cnt % 10 == 0) w *= 0.99;

			cnt++;
			//lambda = lambda - w * 2. * func * df / (2 * df * df - func * ddf);
			number fact = -func / df;
			if(fabs(fact) > hlength) fact = copysign(hlength, fact);
			lambda = lambda + w * fact;
		} while(true);

		//do_bisec = true;

		// bisection method used as fallback if Newton-Rapson fails to converge
		// within 20 iterations;
		if(do_bisec) {
			int cnt2 = 0;

			bool stop = false;
			number lambda1 = -hlength;
			number lambda2 = +hlength;

			num = dry + lambda1 * ny;
			den = drx + lambda1 * nx;
			rho = sqrt(num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if(rho < 0.5) f1 = lambda1 + drn - rho * drho;
			else f1 = lambda1 + drn - 0.5 * drho;

			num = dry + lambda2 * ny;
			den = drx + lambda2 * nx;
			rho = sqrt(num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if(rho < 0.5) f2 = lambda2 + drn - rho * drho;
			else f2 = lambda2 + drn - 0.5 * drho;

			number lambdai;
			number fi;
			do {

				lambdai = 0.5 * (lambda1 + lambda2);
				num = dry + lambdai * ny;
				den = drx + lambdai * nx;
				rho = sqrt(num * num + den * den);
				drho = (num * ny + den * nx) / rho;
				if(rho < 0.5) fi = lambdai + drn - rho * drho;
				else fi = lambdai + drn - 0.5 * drho;

				if(fi * f1 > 0) {
					// same sign;
					lambda1 = lambdai;
					f1 = fi;
				}
				else {
					lambda2 = lambdai;
					f2 = fi;
				}

				if(fabs(fi) < 1.e-6) {
					stop = true;
				}

				cnt2++;

			} while(stop == false);
			func = fi;
			lambda = lambdai;
		}
		// end of fallback bisection algorithm

		sinphi = num / rho;
		cosphi = den / rho;

		//LR_vector T = Dj + p->orientation.v1 * 0.5 * cosphi + p->orientation.v2 * 0.5 * sinphi - (Ci + lambda * u2);
		LR_vector T = Dj + p->orientation.v1 * 0.5 * cosphi + p->orientation.v2 * 0.5 * sinphi - Ci;
		number Tpar = T * u2;
		number Tpersq = T.norm() - Tpar * Tpar;
		if(fabs(Tpar) < hlength && Tpersq < 0.5 * 0.5) {
			//printf ("%s %d\n", __FILE__, __LINE__);
			return true;
		}
	}

	for(int i = -1; i < 2; i += 2) { // face on cylinder 2
		LR_vector Dj = dr + u2 * (i * hlength);
		LR_vector Ci(0., 0., 0.);
		LR_vector Ui = Ci + u1 * ((Dj - Ci) * u1);

		LR_vector tmpv = Dj - Ui;
		number tmpvnorm = tmpv.norm();
		// sphere on disk and cylinder don't overlap
		if(tmpvnorm > (number) 1.f) continue;
		// sphere on disk is in the continuation zone of the cylinder; this case
		// would have been detected in the disk-disk overlap
		if(tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u1) > hlength) continue;

		// center of disk within the cylinder
		if(tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u1) < hlength) return true;

		// early check: does the cylinder axis intersect disk within the cylinder?
		number d = (Dj - Ci) * u2 / u1dotu2;
		if(fabs(d) < hlength) {
			LR_vector Z = Ci + u1 * d;
			if((Z - Dj).norm() < 0.5 * 0.5) {
				//printf ("%s %d\n", __FILE__, __LINE__);
				return true;
			}
		}

		// we minimize sin(phi) * (drx + lambda xn) - cos(phi) (dry + lambda xn)
		// with lambda = lamnda(phi) = -drn + R cos(phi) nx + R sin(phi) ny;
		number drx = (Ci - Dj) * q->orientation.v1;
		number dry = (Ci - Dj) * q->orientation.v2;
		number drn = (Ci - Dj) * u1;
		number nx = u1 * q->orientation.v1;
		number ny = u1 * q->orientation.v2;

		number cosphi, sinphi, func, df;
		number num, den;

		// first, we check wether we have the same sign:
		lambda = -hlength;
		number rho, drho, ddrho, f1, f2;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num * num + den * den);
		drho = (num * ny + den * nx) / rho;
		ddrho = (nx * nx + ny * ny - drho * drho) / rho;
		if(rho < 0.5) f1 = lambda + drn - rho * drho;
		else f1 = lambda + drn - 0.5 * drho;

		lambda = hlength;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num * num + den * den);
		drho = (num * ny + den * nx) / rho;
		ddrho = (nx * nx + ny * ny - drho * drho) / rho;
		if(rho < 0.5) f2 = lambda + drn - rho * drho;
		else f2 = lambda + drn - 0.5 * drho;

		if(f1 * f2 > 0.) {
			// f1 and f2 have the same sign... thus, the minimum distance
			// between axis and disk is not contained in the cylinder. Thus
			// again, if this were an overlap, it would mean that the disk
			// overlaps with one of the faces of the cylinder, which has been
			// checked already. Thus even more, we can safely continue;
			continue;
		}

		bool do_bisec = false;
		lambda = 0.;
		number w = 0.95;
		int cnt = 0;
		do {
			num = dry + lambda * ny;
			den = drx + lambda * nx;
			rho = sqrt(num * num + den * den);
			if(rho < 0.5) {
				func = lambda + drn - num * ny - den * nx;
				df = 1 - nx * nx - ny * ny;
				//ddf = 0.;
			}
			else {
				drho = (num * ny + den * nx) / rho;
				ddrho = (nx * nx + ny * ny - drho * drho) / rho;
				//dddrho = -3. * drho * ddrho / rho;
				func = lambda + drn - 0.5 * drho;
				df = 1. - 0.5 * ddrho;
				//ddf = - 0.5 * dddrho;
			}

			//printf ("%d %g %g %g %g HHH 2\n", cnt, lambda, func, df, ddf);
			if(fabs(func) < 1.e-6) break;

			if(cnt > 20) {
				// fallback to bisection algorithm
				do_bisec = true;
				break;
			}
			if(cnt % 100 == 0 && cnt > 10) fprintf(stderr, "File %s, line %d: slow convergence, %d %g %g\n", __FILE__, __LINE__, cnt, func, df);

			if(cnt > 0 && cnt % 10 == 0) w *= 0.99;

			cnt++;
			number fact = -func / df;
			if(fabs(fact) > hlength) fact = copysign(hlength, fact);
			lambda = lambda + w * fact;
		} while(true);

		// bisection method used as fallback if Newton-Rapson fails to converge
		// within 20 iterations;
		//do_bisec = true;
		if(do_bisec) {
			int cnt2 = 0;

			bool stop = false;
			number lambda1 = -hlength;
			number lambda2 = +hlength;

			num = dry + lambda1 * ny;
			den = drx + lambda1 * nx;
			rho = sqrt(num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if(rho < 0.5) f1 = lambda1 + drn - rho * drho;
			else f1 = lambda1 + drn - 0.5 * drho;

			num = dry + lambda2 * ny;
			den = drx + lambda2 * nx;
			rho = sqrt(num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if(rho < 0.5) f2 = lambda2 + drn - rho * drho;
			else f2 = lambda2 + drn - 0.5 * drho;

			number lambdai;
			number fi;
			do {

				lambdai = 0.5 * (lambda1 + lambda2);
				num = dry + lambdai * ny;
				den = drx + lambdai * nx;
				rho = sqrt(num * num + den * den);
				drho = (num * ny + den * nx) / rho;
				if(rho < 0.5) fi = lambdai + drn - rho * drho;
				else fi = lambdai + drn - 0.5 * drho;

				if(fi * f1 > 0) {
					lambda1 = lambdai;
					f1 = fi;
				}
				else {
					lambda2 = lambdai;
					f2 = fi;
				}

				if(fabs(fi) < 1.e-6) {
					stop = true;
				}

				cnt2++;

			} while(stop == false);
			func = fi;
			lambda = lambdai;
		}
		// end of fallback bisection algorithm

		sinphi = num / rho;
		cosphi = den / rho;

		LR_vector T = Dj + q->orientation.v1 * (0.5 * cosphi) + q->orientation.v2 * (0.5 * sinphi) - Ci;
		number Tpar = T * u1;
		number Tpersq = T.norm() - Tpar * Tpar;
		if(fabs(Tpar) < hlength && Tpersq < 0.5 * 0.5) {
			//printf ("%g %g %g\n", Tpar, Tpersq, lambda);
			//printf ("%s %d\n", __FILE__, __LINE__);
			return true;
		}
	}

	return false;
}

bool InteractionUtils::sphere_spherocylinder_overlap(LR_vector dr, number sphere_radius, LR_vector sc_versor, number sc_length, number sc_radius) {
	// dr goes from sphere to spherocylinder
	number drdotu = dr * sc_versor;

	//if (drdotu <  -_length / 2.) drdotu = -sc_length / 2.;
	//if (drdotu >  -_length / 2.) drdotu =  sc_length / 2.;
	// equivalent to above two lines
	//if (fabs(drdotu) >  sc_length / 2.) drdotu =  sc_length / 2.;
	number lambda = -drdotu;
	if(lambda > sc_length / 2.) lambda = sc_length / 2.;
	if(lambda < -sc_length / 2.) lambda = -sc_length / 2.;

	number sqr_distance = dr.norm() + lambda * lambda + 2. * lambda * drdotu;

	if(sqr_distance < (sphere_radius + sc_radius) * (sphere_radius + sc_radius)) return true;
	else return false;
}

bool InteractionUtils::sphere_box_overlap(LR_vector dr, number sphere_radius, LR_matrix box, number l1, number l2, number l3) {
	// algorithm http://www.idt.mdh.se/personal/tla/publ/sb.pdf
	LR_vector rdr(0., 0., 0.);
	dr = -dr;

	// some early exit test may be implemented here
	number ci = dr * box.v1;
	if(ci > (l1 / 2.)) rdr += box.v1 * (ci - (l1 / 2.));
	if(ci < -(l1 / 2.)) rdr += box.v1 * -(ci + (l1 / 2.));

	ci = (dr * box.v2);
	if(ci > (l2 / 2.)) rdr += box.v2 * (ci - (l2 / 2.));
	if(ci < -(l2 / 2.)) rdr += box.v2 * -(ci + (l2 / 2.));

	ci = (dr * box.v3);
	if(ci > (l3 / 2.)) rdr += box.v3 * (ci - (l3 / 2.));
	if(ci < -(l3 / 2.)) rdr += box.v3 * -(ci + (l3 / 2.));

	if(rdr.norm() <= sphere_radius * sphere_radius) return true;
	else return false;
}

// helper function to solve systems by Gaussian elimination
// TODO: perhaps check that passing by reference is faster
// TODO: perhaps remove the copying of Ca to C
LR_vector InteractionUtils::gauss_elim_solve(LR_matrix &C, LR_vector &rhs) {
	int keys[3] = { -1, -1, -1 };

	// find maximum of the first column
	int maxid = 0;
	if(fabs(C[0][0]) > fabs(C[1][0])) maxid = (fabs(C[0][0]) > fabs(C[2][0]) ? 0 : 2);
	else maxid = ((fabs(C[1][0]) > fabs(C[2][0])) ? 1 : 2);
	keys[0] = maxid;

	int i = (maxid + 1) % 3;
	int j = (maxid + 2) % 3;
	keys[1] = (fabs(C[i][1]) > fabs(C[j][1])) ? i : j;

	keys[2] = 3 - keys[1] - keys[0];

	LR_vector u(C[keys[0]][0], C[keys[0]][1], C[keys[0]][2]);
	LR_vector v(C[keys[1]][0], C[keys[1]][1], C[keys[1]][2]);
	LR_vector w(C[keys[2]][0], C[keys[2]][1], C[keys[2]][2]);
	LR_vector r(rhs[keys[0]], rhs[keys[1]], rhs[keys[2]]);

	number f = (u[0] / v[0]);
	v = u - f * v;
	r[1] = r[0] - f * r[1];
	f = (u[0] / w[0]);
	w = u - f * w;
	r[2] = r[0] - f * r[2];

	f = (v[1] / w[1]);
	w = v - f * w;
	r[2] = r[1] - f * r[2];

	number z = r[2] / w[2];
	number y = (r[1] - v[2] * z) / v[1];
	number x = (r[0] - u[1] * y - u[2] * z) / u[0];

	return LR_vector(x, y, z);
}

// intersection test between a segment and a triangle,
// defined by its three vertexes
bool InteractionUtils::edge_triangle_intersection(LR_vector &S1, LR_vector &S2, LR_vector &P1, LR_vector &P2, LR_vector &P3) {
	LR_vector rhs;
	LR_vector res;

	LR_matrix C(P2.x - P1.x, P3.x - P1.x, S1.x - S2.x, P2.y - P1.y, P3.y - P1.y, S1.y - S2.y, P2.z - P1.z, P3.z - P1.z, S1.z - S2.z);

	rhs = S1 - P1;
	res = InteractionUtils::gauss_elim_solve(C, rhs);
	if(res[0] >= (number) 0. && res[1] >= (number) 0. && res[0] + res[1] <= (number) 1. && res[2] >= (number) 0. && res[2] <= (number) 1.) {
		return true;
	}

	return false;
}

// test for the intersection between two triangles, the first
// has vertexes v1,v2,v3 and the second has vertexes u1,u2,u3.
// TODO: perhaps implement faster test presented
// in the paper http://webee.technion.ac.il/~ayellet/Ps/TroppTalShimshoni.pdf
// DOI 10.1002/cav.115
// TODO: not tested!!!!
bool triangle_intersection(LR_vector * P1, LR_vector * P2, LR_vector * P3, LR_vector * Q1, LR_vector * Q2, LR_vector * Q3) {
	// for each vertex of face 1, we need to see if it intersects face 2
	// this will happen if, by solving the equation
	// P + a1 * e1 + a2 * e2 = Qi + bi * ei

	// faces and edges
	// TODO: memory leaks?
	LR_vector *P = new LR_vector[3];
	LR_vector *Q = new LR_vector[3];
	LR_vector *p = new LR_vector[3];
	LR_vector *q = new LR_vector[3];

	P[0] = *P1;
	P[1] = *P2;
	P[2] = *P3;
	Q[0] = *Q1;
	Q[1] = *Q2;
	Q[2] = *Q3;

	p[0] = *P2 - *P1;
	p[1] = *P3 - *P2;
	p[2] = *P1 - *P3;
	q[0] = *Q2 - *Q1;
	q[1] = *Q3 - *Q2;
	q[2] = *Q1 - *Q3;

	LR_matrix C;
	LR_vector rhs;
	LR_vector res; // solution: will contain a1, a2 and bi in this order

	// check if edges of triangle Q go through triangle P
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			C = LR_matrix(p[i].x, -p[(i + 2) % 3].x, q[j].x, p[i].y, -p[(i + 2) % 3].y, q[j].y, p[i].z, -p[(i + 2) % 3].z, q[j].z);
			rhs = Q[j] - P[i];
			res = InteractionUtils::gauss_elim_solve(C, rhs);
			if(res[0] >= (number) 0. && res[1] >= (number) 0. && res[0] + res[1] <= 1. && res[2] >= (number) 0. && res[2] <= (number) 1.) return true;
		}
	}

	// check if edges of triangle P go through triangle Q
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			C = LR_matrix(q[i].x, -q[(i + 2) % 3].x, p[j].x, q[i].y, -q[(i + 2) % 3].y, p[j].y, q[i].z, -q[(i + 2) % 3].z, p[j].z);
			rhs = P[j] - Q[i];
			res = InteractionUtils::gauss_elim_solve(C, rhs);
			if(res[0] >= (number) 0. && res[1] >= (number) 0. && res[0] + res[1] <= 1. && res[2] >= (number) 0. && res[2] <= (number) 1.) return true;
		}
	}

	// the next line has been added to get rid of a warning
	return false;
}

