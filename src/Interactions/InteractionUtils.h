/*
 * InteractionUtils.h
 *
 *  Created on: 11/Nov/2013
 *      Author: Flavio
 */

#ifndef INTERACTION_UTILS_H_
#define INTERACTION_UTILS_H_

#include "BaseInteraction.h"

/**
 * @brief This header contains functions that are potentially used in more than one interaction, such as box overlaps, cylinder overlaps,...
 *
 */
class InteractionUtils {
public:
	/**
 	 * @brief overlap test between two boxes (cuboids) of the same shape
 	 *
 	 * The overlap test is carried out using the separating axis theorem. It
 	 * can probably be optimized.
 	 *
 	 * @param p Particle to be tested; the orientation matrix is used to
 	 * define the box orientation
 	 * @param q Particle to be tested
 	 * @param dr Distance in between the centers of the particles. The order
 	 * is q - p (from p to q).
 	 * @param l1 box side
 	 * @param l2 box side
 	 * @param l3 box side
 	 */
	template<typename number> static bool box_overlap (BaseParticle<number> * p, BaseParticle<number> * q, LR_vector<number> dr, number l1, number l2, number l3);

	/**
	 * @brief overlap test between two spherocylinder (hard rods)
	 *
	 * This function implements the algorithm of Vega and Lago to compute
	 * the distance between two hard rods. It then tests whether this
	 * distance is smaller than 1. It might be optimized as an overlap test
	 * adding a few early exit statements that are not in there at the moment.
	 *
	 * The radius of the spherocylinders is assumed to be 0.5, and the
	 * spherocylinders are of the same length.
	 *
	 * @param dr between the centers of the two spherocylinders.
	 * @param u1 (normalised) axis identifying one of the spherocylinders.
	 * @param u2 (normalised) axis identifying the other spherocylinder.
	 * @param length length of the spherocylinders
	 */
	template<typename number> static bool spherocylinder_overlap (LR_vector<number> dr, LR_vector<number> u1, LR_vector<number> u2, number length);

	/**
	 * @brief overlap test between a disk and a sphere
	 *
	 * This function computes whether a disk and a sphere overlap.
	 * The disk is assumed to have radius 0.5, while the sphere's radius is
	 * passed as an argument.
	 *
	 * @param dr distance from disk center to sphere center.
	 * @param x x-axis on the plane identified by the disk.
	 * @param y y-axis on the plane identified by the disk.
	 * @param z axis normal to the plane identified by the disk.
	 * @param R radius of the sphere
	 */
	template<typename number> static bool disk_sphere_overlap (LR_vector<number> dr, LR_vector<number> x, LR_vector<number> y, LR_vector<number> z, number R);

	/**
	 * @brief overlap test between a cylinder and a sphere
	 *
	 * This function computes whether a cylinder and a sphere overlap.
	 * The disk is assumed to have radius 0.5, while the sphere's radius is
	 * passed as an argument.
	 *
	 * @param dr distance from cylinder center to sphere center.
	 * @param n normal to the plane on which the disk is
	 * @param height cylinder height
	 * @param R radius of the sphere
	 */
	template<typename number> static bool cylinder_sphere_overlap (LR_vector<number> dr, LR_vector<number> n, number height, number R);

	/**
	 * @brief overlap test between two cylinders
	 *
	 * The cylinder overlap test combines an overlap test in between two
	 * boxes that are built around the spherocylinder and between two
	 * spherocylinders that are build around the cylinders. If both tests
	 * return positively, the overlap is assumed. This algorithm is
	 * original, but it works fine.
	 *
	 * Different combinations might be implemented in the future as they may
	 * be more efficient.
	 *
	 * the radius is assumed to be 0.5, and the cylinders' axes point along
	 * the v3 direction of the orientation matrix of each of the
	 * BaseParticle objects feeded in.
	 *
	 * @param p first cylinder
	 * @param q second cylinder
	 * @param dr distance q-p between the two cylinders
	 * @param l lenght of the two cylinders
	 */
	template<typename number> static bool cylinder_overlap (BaseParticle<number> * p, BaseParticle<number> * q, LR_vector<number> dr, number l);

	/**
	 * @brief Overlap between a sphere and a spherocylinder
	 *
	 * @param dr distance from sphere to spherocylinder
	 * @param R_s radius of the sphere
	 * @param u versor identifying the direction of the spherocylinder
	 * @param sc_length length of the spherocylinder
	 * @param R_sc spherocylinder radius
	 */
	template<typename number> static bool sphere_spherocylinder_overlap (LR_vector<number> dr, number R_s, LR_vector<number> u, number sc_length, number R_sc);

	/**
 	 * @brief Overlap between a sphere and a box
 	 *
 	 * Overlap between a sphere and a box...
	 * idea taken from http://www.idt.mdh.se/personal/tla/publ/sb.pdf
 	 *
	 * @param dr distance from sphere to spherocylinder
	 * @param R_s radius of the sphere
 	 * @param O orientation matrix identfying the orientation of the box
 	 * @param l1 box length
 	 * @param l2 box length
 	 * @param l3 box length
 	 */
	template<typename number> static bool sphere_box_overlap (LR_vector<number> dr, number R_s, LR_matrix<number> O, number l1, number l2, number l3);
};

/// SAT to evaluate whether two boxes overlap
template<typename number>
bool InteractionUtils::box_overlap (BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> dr, number lx, number ly, number lz) {
	// here we compute the 15 potential separating axes
	LR_vector<number> sep[15];

	sep[0] = p->orientation.v1;
	sep[1] = p->orientation.v2;
	sep[2] = p->orientation.v3;

	sep[3] = q->orientation.v1;
	sep[4] = q->orientation.v2;
	sep[5] = q->orientation.v3;

	sep[6] =  p->orientation.v1.cross(q->orientation.v1);
	sep[7] =  p->orientation.v1.cross(q->orientation.v2);
	sep[8] =  p->orientation.v1.cross(q->orientation.v3);
	sep[9] =  p->orientation.v2.cross(q->orientation.v1);
	sep[10] = p->orientation.v2.cross(q->orientation.v2);
	sep[11] = p->orientation.v2.cross(q->orientation.v3);
	sep[12] = p->orientation.v3.cross(q->orientation.v1);
	sep[13] = p->orientation.v3.cross(q->orientation.v2);
	sep[14] = p->orientation.v3.cross(q->orientation.v3);

	for (int k = 6; k < 15; k ++) sep[k].normalize();

	// now we have the separating vectors; we should look for Ra and Rb
	number Ra, Rb;
	for (int k = 0; k < 15; k ++) {
		Ra = fabs(p->orientation.v1 * sep[k]) * lx / 2.+
		     fabs(p->orientation.v2 * sep[k]) * ly / 2.+
		     fabs(p->orientation.v3 * sep[k]) * lz / 2.;

		Rb = fabs(q->orientation.v1 * sep[k]) * lx / 2.+
		     fabs(q->orientation.v2 * sep[k]) * ly / 2.+
		     fabs(q->orientation.v3 * sep[k]) * lz / 2.;
		if (fabs(dr * sep[k]) > (Ra + Rb)) {
			return false;
		}
	}
	return true;
}

/// vega and Lago's function; assumes the lengths of the line segments are the same
/// and that the vectors u1 and u2 are normalized
template<typename number>
bool InteractionUtils::spherocylinder_overlap (LR_vector<number> dr, LR_vector<number> u1, LR_vector<number> u2, number length) {
	
	// this function computes a distance. Some early exit things might be implemented
	number hlength = length / 2.;
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;
	if (cc < 1.e-9 && cc > 1.e-12) {
		lambda = ( drdotu1 - u1dotu2 * drdotu2) / cc;
		mu     = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if (!(fabs (lambda) <= hlength && fabs (mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if (aux1 > aux2) {
				lambda = copysign (hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if (fabs(mu) > hlength) mu = copysign (hlength, mu);
			}
			else {
				mu = copysign (hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
		number dr1 = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;
		number lambda1=lambda;
		number mu1=mu;
		
		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if (fabs(mu) > hlength) mu = copysign(hlength, mu);
		number dr2 = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;

		//printf ("%g %g ##\n", sqrt(dr1), sqrt(dr2));
		if ((dr1 < 1.) != (dr2 < 1.)) {
			printf ("INCONSISTENT (cc=%g - %g)\n", cc, u1dotu2 * u1dotu2);
			printf ("sqrt(dr1) = %g, sqrt(dr2) = %g ##\n", sqrt(dr1), sqrt(dr2));
			printf ("lambda, mu, drdotu1, drdotu2, u1dotu2 %g %g %g %g %g\n", lambda1, mu1, drdotu1, drdotu2, u1dotu2);
			printf ("lambda, mu, drdotu1, drdotu2, u1dotu2 %g %g %g %g %g\n", lambda, mu, drdotu1, drdotu2, u1dotu2);
			printf ("u1 = np.array([% 12.9g, % 12.9g, % 12.9g])\n", u1.x, u1.y, u1.z);
			printf ("u2 = np.array([% 12.9g, % 12.9g, % 12.9g])\n", u2.x, u2.y, u2.z);
			//abort();
		}
	}

	if (cc < 1.e-12) {
		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if (fabs(mu) > hlength) mu = copysign(hlength, mu);
	}
	else {
		// line segments not parallel
		lambda = ( drdotu1 - u1dotu2 * drdotu2) / cc;
		mu     = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if (!(fabs (lambda) <= hlength && fabs (mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if (aux1 > aux2) {
				lambda = copysign (hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if (fabs(mu) > hlength) mu = copysign (hlength, mu);
			}
			else {
				mu = copysign (hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	if (dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1 <= (number) 1.f)
		return true;
	else
		return false;
}

template<typename number>
bool InteractionUtils::cylinder_sphere_overlap (LR_vector<number> dr, LR_vector<number> n, number height, number R) {
	// dr goes from cylinder to sphere; cylinder has radius 0.5
	number drpar = dr * n;
	number drpersq = dr.norm() - drpar * drpar;

	if (fabs (drpar) < (height / 2.)) {
		// sphere and cylindrical tube may overlap
		if (drpersq < (0.5 + R) * (0.5 + R)) return true;
		else return false;
	}
	else {
		// we check if the sphere is "on top"
		if (drpersq < 0.5 * 0.5) {
			if (fabs(drpar) < height / 2. + R) return true;
			else return false;
		}
		
		// here we check the disk/sphere overlap between
		// one face of the cylinder and the sphere;
		// we choose the face closest to the sphere
		LR_vector<number> mydr;
		if (drpar > 0.) mydr = dr - n * (height / 2.);
		else mydr = dr + n * (height / 2.);
		
		// early ejection
		if (mydr.norm() > (0.5 + R) * (0.5 + R)) return false;
		
		// we get two vectors perpendicular to n
		LR_vector<number> x = dr - n * drpar;
		x.normalize();
		LR_vector<number> y = n.cross(x);
		return disk_sphere_overlap (mydr, x, y, n, R);
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
		LR_vector<number> mydr;
		if (drpar > 0.) mydr = dr - n * (height / 2.);
		else mydr = dr + n * (height / 2.);
		
		// early ejection
		if (mydr.norm() > (0.5 + R) * (0.5 + R)) return false;
		
		// we get two vectors perpendicular to n
		LR_vector<number> x = dr - n * drpar;
		x.normalize();
		LR_vector<number> y = n.cross(x);
		return disk_sphere_overlap (mydr, x, y, n, R);
	}
	*/
}

template<typename number>
bool InteractionUtils::disk_sphere_overlap (LR_vector<number> dr, LR_vector<number> x, LR_vector<number> y, LR_vector<number> z, number R) {
	// assumes disk is in 0, 0, 0, with radius is 0.5;
	// dr goes from disk to sphere

	// sphere inside disk, but disk perimeter does not intersect sphere (like a ball in a glass)
	// Closest point from sphere's center to the plane identified by disk
	number drn = dr * z;
	LR_vector<number> P = dr - z * drn;
	if (fabs(drn) < R && P.norm() < 0.25) return true;

	// we now check if the perimeter intersects the sphere
	number phi = atan2 (dr * y, dr * x);

	// two possible values for phi, one corresponding to the maximum distance and one to the
	// minimum distance
	P = (x * 0.5 * cos(phi)) + (y  * 0.5 * sin(phi));
	if ((dr - P).norm() < R * R) return true;
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
template<typename number>
bool InteractionUtils::cylinder_overlap (BaseParticle<number> *p, BaseParticle<number> * q, LR_vector<number> dr, number length) {
	// one more thing that can be optimized is, instead of checking
	// (box && box) && (spherocylinder && spherocylinder)
	// we may want to check
	// (box && spherocylinder) && (box && spherocylinder)
	// depending on which is faster... box && spherocylinder has not been written yet...

	if (spherocylinder_overlap (dr, p->orientation.v3, q->orientation.v3, length)) return box_overlap (p, q, dr, (number)1.f, (number)1.f, length);
	else return false;
}*/
template<typename number>
bool InteractionUtils::cylinder_overlap (BaseParticle<number> *p, BaseParticle<number> * q, LR_vector<number> dr, number length) {

	number hlength = length / 2.;
	// dr = pos_1 - pos_2 with PBC; thus, we invert it;
	dr = -dr;
	LR_vector<number> u1 = p->orientation.v3;
	LR_vector<number> u2 = q->orientation.v3;
	u1.normalize();
	u2.normalize();
	number drdotu1 = dr * u1;
	number drdotu2 = dr * u2;
	number u1dotu2 = u1 * u2;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;
	if (cc < 1.e-6) {
		// special case of (almost) parallel line segments
		if (drdotu1 != (number) 0.f) {
			// parallel line segments, on different lines
			lambda = copysign (hlength, drdotu1);
			mu = lambda * u1dotu2 - drdotu2;
			if (fabs(mu) > hlength) mu = copysign (hlength, mu);
		}
		else {
			// parallel line segments, along the same line
			lambda = (number) 0.f;
			mu = (number) 0.f;
		}
	}
	else {
		// line segments not parallel
		lambda = ( drdotu1 - u1dotu2 * drdotu2) / cc;
		mu =     (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if (!(fabs (lambda) <= hlength && fabs (mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if (aux1 > aux2) {
				lambda = copysign (hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if (fabs(mu) > hlength) mu = copysign (hlength, mu);
			}
			else {
				mu = copysign (hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}

	// now the distance is dr + mu * u2 - lambda * u1
	LR_vector<number> dist = dr + mu * u2 - lambda * u1;
	number distnormsq = dr.norm() + lambda * lambda + mu * mu - 2.f * lambda * mu * u1dotu2 + 2.f * mu * drdotu2 - 2.f * lambda * drdotu1;
	// if no overlap between spherocylinders we eject;
	if (!(distnormsq <= (number) 1.f)) return false;

	// if the distance goes through both cylindrical rims, we can eject now
	// and report with an overlap
	if (fabs(lambda + (dist * u1)) < hlength && fabs(mu + (-dist * u2)) < hlength ) {
		//printf ("%s %d\n", __FILE__, __LINE__);
		return true;
	}

	//if (box_overlap(p, q, dr, (number)1.f, (number)1.f, length) == false) return false;

	// We now handle disk-disk overlaps.
	LR_vector<number> v = u1.cross(u2);
	v.normalize();
	for (int i = -1; i < 2; i += 2) { // face on cylinder 1
		for (int j = -1; j < 2; j += 2) { // face on cylinder 2
			// we put p in (0, 0, 0) and q in (drx, dry, drz);
			LR_vector<number> D1 = i * hlength * u1;
			LR_vector<number> D2 = dr + j * hlength * u2;

			// we preempt the check by checking if spheres overlap first
			if ((D1 - D2).norm() > (number)1.f) continue;

			// line in between planes: g = g0 * eta * v
			LR_vector<number> g0;

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
			LR_vector<number> tmpv = g0 - D1;
			number eta = - tmpv * v;
			LR_vector<number> P1 = D1 + tmpv + eta * v;

			tmpv = g0 - D2;
			eta = - tmpv * v;
			LR_vector<number> P2 = D2 + tmpv + eta * v;

			number P1D1norm = (P1 - D1).norm();
			number P2D2norm = (P2 - D2).norm();

			// no overlap in this case
			if (P1D1norm > 0.25 || P2D2norm > 0.25) continue;

			number S1 = sqrt(0.25 - P1D1norm);
			number S2 = sqrt(0.25 - P2D2norm);

			// overlap?
			if ((P1 - P2).module() < S1 + S2) {
				//printf ("%s %d\n", __FILE__, __LINE__);
				return true;
			}
		}
	}

	for (int i = -1; i < 2; i += 2) { // face on cylinder 1
		LR_vector<number> Dj = i * hlength * u1;
		LR_vector<number> Ci = dr;
		LR_vector<number> Ui = Ci + u2 * ((Dj - Ci) * u2);

		LR_vector<number> tmpv = Dj - Ui;
		number tmpvnorm = tmpv.norm();
		// sphere on disk and cylinder don't overlap
		if (tmpvnorm > (number) 1.f) continue;

		// sphere on disk is in the continuation zone of the cylinder; this case
		// would have been detected in the disk-disk overlap
 		if (tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u2) > hlength) continue;
		// center of disk within the cylinder
 		if (tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u2) < hlength) return true;

		// early check: does the cylinder axis intersect disk within the cylinder?
		// if so, we can return immediately if the intersection point is inside the
		// cylinder.
		number d = (Dj - Ci) * u1 / u1dotu2;
		if (fabs(d) < hlength) {
			LR_vector<number> Z = Ci + u2 * d;
			if ((Z - Dj).norm() < 0.5 * 0.5) {
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
		rho = sqrt(num*num + den*den);
		drho = (num * ny + den * nx) / rho;
		if (rho < 0.5) f1 = lambda + drn - rho * drho;
		else f1 = lambda + drn - 0.5 * drho;

		lambda = hlength;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num*num + den*den);
		drho = (num * ny + den * nx) / rho;
		if (rho < 0.5) f2 = lambda + drn - rho * drho;
		else f2 = lambda + drn - 0.5 * drho;

		if (f1 * f2 > 0.) {
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
			rho = sqrt(num*num + den*den);
			if (rho < 0.5) {
				func = lambda + drn - num * ny - den * nx;
				df = 1 - nx*nx - ny*ny;
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

			if (fabs(func) < 1.e-6) break;

			//printf ("%d %g %g %g %g HHH\n", cnt, lambda, func, df, ddf);
			if (cnt > 20) {
				//fprintf (stderr, "Too slow convergence. Aborting at line %d of file %s\n", __LINE__, __FILE__);
				do_bisec = true;
				break;
			}
			//if (cnt % 100 == 0 && cnt > 10) fprintf (stderr, "File %s, line %d: slow convergence, %d %g %g\n", __FILE__, __LINE__, cnt, func, df);

			if (cnt > 0 && cnt % 10 == 0) w*= 0.99;

			cnt ++;
			//lambda = lambda - w * 2. * func * df / (2 * df * df - func * ddf);
			number fact = - func / df;
			if (fabs(fact) > hlength) fact = copysign(hlength, fact);
			lambda = lambda + w * fact;
		} while (true);
		
		//do_bisec = true;
		
		// bisection method used as fallback if Newton-Rapson fails to converge
		// within 20 iterations;
		if (do_bisec) {
			int cnt2 = 0;
			
			bool stop = false;
			number lambda1 = -hlength;
			number lambda2 = +hlength;
			
			num = dry + lambda1 * ny;
			den = drx + lambda1 * nx;
			rho = sqrt (num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if (rho < 0.5) f1 = lambda1 + drn - rho * drho;
			else f1 = lambda1 + drn - 0.5 * drho;
			
			num = dry + lambda2 * ny;
			den = drx + lambda2 * nx;
			rho = sqrt (num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if (rho < 0.5) f2 = lambda2 + drn - rho * drho;
			else f2 = lambda2 + drn - 0.5 * drho;
			
			number lambdai;
			number fi;
			do {
				
				lambdai = 0.5 * (lambda1 + lambda2);
				num = dry + lambdai * ny;
				den = drx + lambdai * nx;
				rho = sqrt (num * num + den * den);
				drho = (num * ny + den * nx) / rho;
				if (rho < 0.5) fi = lambdai + drn - rho * drho;
				else fi = lambdai + drn - 0.5 * drho;

				if (fi * f1 > 0) {
					// same sign;
					lambda1 = lambdai;
					f1 = fi;
				}
				else {
					lambda2 = lambdai;
					f2 = fi;
				}

				if (fabs(fi) < 1.e-6) {
					stop = true;
				}

				cnt2 ++;

			} while (stop == false);
			func = fi;
			lambda = lambdai;
		}
		// end of fallback bisection algorithm 
		
		sinphi = num / rho;
		cosphi = den / rho;

		//LR_vector<number> T = Dj + p->orientation.v1 * 0.5 * cosphi + p->orientation.v2 * 0.5 * sinphi - (Ci + lambda * u2);
		LR_vector<number> T = Dj + p->orientation.v1 * 0.5 * cosphi + p->orientation.v2 * 0.5 * sinphi - Ci;
		number Tpar = T * u2;
		number Tpersq = T.norm() - Tpar * Tpar;
		if (fabs(Tpar) < hlength && Tpersq < 0.5 * 0.5) {
			//printf ("%s %d\n", __FILE__, __LINE__);
			return true;
		}
	}

	for (int i = -1; i < 2; i += 2) { // face on cylinder 2
		LR_vector<number> Dj = dr + u2 * (i * hlength);
		LR_vector<number> Ci (0., 0., 0.);
		LR_vector<number> Ui = Ci + u1 * ((Dj - Ci) * u1);

		LR_vector<number> tmpv = Dj - Ui;
		number tmpvnorm = tmpv.norm();
		// sphere on disk and cylinder don't overlap
		if (tmpvnorm > (number) 1.f) continue;
		// sphere on disk is in the continuation zone of the cylinder; this case
		// would have been detected in the disk-disk overlap
 		if (tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u1) > hlength) continue;

		// center of disk within the cylinder
 		if (tmpvnorm < 0.5 * 0.5 && fabs((Dj - Ci) * u1) < hlength) return true;

		// early check: does the cylinder axis intersect disk within the cylinder?
		number d = (Dj - Ci) * u2 / u1dotu2;
		if (fabs(d) < hlength) {
			LR_vector<number> Z = Ci + u1 * d;
			if ((Z - Dj).norm() < 0.5 * 0.5) {
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
		rho = sqrt(num*num + den*den);
		drho = (num * ny + den * nx) / rho;
		ddrho = (nx*nx + ny*ny - drho*drho) / rho;
		if (rho < 0.5) f1 = lambda + drn - rho * drho;
		else f1 = lambda + drn - 0.5 * drho;

		lambda = hlength;
		num = dry + lambda * ny;
		den = drx + lambda * nx;
		rho = sqrt(num*num + den*den);
		drho = (num * ny + den * nx) / rho;
		ddrho = (nx*nx + ny*ny - drho*drho) / rho;
		if (rho < 0.5) f2 = lambda + drn - rho * drho;
		else f2 = lambda + drn - 0.5 * drho;

		if (f1 * f2 > 0.) {
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
			rho = sqrt(num*num + den*den);
			if (rho < 0.5) {
				func = lambda + drn - num * ny - den * nx;
				df = 1 - nx*nx - ny*ny;
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
			if (fabs(func) < 1.e-6) break;

			if (cnt > 20) {
				// fallback to bisection algorithm
				do_bisec = true;
				break;
			}
			if (cnt % 100 == 0 && cnt > 10) fprintf (stderr, "File %s, line %d: slow convergence, %d %g %g\n", __FILE__, __LINE__, cnt, func, df);

			if (cnt > 0 && cnt % 10 == 0) w*= 0.99;

			cnt ++;
			number fact = - func / df;
			if (fabs(fact) > hlength) fact = copysign(hlength, fact);
			lambda = lambda + w * fact;
		} while (true);
		
		// bisection method used as fallback if Newton-Rapson fails to converge
		// within 20 iterations;
		//do_bisec = true;
		if (do_bisec) {
			int cnt2 = 0;
			
			bool stop = false;
			number lambda1 = -hlength;
			number lambda2 = +hlength;
			
			num = dry + lambda1 * ny;
			den = drx + lambda1 * nx;
			rho = sqrt (num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if (rho < 0.5) f1 = lambda1 + drn - rho * drho;
			else f1 = lambda1 + drn - 0.5 * drho;
			
			num = dry + lambda2 * ny;
			den = drx + lambda2 * nx;
			rho = sqrt (num * num + den * den);
			drho = (num * ny + den * nx) / rho;
			if (rho < 0.5) f2 = lambda2 + drn - rho * drho;
			else f2 = lambda2 + drn - 0.5 * drho;
			
			number lambdai;
			number fi;
			do {
				
				lambdai = 0.5 * (lambda1 + lambda2);
				num = dry + lambdai * ny;
				den = drx + lambdai * nx;
				rho = sqrt (num * num + den * den);
				drho = (num * ny + den * nx) / rho;
				if (rho < 0.5) fi = lambdai + drn - rho * drho;
				else fi = lambdai + drn - 0.5 * drho;

				if (fi * f1 > 0) {
					lambda1 = lambdai;
					f1 = fi;
				}
				else {
					lambda2 = lambdai;
					f2 = fi;
				}

				if (fabs(fi) < 1.e-6) {
					stop = true;
				}

				cnt2 ++;

			} while (stop == false);
			func = fi;
			lambda = lambdai;
		}
		// end of fallback bisection algorithm 

		sinphi = num / rho;
		cosphi = den / rho;

		LR_vector<number> T = Dj + q->orientation.v1 * (0.5 * cosphi) + q->orientation.v2 * (0.5 * sinphi) - Ci;
		number Tpar = T * u1;
		number Tpersq = T.norm() - Tpar * Tpar;
		if (fabs(Tpar) < hlength && Tpersq < 0.5 * 0.5) {
			//printf ("%g %g %g\n", Tpar, Tpersq, lambda);
			//printf ("%s %d\n", __FILE__, __LINE__);
			return true;
		}
	}

	return false;
}


template<typename number>
bool InteractionUtils::sphere_spherocylinder_overlap (LR_vector<number> dr, number sphere_radius, LR_vector<number> sc_versor, number sc_length, number sc_radius) {
	// dr goes from sphere to spherocylinder
	number drdotu = dr * sc_versor;

	//if (drdotu <  -_length / 2.) drdotu = -sc_length / 2.;
	//if (drdotu >  -_length / 2.) drdotu =  sc_length / 2.;
	// equivalent to above two lines
	//if (fabs(drdotu) >  sc_length / 2.) drdotu =  sc_length / 2.;
	number lambda = - drdotu;
	if (lambda >  sc_length / 2.) lambda =  sc_length / 2.;
	if (lambda < -sc_length / 2.) lambda = -sc_length / 2.;

	number sqr_distance = dr.norm() + lambda * lambda + 2. * lambda * drdotu;

	if (sqr_distance < (sphere_radius + sc_radius) * (sphere_radius + sc_radius)) return true;
	else return false;
}

template<typename number>
bool InteractionUtils::sphere_box_overlap (LR_vector<number> dr, number sphere_radius, LR_matrix<number> box, number l1, number l2, number l3) {
	// algorithm http://www.idt.mdh.se/personal/tla/publ/sb.pdf
	LR_vector<number> rdr (0., 0., 0.);
	dr = -dr;

	// some early exit test may be implemented here
	number ci = dr * box.v1;
	if (ci >  (l1 / 2.)) rdr += box.v1 *  (ci - (l1 / 2.));
	if (ci < -(l1 / 2.)) rdr += box.v1 * -(ci + (l1 / 2.));

	ci = (dr * box.v2);
	if (ci >  (l2 / 2.)) rdr += box.v2 *  (ci - (l2 / 2.));
	if (ci < -(l2 / 2.)) rdr += box.v2 * -(ci + (l2 / 2.));

	ci = (dr * box.v3);
	if (ci >  (l3 / 2.)) rdr += box.v3 *  (ci - (l3 / 2.));
	if (ci < -(l3 / 2.)) rdr += box.v3 * -(ci + (l3 / 2.));

	if (rdr.norm() <= sphere_radius * sphere_radius) return true;
	else return false;
}

#endif /* INTERACTION_UTILS_H_ */

