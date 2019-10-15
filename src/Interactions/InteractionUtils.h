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
	static bool box_overlap(BaseParticle * p, BaseParticle * q, LR_vector dr, number l1, number l2, number l3);

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
	static bool spherocylinder_overlap(LR_vector dr, LR_vector u1, LR_vector u2, number length);

	/**
	 * @brief distance between two spherocylinders.
	 *
	 * All the options are the same as the above function
	 */
	static number spherocylinder_distance(LR_vector dr, LR_vector u1, LR_vector u2, number length);
	static LR_vector spherocylinder_vector(LR_vector dr, LR_vector u1, LR_vector u2, number length);

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
	static bool disk_sphere_overlap(LR_vector dr, LR_vector x, LR_vector y, LR_vector z, number R);

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
	static bool cylinder_sphere_overlap(LR_vector dr, LR_vector n, number height, number R);

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
	static bool cylinder_overlap(BaseParticle * p, BaseParticle * q, LR_vector dr, number l);

	/**
	 * @brief Overlap between a sphere and a spherocylinder
	 *
	 * @param dr distance from sphere to spherocylinder
	 * @param R_s radius of the sphere
	 * @param u versor identifying the direction of the spherocylinder
	 * @param sc_length length of the spherocylinder
	 * @param R_sc spherocylinder radius
	 */
	static bool sphere_spherocylinder_overlap(LR_vector dr, number R_s, LR_vector u, number sc_length, number R_sc);

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
	static bool sphere_box_overlap(LR_vector dr, number R_s, LR_matrix O, number l1, number l2, number l3);

	/*
	 * @brief Solves a 3 eq in 3 unknowns system
	 *
	 * Function to solve a system of 3 equations in 3 unknowns.
	 *
	 * @param C matrix defining the system
	 * @param rhs right-hand side of the system
	 * @param sol pointer to an LR_vector where the solution will be stored
	 */
	// TODO: perhaps make it return the vector instead of using the pointer?
	static LR_vector gauss_elim_solve(LR_matrix &C, LR_vector &rhs);

	/*
	 * @brief checks wheter a segment intersects a triangle in 3D
	 *
	 * This function checks if a segment intersects a triangle in 3D.
	 *
	 * @param S1 pointer to vector defining the start of the segment
	 * @param S2 pointer to vector defining the end of the segment
	 * @param P1 pointer to vector defining the first vertex of the triangle
	 * @param P2 pointer to vector defining the second vertex of the triangle
	 * @param P3 pointer to vector defining the third vertex of the triangle
	 */
	static bool edge_triangle_intersection(LR_vector &S1, LR_vector &S2, LR_vector &P1, LR_vector &P2, LR_vector &P3);
};

#endif /* INTERACTION_UTILS_H_ */
