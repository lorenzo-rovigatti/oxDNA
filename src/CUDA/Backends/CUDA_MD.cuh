/* Verlet related constants */
__constant__ float MD_sqr_verlet_skin[1];
/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_dt[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

// this function modifies L
// Conversion of the R matrix in _get_updated_orientation above into a quaternion simplifies into the qR quaternion used in _get_updated_orientation. This is then multiplied by the original orientation quaternion to produce the updated quaternion. Quaternion multiplication is the cross_product - dot_product and is therefore non-commutative
__device__ GPU_quat _get_updated_orientation(c_number4 &L, GPU_quat &old_o) {
	c_number norm = _module(L);
	L.x /= norm;
	L.y /= norm;
	L.z /= norm;

	c_number sintheta, costheta;
	sincos(MD_dt[0] * norm, &sintheta, &costheta);
	c_number qw = (c_number) 0.5f * sqrtf(fmaxf((c_number) 0.f, (c_number) 2.f + 2.f*costheta));
	c_number winv = (c_number)1.f /qw;
	GPU_quat R = {(c_number) 0.5f*L.x*sintheta*winv, (c_number) 0.5f*L.y*sintheta*winv, (c_number) 0.5f*L.z*sintheta*winv, qw};

	return quat_multiply(old_o, R);
}

__global__ void first_step(c_number4 *poss, GPU_quat *orientations, c_number4 *list_poss, c_number4 *vels, c_number4 *Ls, c_number4 *forces, c_number4 *torques, bool *are_lists_old) {
	if(IND >= MD_N[0]) return;

	const c_number4 F = forces[IND];

	c_number4 r = poss[IND];
	c_number4 v = vels[IND];

	v.x += F.x * (MD_dt[0] * (c_number) 0.5f);
	v.y += F.y * (MD_dt[0] * (c_number) 0.5f);
	v.z += F.z * (MD_dt[0] * (c_number) 0.5f);

	r.x += v.x * MD_dt[0];
	r.y += v.y * MD_dt[0];
	r.z += v.z * MD_dt[0];

	vels[IND] = v;
	poss[IND] = r;

	const c_number4 T = torques[IND];
	c_number4 L = Ls[IND];

	L.x += T.x * (MD_dt[0] * (c_number) 0.5f);
	L.y += T.y * (MD_dt[0] * (c_number) 0.5f);
	L.z += T.z * (MD_dt[0] * (c_number) 0.5f);

	Ls[IND] = L;


	GPU_quat qold_o = orientations[IND];
	orientations[IND] = _get_updated_orientation(L, qold_o);

	// do verlet lists need to be updated?
	if(quad_distance(r, list_poss[IND]) > MD_sqr_verlet_skin[0]) are_lists_old[0] = true;
}

__global__ void compute_molecular_coms(c_number4 *mol_coms, int *particles_to_mols, int *mol_sizes, c_number4 *poss) {
	if(IND >= MD_N[0]) {
		return;
	}

	int mol_id = particles_to_mols[IND];
	c_number4 p_contrib = poss[IND] / (c_number) mol_sizes[mol_id];

	LR_atomicAddXYZ(&(mol_coms[mol_id]), p_contrib);
}

__global__ void rescale_molecular_positions(c_number4 *mol_coms, int *particles_to_mols, c_number4 *poss, c_number4 shift_factor) {
	if(IND >= MD_N[0]) {
		return;
	}

	c_number4 ppos = poss[IND];
	c_number4 mol_com = mol_coms[particles_to_mols[IND]];
	ppos.x += mol_com.x * shift_factor.x;
	ppos.y += mol_com.y * shift_factor.y;
	ppos.z += mol_com.z * shift_factor.z;
	poss[IND] = ppos;
}

__global__ void rescale_positions(c_number4 *poss, c_number4 ratio) {
	if(IND >= MD_N[0]) {
		return;
	}
	c_number4 ppos = poss[IND];
	ppos.x *= ratio.x;
	ppos.y *= ratio.y;
	ppos.z *= ratio.z;
	poss[IND] = ppos;
}

__global__ void set_external_forces(c_number4 *poss, GPU_quat *orientations, CUDA_trap *ext_forces, c_number4 *forces, c_number4 *torques, llint step, int max_ext_forces, CUDABox *box) {
	if(IND >= MD_N[0]) return;
	// if there are no external forces then just put the force to 0
	if(max_ext_forces == 0) {
		forces[IND] = make_c_number4(0, 0, 0, 0);
		torques[IND] = make_c_number4(0, 0, 0, 0);
		return;
	}

	c_number4 ppos = poss[IND];
	c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 T = make_c_number4(0, 0, 0, 0);

	for(int i = 0; i < max_ext_forces; i++) {
		// coalesced
		CUDA_trap extF = ext_forces[MD_N[0] * i + IND];
		switch (extF.type) {
			case CUDA_TRAP_CONSTANT: {
				c_number intensity = extF.constant.F0 + step*extF.constant.rate;
				float3 dir = make_float3(extF.constant.x, extF.constant.y, extF.constant.z);
				if(extF.constant.dir_as_centre) {
					dir.x -= ppos.x;
					dir.y -= ppos.y;
					dir.z -= ppos.z;
					c_number dir_mod = _module(dir);
					dir.x /= dir_mod;
					dir.y /= dir_mod;
					dir.z /= dir_mod;
				}
				F.x += dir.x*intensity;
				F.y += dir.y*intensity;
				F.z += dir.z*intensity;
				break;
			}
			case CUDA_TRAP_MUTUAL: {
				c_number4 qpos = poss[extF.mutual.p_ind];

				c_number4 dr = (extF.mutual.PBC) ? box->minimum_image(ppos, qpos) : qpos - ppos;
				c_number dr_abs = _module(dr);

				c_number4 force = dr * ((dr_abs - (extF.mutual.r0 + extF.mutual.rate * step)) * (extF.mutual.stiff + (step * extF.mutual.stiff_rate)) / dr_abs);

				F.x += force.x;
				F.y += force.y;
				F.z += force.z;
				break;
			}
			case CUDA_TRAP_MOVING: {
				c_number tx = extF.moving.pos0.x + extF.moving.rate * step * extF.moving.dir.x;
				c_number ty = extF.moving.pos0.y + extF.moving.rate * step * extF.moving.dir.y;
				c_number tz = extF.moving.pos0.z + extF.moving.rate * step * extF.moving.dir.z;

				F.x += - extF.moving.stiff * (ppos.x - tx);
				F.y += - extF.moving.stiff * (ppos.y - ty);
				F.z += - extF.moving.stiff * (ppos.z - tz);
				break;
			}
			case CUDA_TRAP_MOVING_LOWDIM: {
				c_number tx = extF.lowdim.pos0.x + extF.lowdim.rate * step * extF.lowdim.dir.x;
				c_number ty = extF.lowdim.pos0.y + extF.lowdim.rate * step * extF.lowdim.dir.y;
				c_number tz = extF.lowdim.pos0.z + extF.lowdim.rate * step * extF.lowdim.dir.z;

				if(extF.lowdim.visX) F.x += - extF.lowdim.stiff * (ppos.x - tx);
				if(extF.lowdim.visY) F.y += - extF.lowdim.stiff * (ppos.y - ty);
				if(extF.lowdim.visZ) F.z += - extF.lowdim.stiff * (ppos.z - tz);
				break;
			}
			case CUDA_REPULSION_PLANE: {
				c_number distance = extF.repulsionplane.dir.x*ppos.x + extF.repulsionplane.dir.y*ppos.y + extF.repulsionplane.dir.z*ppos.z + extF.repulsionplane.position;
				if(distance < 0.0f) {
					F.x += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.x;
					F.y += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.y;
					F.z += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.z;
				}
				break;
			}
			case CUDA_REPULSION_PLANE_MOVING: {
				for(int idx = extF.repulsionplanemoving.low_idx; idx <= extF.repulsionplanemoving.high_idx; idx++) {
					c_number4 qpos = poss[idx];
					c_number distance = extF.repulsionplanemoving.dir.x*(ppos.x-qpos.x) + extF.repulsionplanemoving.dir.y*(ppos.y-qpos.y) + extF.repulsionplanemoving.dir.z*(ppos.z-qpos.z);
					if(distance < 0.0f) {
						F.x += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.x;
						F.y += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.y;
						F.z += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.z;
					}
				}
				break;
			}
			case CUDA_REPULSIVE_SPHERE: {
				c_number4 centre = make_c_number4(extF.repulsivesphere.centre.x, extF.repulsivesphere.centre.y, extF.repulsivesphere.centre.z, 0.);
				c_number4 dist = box->minimum_image(centre, ppos);
				c_number mdist = _module(dist);
				c_number radius = extF.repulsivesphere.r0 + extF.repulsivesphere.rate*(c_number) step;
				c_number radius_ext = extF.repulsivesphere.r_ext;
				if(mdist > radius && mdist < radius_ext) {
					c_number force_mod = extF.repulsivesphere.stiff*((c_number)1.f - radius/mdist);
					F.x += -dist.x*force_mod;
					F.y += -dist.y*force_mod;
					F.z += -dist.z*force_mod;
				}
				break;
			}
			case CUDA_REPULSIVE_SPHERE_SMOOTH: {
				c_number4 centre = make_c_number4(extF.repulsivespheresmooth.centre.x, extF.repulsivespheresmooth.centre.y, extF.repulsivespheresmooth.centre.z, 0.);
				c_number4 dist = box->minimum_image(centre, ppos);
				c_number mdist = _module(dist);
				c_number r0 = extF.repulsivespheresmooth.r0;
				c_number r_ext = extF.repulsivespheresmooth.r_ext;
				c_number smooth = extF.repulsivespheresmooth.smooth;
				c_number alpha = extF.repulsivespheresmooth.alpha;
				c_number stiff = extF.repulsivespheresmooth.stiff;

				if(mdist > r0 && mdist < r_ext) {
					c_number force_mod = 0.f;
					if(mdist >= alpha && mdist <= r_ext) {
						force_mod = -(stiff * 0.5f * expf((mdist - alpha) / smooth)) / mdist;
					}
					else {
						if(mdist < alpha && mdist >= r0) {
							force_mod = -(stiff * mdist - stiff * 0.5f * expf(-(mdist - alpha) / smooth)) / mdist;
						}
					}

					F.x += dist.x*force_mod;
					F.y += dist.y*force_mod;
					F.z += dist.z*force_mod;
				}

				break;
			}
			case CUDA_LJ_WALL: {
				c_number distance = CUDA_DOT(extF.ljwall.dir, ppos) + extF.ljwall.position;
				c_number rel_distance = distance/extF.ljwall.sigma;
				if(rel_distance < extF.ljwall.cutoff) {
					c_number lj_part = powf(rel_distance, -extF.ljwall.n);
					c_number factor = 4.f*extF.ljwall.n*extF.ljwall.stiff*(2.f*SQR(lj_part) - lj_part)/distance;
					F.x += factor*extF.ljwall.dir.x;
					F.y += factor*extF.ljwall.dir.y;
					F.z += factor*extF.ljwall.dir.z;
				}
				break;
			}
			case CUDA_CONSTANT_RATE_TORQUE: {
				c_number t = (extF.constantratetorque.F0 + extF.constantratetorque.rate * (c_number) step);

				c_number sintheta = sinf(t);
				c_number costheta = cosf(t);
				c_number olcos = (c_number)1.f - costheta;

				c_number xyo = extF.constantratetorque.axis.x * extF.constantratetorque.axis.y * olcos;
				c_number xzo = extF.constantratetorque.axis.x * extF.constantratetorque.axis.z * olcos;
				c_number yzo = extF.constantratetorque.axis.y * extF.constantratetorque.axis.z * olcos;
				c_number xsin = extF.constantratetorque.axis.x * sintheta;
				c_number ysin = extF.constantratetorque.axis.y * sintheta;
				c_number zsin = extF.constantratetorque.axis.z * sintheta;

				float3 to_rotate = extF.constantratetorque.pos0 - extF.constantratetorque.center;
				float3 postrap = make_float3(
						(extF.constantratetorque.axis.x * extF.constantratetorque.axis.x * olcos + costheta)*to_rotate.x + (xyo - zsin)*to_rotate.y + (xzo + ysin)*to_rotate.z,
						(xyo + zsin)*to_rotate.x + (extF.constantratetorque.axis.y * extF.constantratetorque.axis.y * olcos + costheta)*to_rotate.y + (yzo - xsin)*to_rotate.z,
						(xzo - ysin)*to_rotate.x + (yzo + xsin)*to_rotate.y + (extF.constantratetorque.axis.z * extF.constantratetorque.axis.z * olcos + costheta)*to_rotate.z
				);
				postrap += extF.constantratetorque.center;

				// we "mask" the resulting vector;
				F.x += -extF.constantratetorque.stiff*(ppos.x - postrap.x) * extF.constantratetorque.mask.x;
				F.y += -extF.constantratetorque.stiff*(ppos.y - postrap.y) * extF.constantratetorque.mask.y;
				F.z += -extF.constantratetorque.stiff*(ppos.z - postrap.z) * extF.constantratetorque.mask.z;
				break;
			}
			case CUDA_GENERIC_CENTRAL_FORCE: {
				c_number strength = extF.genericconstantforce.F0;
				c_number4 dir = make_c_number4(
						extF.genericconstantforce.x - ppos.x,
						extF.genericconstantforce.y - ppos.y,
						extF.genericconstantforce.z - ppos.z, (c_number) 0.f);

				c_number dist_sqr = CUDA_DOT(dir, dir);
				if(dist_sqr >= extF.genericconstantforce.inner_cut_off_sqr) {
					if(extF.genericconstantforce.outer_cut_off_sqr == 0.f || dist_sqr <= extF.genericconstantforce.outer_cut_off_sqr) {
						c_number dist = sqrtf(dist_sqr);
						dir.x /= dist;
						dir.y /= dist;
						dir.z /= dist;

						F.x += dir.x*strength;
						F.y += dir.y*strength;
						F.z += dir.z*strength;
					}
				}
				break;
			}
			case CUDA_LJ_CONE: {
				float3 v_from_apex = make_float3(
					ppos.x - extF.ljcone.pos0.x,
					ppos.y - extF.ljcone.pos0.y,
					ppos.z - extF.ljcone.pos0.z
				);
				c_number d_from_apex = _module(v_from_apex);

				c_number d_along_axis = CUDA_DOT(v_from_apex, extF.ljcone.dir);
				float3 v_along_axis = extF.ljcone.dir*d_along_axis;

				c_number beta = CUDA_LRACOS(CUDA_DOT(v_from_apex, v_along_axis)/(d_from_apex*d_along_axis));
				c_number gamma = extF.ljcone.alpha - beta;

				c_number d_from_cone = d_from_apex*sinf(gamma);
				c_number rel_distance = d_from_cone/extF.ljcone.sigma; // distance from the plane in units of _sigma

				if(rel_distance < extF.ljcone.cutoff && extF.ljcone.stiff > 0.f) {
					// now we compute the normal to the cone
					float3 v_from_axis = v_along_axis - v_from_apex;
					float3 normal = v_from_axis + extF.ljcone.dir*(_module(v_from_axis)/extF.ljcone.sin_alpha);

					c_number lj_part = powf(rel_distance, -extF.ljcone.n);
					c_number factor = 4*extF.ljcone.n*extF.ljcone.stiff*(2*SQR(lj_part) - lj_part) / (d_from_cone * _module(normal));

					F.x += factor*normal.x;
					F.y += factor*normal.y;
					F.z += factor*normal.z;
				}
				break;
			}
			case CUDA_REPULSIVE_ELLIPSOID: {
				c_number4 centre = make_c_number4(extF.repulsiveellipsoid.centre.x, extF.repulsiveellipsoid.centre.y, extF.repulsiveellipsoid.centre.z, 0.);
				c_number4 r_1 = make_c_number4(extF.repulsiveellipsoid.r_1.x, extF.repulsiveellipsoid.r_1.y, extF.repulsiveellipsoid.r_1.z, 0.);
				c_number4 r_2 = make_c_number4(extF.repulsiveellipsoid.r_2.x, extF.repulsiveellipsoid.r_2.y, extF.repulsiveellipsoid.r_2.z, 0.);
				c_number4 dist = box->minimum_image(centre, ppos);

				c_number internal_cut = SQR(dist.x) / SQR(r_2.x) + SQR(dist.y) / SQR(r_2.y) + SQR(dist.z) / SQR(r_2.z);
				c_number external_cut = SQR(dist.x) / SQR(r_1.x) + SQR(dist.y) / SQR(r_1.y) + SQR(dist.z) / SQR(r_1.z);

				if(internal_cut >= 1. || external_cut <= 1.) {
					c_number mdist = _module(dist);

					c_number force_part = extF.repulsiveellipsoid.stiff / mdist;

					F.x += -dist.x * force_part;
					F.y += -dist.y * force_part;
					F.z += -dist.z * force_part;
				}

				c_number mdist = _module(dist);
				c_number radius = extF.repulsivesphere.r0 + extF.repulsivesphere.rate*(c_number) step;
				c_number radius_ext = extF.repulsivesphere.r_ext;
				if(mdist > radius && mdist < radius_ext) {
					c_number force_mod = extF.repulsivesphere.stiff*((c_number)1.f - radius/mdist);
					F.x += -dist.x*force_mod;
					F.y += -dist.y*force_mod;
					F.z += -dist.z*force_mod;
				}
				break;
			}
			case CUDA_COM_FORCE: {
				c_number4 com = make_c_number4(0., 0., 0., 0.);
				for(int index = 0; index < extF.comforce.n_com; index++){
					int p_idx = extF.comforce.com_indexes[index];
					com += poss[p_idx];
				}
				com.x /= extF.comforce.n_com;
				com.y /= extF.comforce.n_com;
				com.z /= extF.comforce.n_com;

				c_number4 ref = make_c_number4(0., 0., 0., 0.);
				for(int index = 0; index < extF.comforce.n_ref; index++){
					int p_idx = extF.comforce.ref_indexes[index];
					ref += poss[p_idx];
				}
				ref.x /= extF.comforce.n_ref;
				ref.y /= extF.comforce.n_ref;
				ref.z /= extF.comforce.n_ref;

				c_number4 dr = ref - com;
				c_number dr_abs = _module(dr);
				c_number4 force = dr * ((dr_abs - (extF.comforce.r0 + extF.comforce.rate * step)) * extF.comforce.stiff / dr_abs) / extF.comforce.n_com;

				F.x += force.x;
				F.y += force.y;
				F.z += force.z;

				break;
			}
			case CUDA_LR_COM_TRAP: {
				c_number4 p1a_vec({0.f, 0.f, 0.f, 0.f});
				c_number4 p2a_vec({0.f, 0.f, 0.f, 0.f});
				for(int will = 0; will < extF.ltcomtrap.p1a_size; will++) {
					p1a_vec.x += poss[extF.ltcomtrap.p1a[will]].x / extF.ltcomtrap.p1a_size;
					p1a_vec.y += poss[extF.ltcomtrap.p1a[will]].y / extF.ltcomtrap.p1a_size;
					p1a_vec.z += poss[extF.ltcomtrap.p1a[will]].z / extF.ltcomtrap.p1a_size;
				}

				for(int will = 0; will < extF.ltcomtrap.p2a_size; will++) {
					p2a_vec.x += poss[extF.ltcomtrap.p2a[will]].x / extF.ltcomtrap.p2a_size;
					p2a_vec.y += poss[extF.ltcomtrap.p2a[will]].y / extF.ltcomtrap.p2a_size;
					p2a_vec.z += poss[extF.ltcomtrap.p2a[will]].z / extF.ltcomtrap.p2a_size;
				}

				c_number4 dr = p1a_vec - p2a_vec;
				c_number dr_mod = _module(dr);

				int ix_left = (int) ((dr_mod - extF.ltcomtrap.xmin) / extF.ltcomtrap.dX);
				int ix_right = ix_left + 1;

				// make sure that we never get out of boundaries. If we are then we keep the force 0
				c_number meta_Fx = 0.f;
				if(ix_left >= 0 && ix_right <= (extF.ltcomtrap.N_grid - 1)) {
					meta_Fx = -(extF.ltcomtrap.potential_grid[ix_right] - extF.ltcomtrap.potential_grid[ix_left]) / extF.ltcomtrap.dX;
				}

				c_number4 force = (dr) * (meta_Fx / dr_mod);

				if(extF.ltcomtrap.mode == 1) {
					F.x += force.x / extF.ltcomtrap.p1a_size;
					F.y += force.y / extF.ltcomtrap.p1a_size;
					F.z += force.z / extF.ltcomtrap.p1a_size;
				}
				else if(extF.ltcomtrap.mode == 2) {
					F.x -= force.x / extF.ltcomtrap.p2a_size;
					F.y -= force.y / extF.ltcomtrap.p2a_size;
					F.z -= force.z / extF.ltcomtrap.p2a_size;
				}
				break;
			}
			default: {
				break;
			}
		}
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void second_step(c_number4 *vels, c_number4 *Ls, c_number4 *forces, c_number4 *torques) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 v = vels[IND];

	v.x += F.x * MD_dt[0] * (c_number) 0.5f;
	v.y += F.y * MD_dt[0] * (c_number) 0.5f;
	v.z += F.z * MD_dt[0] * (c_number) 0.5f;
	v.w = (v.x*v.x + v.y*v.y + v.z*v.z) * (c_number) 0.5f;

	vels[IND] = v;

	c_number4 T = torques[IND];
	c_number4 L = Ls[IND];

	L.x += T.x * MD_dt[0] * (c_number) 0.5f;
	L.y += T.y * MD_dt[0] * (c_number) 0.5f;
	L.z += T.z * MD_dt[0] * (c_number) 0.5f;
	L.w = (L.x*L.x + L.y*L.y + L.z*L.z) * (c_number) 0.5f;

	Ls[IND] = L;
}
