/* Verlet related constants */
__constant__ float MD_sqr_verlet_skin[1];
/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_dt[1];
__constant__ float MD_box_side[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

// this function modifies L
// Conversion of the R matrix in _get_updated_orientation above into a quaternion simplifies into the qR quaternion used in _get_updated_orientation. This is then multiplied by the original orientation quaternion to produce the updated quaternion. Quaternion multiplication is the cross_product - dot_product and is therefore non-commutative
template <typename number, typename number4>
__device__ GPU_quat<number> _get_updated_orientation(number4 &L, GPU_quat<number> &old_o) {
	number norm = _module<number, number4>(L);
	L.x /= norm;
	L.y /= norm;
	L.z /= norm;

	number sintheta, costheta;
	sincos(MD_dt[0] * norm, &sintheta, &costheta);
	number qw = (number) 0.5f * sqrtf(fmaxf((number) 0.f, (number) 2.f + 2.f*costheta));
	number winv = (number)1.f /qw;
	GPU_quat<number> R = {(number) 0.5f*L.x*sintheta*winv, (number) 0.5f*L.y*sintheta*winv, (number) 0.5f*L.z*sintheta*winv, qw};

	return quat_multiply(old_o, R);
}

template <typename number, typename number4>
__global__ void first_step(number4 *poss, GPU_quat<number> *orientations, number4 *list_poss, number4 *vels, number4 *Ls, number4 *forces, number4 *torques, bool *are_lists_old) {
	if(IND >= MD_N[0]) return;

	const number4 F = forces[IND];

	number4 r = poss[IND];
	number4 v = vels[IND];

	v.x += F.x * (MD_dt[0] * (number) 0.5f);
	v.y += F.y * (MD_dt[0] * (number) 0.5f);
	v.z += F.z * (MD_dt[0] * (number) 0.5f);

	r.x += v.x * MD_dt[0];
	r.y += v.y * MD_dt[0];
	r.z += v.z * MD_dt[0];

	vels[IND] = v;
	poss[IND] = r;

	const number4 T = torques[IND];
	number4 L = Ls[IND];

	L.x += T.x * (MD_dt[0] * (number) 0.5f);
	L.y += T.y * (MD_dt[0] * (number) 0.5f);
	L.z += T.z * (MD_dt[0] * (number) 0.5f);

	Ls[IND] = L;


	GPU_quat<number> qold_o = orientations[IND];
	orientations[IND] = _get_updated_orientation<number>(L, qold_o);

	// do verlet lists need to be updated?
	if(quad_distance<number, number4>(r, list_poss[IND]) > MD_sqr_verlet_skin[0]) are_lists_old[0] = true;
}

template <typename number, typename number4>
__global__ void set_external_forces(number4 *poss, GPU_quat<number> *orientations, CUDA_trap<number> *ext_forces, number4 *forces, number4 *torques, llint step, int max_ext_forces) {
	if(IND >= MD_N[0]) return;
	// if there are no external forces then just put the force to 0
	if(max_ext_forces == 0) {
		forces[IND] = make_number4<number, number4>(0, 0, 0, 0);
		torques[IND] = make_number4<number, number4>(0, 0, 0, 0);
		return;
	}

	number4 ppos = poss[IND];
	number4 F = make_number4<number, number4>(0, 0, 0, 0);
	number4 T = make_number4<number, number4>(0, 0, 0, 0);

	for(int i = 0; i < max_ext_forces; i++) {
		// coalesced
		CUDA_trap<number> extF = ext_forces[MD_N[0] * i + IND];
		switch (extF.type) {
			case CUDA_TRAP_CONSTANT: {
				number intensity = extF.constant.F0 + step * extF.constant.rate;
				F.x += extF.constant.x * intensity;
				F.y += extF.constant.y * intensity;
				F.z += extF.constant.z * intensity;
				break;
			}
			case CUDA_TRAP_MUTUAL: {
				number4 qpos = poss[extF.mutual.p_ind];

				number dx = qpos.x - ppos.x;
				number dy = qpos.y - ppos.y;
				number dz = qpos.z - ppos.z;

				if (extF.mutual.PBC) {
					dx -= floorf(dx/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
					dy -= floorf(dy/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
					dz -= floorf(dz/MD_box_side[0] + (number) 0.5f) * MD_box_side[0];
				}

				number4 dr = make_number4<number, number4>(dx, dy, dz, (number) 0.f);
				number dr_abs = _module<number, number4>(dr);

				number max_dr = (number) 5.f;
				if (dr_abs > max_dr) {
					dr = dr * (max_dr / dr_abs);
					dr_abs = max_dr;
				}

				number4 force = dr / dr_abs * (dr_abs - extF.mutual.r0) * extF.mutual.stiff;

				F.x += force.x;
				F.y += force.y;
				F.z += force.z;
				break;
			}
			case CUDA_TRAP_MOVING: {
				number tx = extF.moving.pos0.x + extF.moving.rate * step * extF.moving.dir.x;
				number ty = extF.moving.pos0.y + extF.moving.rate * step * extF.moving.dir.y;
				number tz = extF.moving.pos0.z + extF.moving.rate * step * extF.moving.dir.z;

				F.x += - extF.moving.stiff * (ppos.x - tx);
				F.y += - extF.moving.stiff * (ppos.y - ty);
				F.z += - extF.moving.stiff * (ppos.z - tz);
				break;
			}
			case CUDA_TRAP_MOVING_LOWDIM: {
				number tx = extF.lowdim.pos0.x + extF.lowdim.rate * step * extF.lowdim.dir.x;
				number ty = extF.lowdim.pos0.y + extF.lowdim.rate * step * extF.lowdim.dir.y;
				number tz = extF.lowdim.pos0.z + extF.lowdim.rate * step * extF.lowdim.dir.z;

				if(extF.lowdim.visX) F.x += - extF.lowdim.stiff * (ppos.x - tx);
				if(extF.lowdim.visY) F.y += - extF.lowdim.stiff * (ppos.y - ty);
				if(extF.lowdim.visZ) F.z += - extF.lowdim.stiff * (ppos.z - tz);
				break;
			}
			case CUDA_REPULSION_PLANE: {
				number distance = extF.repulsionplane.dir.x*ppos.x + extF.repulsionplane.dir.y*ppos.y + extF.repulsionplane.dir.z*ppos.z + extF.repulsionplane.position;
				if(distance < 0.0f) {
					F.x += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.x;
					F.y += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.y;
					F.z += -distance * extF.repulsionplane.stiff * extF.repulsionplane.dir.z;
				}
				break;
			}
			case CUDA_REPULSION_PLANE_MOVING: {
				number4 qpos = poss[extF.repulsionplanemoving.p_ind];
				number distance = extF.repulsionplanemoving.dir.x*(ppos.x-qpos.x) + extF.repulsionplanemoving.dir.y*(ppos.y-qpos.y) + extF.repulsionplanemoving.dir.z*(ppos.z-qpos.z);
				if(distance < 0.0f) {
					F.x += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.x;
					F.y += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.y;
					F.z += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.z;
				}
				break;
			}
			case CUDA_LJ_WALL: {
				number distance = extF.ljwall.dir.x*ppos.x + extF.ljwall.dir.y*ppos.y + extF.ljwall.dir.z*ppos.z + extF.ljwall.position;
				number rel_distance = distance/extF.ljwall.sigma;
				if(rel_distance < extF.ljwall.cutoff) {
					number lj_part = powf(rel_distance, -extF.ljwall.n);
					number factor = 4.f*extF.ljwall.n*extF.ljwall.stiff*(2.f*SQR(lj_part) - lj_part)/distance;
					F.x += factor*extF.ljwall.dir.x;
					F.y += factor*extF.ljwall.dir.y;
					F.z += factor*extF.ljwall.dir.z;
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

template <typename number, typename number4>
__global__ void second_step(number4 *vels, number4 *Ls, number4 *forces, number4 *torques) {
	if(IND >= MD_N[0]) return;

	number4 F = forces[IND];
	number4 v = vels[IND];

	v.x += F.x * MD_dt[0] * (number) 0.5f;
	v.y += F.y * MD_dt[0] * (number) 0.5f;
	v.z += F.z * MD_dt[0] * (number) 0.5f;
	v.w = (v.x*v.x + v.y*v.y + v.z*v.z) * (number) 0.5f;

	vels[IND] = v;

	number4 T = torques[IND];
	number4 L = Ls[IND];

	L.x += T.x * MD_dt[0] * (number) 0.5f;
	L.y += T.y * MD_dt[0] * (number) 0.5f;
	L.z += T.z * MD_dt[0] * (number) 0.5f;
	L.w = (L.x*L.x + L.y*L.y + L.z*L.z) * (number) 0.5f;

	Ls[IND] = L;
}
