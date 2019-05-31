/* Verlet related constants */
__constant__ float MD_sqr_verlet_skin[1];
/* System constants */
__constant__ int MD_N[1];
__constant__ float MD_dt[1];

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

template<typename number, typename number4>
__global__ void rescale_positions(number4 *poss, number4 ratio) {
	if(IND >= MD_N[0]) return;
	number4 ppos = poss[IND];
	ppos.x *= ratio.x;
	ppos.y *= ratio.y;
	ppos.z *= ratio.z;
	poss[IND] = ppos;
}

template<typename number, typename number4>
__global__ void set_external_forces(number4 *poss, GPU_quat<number> *orientations, CUDA_trap<number> *ext_forces, number4 *forces, number4 *torques, llint step, int max_ext_forces, CUDABox<number, number4> *box) {
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
				number intensity = extF.constant.F0 + step*extF.constant.rate;
				float3 dir = make_float3(extF.constant.x, extF.constant.y, extF.constant.z);
				if(extF.constant.dir_as_centre) {
					dir.x -= ppos.x;
					dir.y -= ppos.y;
					dir.z -= ppos.z;
					number dir_mod = _module<number>(dir);
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
				number4 qpos = poss[extF.mutual.p_ind];

				number4 dr = (extF.mutual.PBC) ? box->minimum_image(ppos, qpos) : qpos - ppos;
				number dr_abs = _module<number, number4>(dr);

				number4 force = dr * ((dr_abs - (extF.mutual.r0 + extF.mutual.rate * step)) * extF.mutual.stiff / dr_abs);

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
				for(int idx = extF.repulsionplanemoving.low_idx; idx <= extF.repulsionplanemoving.high_idx; idx++) {
					number4 qpos = poss[idx];
					number distance = extF.repulsionplanemoving.dir.x*(ppos.x-qpos.x) + extF.repulsionplanemoving.dir.y*(ppos.y-qpos.y) + extF.repulsionplanemoving.dir.z*(ppos.z-qpos.z);
					if(distance < 0.0f) {
						F.x += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.x;
						F.y += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.y;
						F.z += -distance * extF.repulsionplanemoving.stiff * extF.repulsionplanemoving.dir.z;
					}
				}
				break;
			}
			case CUDA_REPULSIVE_SPHERE: {
				number4 centre = make_number4<number, number4>(extF.repulsivesphere.centre.x, extF.repulsivesphere.centre.y, extF.repulsivesphere.centre.z, 0.);
				number4 dist = box->minimum_image(centre, ppos);
				number mdist = _module<number, number4>(dist);
				number radius = extF.repulsivesphere.r0 + extF.repulsivesphere.rate*(number) step;
				number radius_ext = extF.repulsivesphere.r_ext;
				if(mdist > radius && mdist < radius_ext) {
					number force_mod = extF.repulsivesphere.stiff*((number)1.f - radius/mdist);
					F.x += -dist.x*force_mod;
					F.y += -dist.y*force_mod;
					F.z += -dist.z*force_mod;
				}
				break;
			}
			case CUDA_REPULSIVE_SPHERE_SMOOTH: {
				number4 centre = make_number4<number, number4>(extF.repulsivespheresmooth.centre.x, extF.repulsivespheresmooth.centre.y, extF.repulsivespheresmooth.centre.z, 0.);
				number4 dist = box->minimum_image(centre, ppos);
				number mdist = _module<number, number4>(dist);
				number r0 = extF.repulsivespheresmooth.r0;
				number r_ext = extF.repulsivespheresmooth.r_ext;
				number smooth = extF.repulsivespheresmooth.smooth;
				number alpha = extF.repulsivespheresmooth.alpha;
				number stiff = extF.repulsivespheresmooth.stiff;

				if(mdist > r0 && mdist < r_ext) {
					number force_mod = 0.f;
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
				number distance = CUDA_DOT(extF.ljwall.dir, ppos) + extF.ljwall.position;
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
			case CUDA_CONSTANT_RATE_TORQUE: {
				number t = (extF.constantratetorque.F0 + extF.constantratetorque.rate * (number) step);

				number sintheta = sinf(t);
				number costheta = cosf(t);
				number olcos = (number)1.f - costheta;

				number xyo = extF.constantratetorque.axis.x * extF.constantratetorque.axis.y * olcos;
				number xzo = extF.constantratetorque.axis.x * extF.constantratetorque.axis.z * olcos;
				number yzo = extF.constantratetorque.axis.y * extF.constantratetorque.axis.z * olcos;
				number xsin = extF.constantratetorque.axis.x * sintheta;
				number ysin = extF.constantratetorque.axis.y * sintheta;
				number zsin = extF.constantratetorque.axis.z * sintheta;

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
				number strength = extF.genericconstantforce.F0;
				number4 dir = make_number4<number, number4>(
						extF.genericconstantforce.x - ppos.x,
						extF.genericconstantforce.y - ppos.y,
						extF.genericconstantforce.z - ppos.z, (number) 0.f);

				number dist_sqr = CUDA_DOT(dir, dir);
				if(dist_sqr >= extF.genericconstantforce.inner_cut_off_sqr) {
					if(extF.genericconstantforce.outer_cut_off_sqr == 0.f || dist_sqr <= extF.genericconstantforce.outer_cut_off_sqr) {
						number dist = sqrtf(dist_sqr);
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
				number d_from_apex = _module<number>(v_from_apex);

				number d_along_axis = CUDA_DOT(v_from_apex, extF.ljcone.dir);
				float3 v_along_axis = extF.ljcone.dir*d_along_axis;

				number beta = CUDA_LRACOS(CUDA_DOT(v_from_apex, v_along_axis)/(d_from_apex*d_along_axis));
				number gamma = extF.ljcone.alpha - beta;

				number d_from_cone = d_from_apex*sinf(gamma);
				number rel_distance = d_from_cone/extF.ljcone.sigma; // distance from the plane in units of _sigma

				if(rel_distance < extF.ljcone.cutoff && extF.ljcone.stiff > 0.f) {
					// now we compute the normal to the cone
					float3 v_from_axis = v_along_axis - v_from_apex;
					float3 normal = v_from_axis + extF.ljcone.dir*(_module<number>(v_from_axis)/extF.ljcone.sin_alpha);

					number lj_part = powf(rel_distance, -extF.ljcone.n);
					number factor = 4*extF.ljcone.n*extF.ljcone.stiff*(2*SQR(lj_part) - lj_part) / (d_from_cone * _module<number>(normal));

					F.x += factor*normal.x;
					F.y += factor*normal.y;
					F.z += factor*normal.z;
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
