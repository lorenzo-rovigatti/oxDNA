/*
 * CUDA_sort.cuh
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#ifndef CUDA_SORT_H_
#define CUDA_SORT_H_

#include "CUDAUtils.h"
//#include <cutil.h>


__global__ void hilbert_curve(const c_number4 *pos, int *hindex);


__global__ void permute_particles(int *sorted_hindex, c_number4 *poss, c_number4 *vels, c_number4 *buff_poss, c_number4 *buff_vels);


__global__ void permute_particles(int *sorted_hindex, c_number4 *poss, c_number4 *buff_poss);


__global__ void permute_particles(int *sorted_hindex, int *inv_sorted_hindex, c_number4 *poss, c_number4 *vels, c_number4 *Ls, GPU_quat *orientations, LR_bonds *bonds, int *particles_to_mols, c_number4 *buff_poss, c_number4 *buff_vels, c_number4 *buff_Ls, GPU_quat *buff_orientations, LR_bonds *buff_bonds, int *buff_particles_to_mols);
__global__ void get_inverted_sorted_hindex(int *sorted_hindex, int *inv_sorted_hindex);
__global__ void reset_sorted_hindex(int *sorted_hindex);

void init_hilb_symbols(int N, int N_unsortable, int depth, float box_side);

#endif /* CUDA_SORT_H_ */
