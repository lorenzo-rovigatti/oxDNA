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

template<typename number4>
__global__ void hilbert_curve(const number4 *pos, int *hindex);

template<typename number, typename number4>
__global__ void permute_particles(int *sorted_hindex, number4 *poss, number4 *vels, number4 *buff_poss, number4 *buff_vels);

template<typename number, typename number4>
__global__ void permute_particles(int *sorted_hindex, number4 *poss, number4 *buff_poss);

template<typename number, typename number4>
__global__ void permute_particles(int *sorted_hindex, int *inv_sorted_hindex, number4 *poss, number4 *vels, number4 *Ls, GPU_quat<number> *orientations, LR_bonds *bonds, number4 *buff_poss, number4 *buff_vels, number4 *buff_Ls, GPU_quat<number> *buff_orientations, LR_bonds *buff_bonds);
__global__ void get_inverted_sorted_hindex(int *sorted_hindex, int *inv_sorted_hindex);
__global__ void reset_sorted_hindex(int *sorted_hindex);

void init_hilb_symbols(int N, int N_unsortable, int depth, float box_side);

#endif /* CUDA_SORT_H_ */
