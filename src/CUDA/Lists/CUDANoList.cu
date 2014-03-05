/*
 * CUDANoList.cpp
 *
 *  Created on: 29/set/2010
 *      Author: lorenzo
 */

#include "CUDANoList.h"

template<typename number, typename number4>
CUDANoList<number, number4>::CUDANoList() {

}

template<typename number, typename number4>
CUDANoList<number, number4>::~CUDANoList() {

}

template class CUDANoList<float, float4>;
template class CUDANoList<double, LR_double4>;
