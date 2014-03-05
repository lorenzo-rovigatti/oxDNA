/*
 * cuda_device_info.h
 *
 *  Created on: 30/lug/2009
 *      Author: lorenzo
 */

#ifndef CUDA_DEVICE_INFO_H_
#define CUDA_DEVICE_INFO_H_

int get_device_count();
void check_device_existance(int device);
cudaDeviceProp get_current_device_prop();
cudaDeviceProp get_device_prop(int device);
cudaError_t set_device(int device);

#endif /* CUDA_DEVICE_INFO_H_ */
