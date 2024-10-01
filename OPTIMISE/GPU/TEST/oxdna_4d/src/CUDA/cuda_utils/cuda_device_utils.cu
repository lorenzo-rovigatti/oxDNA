#include "cuda_device_utils.h"

#include <cstdlib>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime_api.h>

int get_device_count() {
	int deviceCount = 0;
	if(cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
		fprintf(stderr, "cudaGetDeviceCount FAILED, CUDA Driver and Runtime CUDA Driver and Runtime version may be mismatched, exiting.\n");
		exit(1);
	}

	return deviceCount;
}

void check_device_existance(int device) {
	if(device >= get_device_count()) {
		fprintf(stderr, "The selected device doesn't exist, exiting.\n");
		exit(1);
	}
}

cudaDeviceProp get_current_device_prop() {
	int curr_dev;
	cudaGetDevice(&curr_dev);
	return get_device_prop(curr_dev);
}

cudaDeviceProp get_device_prop(int device) {
	check_device_existance(device);

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);

	return deviceProp;
}
