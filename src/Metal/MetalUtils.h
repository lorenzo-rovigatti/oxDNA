/**
 * @file    MetalUtils.h
 * @brief   Utility functions for Metal GPU operations
 *
 * Metal equivalent of CUDAUtils.h
 */

#ifndef METALUTILS_H_
#define METALUTILS_H_

#include "../defs.h"
#include "metal_defs.h"

/**
 * @brief Static class for Metal GPU utility functions
 */
class MetalUtils {
protected:
    static size_t _allocated_dev_mem;

public:
    template<typename T> static void print_device_array(T *array, int count);
    template<typename T> static void check_device_threshold(T *array, int count, int threshold);

    /**
     * @brief Sum the 4th component of m_number4 vectors
     */
    static m_number sum_4th_comp(m_number4 *v, int N) {
        m_number res = 0;
        for(int i = 0; i < N; i++)
            res += v[i].w;
        return res;
    }

    /**
     * @brief Sum m_number4 vectors on GPU using reduction
     */
    static m_number4 sum_m_number4_on_GPU(id<MTLBuffer> buffer, int N, id<MTLDevice> device, id<MTLCommandQueue> queue);

    /**
     * @brief Sum m_number4 to double on GPU
     */
    static double sum_m_number4_to_double_on_GPU(id<MTLBuffer> buffer, int N, id<MTLDevice> device, id<MTLCommandQueue> queue);

    /**
     * @brief Convert int to float representation (bit-level)
     */
    static float int_as_float(const int a) {
        union {
            int a;
            float b;
        } u;
        u.a = a;
        return u.b;
    }

    /**
     * @brief Convert float to int representation (bit-level)
     */
    static int float_as_int(const float a) {
        union {
            float a;
            int b;
        } u;
        u.a = a;
        return u.b;
    }

    /**
     * @brief Get total allocated GPU memory in bytes
     */
    static size_t get_allocated_mem() {
        return _allocated_dev_mem;
    }

    /**
     * @brief Get total allocated GPU memory in MB
     */
    static double get_allocated_mem_mb() {
        return get_allocated_mem() / 1048576.;
    }

    /**
     * @brief Allocate Metal buffer and track memory usage
     */
    template<typename T>
    static id<MTLBuffer> allocate_buffer(id<MTLDevice> device, size_t count, MTLResourceOptions options = MTLResourceStorageModeShared);

    /**
     * @brief Reset memory tracking counter
     */
    static void reset_allocated_mem() {
        _allocated_dev_mem = 0;
    }

    /**
     * @brief Create Metal texture from buffer
     */
    template<typename T>
    static id<MTLTexture> create_texture_from_buffer(id<MTLDevice> device, id<MTLBuffer> buffer, size_t count, MTLPixelFormat format);

    /**
     * @brief Copy data from host to device buffer
     */
    template<typename T>
    static void copy_to_device(id<MTLBuffer> buffer, const T* host_data, size_t count);

    /**
     * @brief Copy data from device buffer to host
     */
    template<typename T>
    static void copy_from_device(T* host_data, id<MTLBuffer> buffer, size_t count);
};

template<typename T>
id<MTLBuffer> MetalUtils::allocate_buffer(id<MTLDevice> device, size_t count, MTLResourceOptions options) {
    size_t size = count * sizeof(T);
    OX_DEBUG("Allocating %zu bytes (%.2f MB) on the GPU", size, size / 1000000.0);

    MetalUtils::_allocated_dev_mem += size;

    id<MTLBuffer> buffer = [device newBufferWithLength:size options:options];
    METAL_CHECK_ERROR(buffer, "Failed to allocate Metal buffer");

    return buffer;
}

template<typename T>
void MetalUtils::copy_to_device(id<MTLBuffer> buffer, const T* host_data, size_t count) {
    if (!buffer || !host_data) {
        return;
    }

    T* device_ptr = static_cast<T*>([buffer contents]);
    memcpy(device_ptr, host_data, count * sizeof(T));
}

template<typename T>
void MetalUtils::copy_from_device(T* host_data, id<MTLBuffer> buffer, size_t count) {
    if (!buffer || !host_data) {
        return;
    }

    T* device_ptr = static_cast<T*>([buffer contents]);
    memcpy(host_data, device_ptr, count * sizeof(T));
}

template<typename T>
void MetalUtils::print_device_array(T *array, int count) {
    if (!array) return;

    printf("Device array contents (%d elements):\n", count);
    for(int i = 0; i < count && i < 100; i++) {
        printf("[%d] = %g\n", i, (double)array[i]);
    }
}

template<typename T>
void MetalUtils::check_device_threshold(T *array, int count, int threshold) {
    if (!array) return;

    for(int i = 0; i < count; i++) {
        if(array[i] > threshold) {
            printf("Warning: Element %d exceeds threshold: %g > %d\n", i, (double)array[i], threshold);
        }
    }
}

#endif /* METALUTILS_H_ */
