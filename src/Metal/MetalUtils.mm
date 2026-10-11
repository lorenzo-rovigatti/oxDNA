/**
 * @file    MetalUtils.mm
 * @brief   Implementation of Metal utility functions
 */

#include "MetalUtils.h"
#include <Foundation/Foundation.h>

size_t MetalUtils::_allocated_dev_mem = 0;

m_number4 MetalUtils::sum_m_number4_on_GPU(id<MTLBuffer> buffer, int N, id<MTLDevice> device, id<MTLCommandQueue> queue) {
    // For now, implement CPU-side reduction
    // TODO: Implement proper GPU reduction kernel
    m_number4 *data = static_cast<m_number4*>([buffer contents]);
    m_number4 result = {0, 0, 0, 0};

    for(int i = 0; i < N; i++) {
        result.x += data[i].x;
        result.y += data[i].y;
        result.z += data[i].z;
        result.w += data[i].w;
    }

    return result;
}

double MetalUtils::sum_m_number4_to_double_on_GPU(id<MTLBuffer> buffer, int N, id<MTLDevice> device, id<MTLCommandQueue> queue) {
    // For now, implement CPU-side reduction
    // TODO: Implement proper GPU reduction kernel
    m_number4 *data = static_cast<m_number4*>([buffer contents]);
    double result = 0.0;

    for(int i = 0; i < N; i++) {
        result += (double)(data[i].x + data[i].y + data[i].z + data[i].w);
    }

    return result;
}

template<typename T>
id<MTLTexture> MetalUtils::create_texture_from_buffer(id<MTLDevice> device, id<MTLBuffer> buffer, size_t count, MTLPixelFormat format) {
    MTLTextureDescriptor *descriptor = [MTLTextureDescriptor new];
    descriptor.pixelFormat = format;
    descriptor.width = count;
    descriptor.height = 1;
    descriptor.textureType = MTLTextureType1D;
    descriptor.usage = MTLTextureUsageShaderRead;

    id<MTLTexture> texture = [buffer newTextureWithDescriptor:descriptor
                                                        offset:0
                                                   bytesPerRow:count * sizeof(T)];

    METAL_CHECK_ERROR(texture, "Failed to create texture from buffer");
    return texture;
}
