/**
 * @file    MetalBackendFactory.mm
 * @brief   Factory function to create Metal backend
 *
 * This wrapper avoids including Objective-C++ headers in C++ files
 */

#include "Backends/MD_MetalBackend.h"
#include "../Backends/SimBackend.h"

// Factory function to create Metal backend
SimBackend* create_metal_backend() {
    return new MD_MetalBackend();
}
