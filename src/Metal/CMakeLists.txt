# CMakeLists.txt for Metal GPU backend
# Requires macOS and Metal-capable GPU (Apple Silicon or compatible)

IF(APPLE)
    MESSAGE(STATUS "Building Metal GPU support for Apple Silicon")

    # Find Metal framework
    FIND_LIBRARY(METAL_FRAMEWORK Metal REQUIRED)
    FIND_LIBRARY(FOUNDATION_FRAMEWORK Foundation REQUIRED)
    FIND_LIBRARY(COREGRAPHICS_FRAMEWORK CoreGraphics REQUIRED)

    IF(NOT METAL_FRAMEWORK)
        MESSAGE(FATAL_ERROR "Metal framework not found. Metal requires macOS 10.11 or later with Metal-capable GPU.")
    ENDIF()

    MESSAGE(STATUS "Metal framework: ${METAL_FRAMEWORK}")
    MESSAGE(STATUS "Foundation framework: ${FOUNDATION_FRAMEWORK}")

    # Set Objective-C++ compiler flags
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fobjc-arc")

    # Metal shader compilation
    SET(METAL_SHADERS
        Shaders/common.metal
        Shaders/md_kernels.metal
        Shaders/list_kernels.metal
        Shaders/dna_kernels.metal
        Shaders/thermostat_kernels.metal
    )

    # Compile Metal shaders to a library
    SET(METAL_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/shaders.metallib)

    # Find Metal compiler and linker explicitly using xcrun
    EXECUTE_PROCESS(COMMAND xcrun -f metal OUTPUT_VARIABLE METAL_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND xcrun -f metallib OUTPUT_VARIABLE METALLIB_LINKER OUTPUT_STRIP_TRAILING_WHITESPACE)
    
    MESSAGE(STATUS "Found Metal compiler: ${METAL_COMPILER}")

    IF(NOT EXISTS "${METAL_COMPILER}")
        MESSAGE(WARNING "Metal compiler not found via xcrun. Shader compilation may fail.")
        SET(METAL_COMPILER xcrun -sdk macosx metal)
        SET(METALLIB_LINKER xcrun -sdk macosx metallib)
    ENDIF()

    # Custom command to compile Metal shaders
    ADD_CUSTOM_COMMAND(
        OUTPUT ${METAL_LIBRARY}
        COMMAND ${METAL_COMPILER} -c -I ${CMAKE_CURRENT_SOURCE_DIR}/Shaders ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/common.metal -o ${CMAKE_CURRENT_BINARY_DIR}/common.air
        COMMAND ${METAL_COMPILER} -c -I ${CMAKE_CURRENT_SOURCE_DIR}/Shaders ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/md_kernels.metal -o ${CMAKE_CURRENT_BINARY_DIR}/md_kernels.air
        COMMAND ${METAL_COMPILER} -c -I ${CMAKE_CURRENT_SOURCE_DIR}/Shaders ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/list_kernels.metal -o ${CMAKE_CURRENT_BINARY_DIR}/list_kernels.air
        COMMAND ${METAL_COMPILER} -c -I ${CMAKE_CURRENT_SOURCE_DIR}/Shaders ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/dna_kernels.metal -o ${CMAKE_CURRENT_BINARY_DIR}/dna_kernels.air
        COMMAND ${METAL_COMPILER} -c -I ${CMAKE_CURRENT_SOURCE_DIR}/Shaders ${CMAKE_CURRENT_SOURCE_DIR}/Shaders/thermostat_kernels.metal -o ${CMAKE_CURRENT_BINARY_DIR}/thermostat_kernels.air
        COMMAND ${METALLIB_LINKER} ${CMAKE_CURRENT_BINARY_DIR}/common.air ${CMAKE_CURRENT_BINARY_DIR}/md_kernels.air ${CMAKE_CURRENT_BINARY_DIR}/list_kernels.air ${CMAKE_CURRENT_BINARY_DIR}/dna_kernels.air ${CMAKE_CURRENT_BINARY_DIR}/thermostat_kernels.air -o ${METAL_LIBRARY}
        DEPENDS ${METAL_SHADERS}
        COMMENT "Compiling Metal shaders"
        VERBATIM
    )

    ADD_CUSTOM_TARGET(metal_shaders ALL DEPENDS ${METAL_LIBRARY})

    # Copy the Metal library to the binary directory so it can be found by the executable
    ADD_CUSTOM_COMMAND(
        TARGET metal_shaders POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${METAL_LIBRARY}
                ${CMAKE_BINARY_DIR}/bin
    )

    # Metal backend source files
    SET(METAL_SOURCES
        MetalUtils.mm
        Backends/MetalBaseBackend.mm
        Backends/MD_MetalBackend.mm
        Lists/MetalNoList.mm
        Lists/MetalSimpleVerletList.mm
        Lists/MetalEdgeList.mm
        Lists/MetalListFactory.mm
        Interactions/MetalBaseInteraction.mm
        Interactions/MetalCPUForceFallback.mm
        Interactions/MetalDNAInteraction.mm
        Interactions/MetalLJInteraction.mm
        Interactions/MetalRNAInteraction.mm
        Interactions/MetalPatchyInteraction.mm
        Interactions/MetalTEPInteraction.mm
        Interactions/MetalInteractionFactory.mm
        Thermostats/MetalBaseThermostat.mm
        Thermostats/MetalBrownianThermostat.mm
        Thermostats/MetalThermostatFactory.mm
        MetalBackendFactory.mm
        ../Utilities/Timings.cpp
        # Add more backend implementations here as needed
    )

    # Prepend absolute path to source files so they can be found by the parent CMakeLists.txt
    LIST(TRANSFORM METAL_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

    # Add Metal sources to the common library
    # Note: These will be added to the main project's common library
    SET(metal_common_SOURCES
        ${METAL_SOURCES}
        PARENT_SCOPE
    )

    # Note: there is intentionally no double-precision switch here. Apple GPUs
    # have no hardware double precision, so shader storage is always float32.
    # Accuracy tiers are selected at run time with `backend_precision`.

    # Include directories
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    # Framework linking will be handled in the main CMakeLists.txt
    SET(METAL_FRAMEWORKS
        ${METAL_FRAMEWORK}
        ${FOUNDATION_FRAMEWORK}
        ${COREGRAPHICS_FRAMEWORK}
        PARENT_SCOPE
    )

    MESSAGE(STATUS "Metal backend configured successfully")

    # Install Metal shader library
    INSTALL(FILES ${METAL_LIBRARY} DESTINATION bin)

ELSE()
    MESSAGE(STATUS "Metal backend is only available on macOS - skipping")
ENDIF(APPLE)
