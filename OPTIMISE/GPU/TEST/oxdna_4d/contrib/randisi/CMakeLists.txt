cmake_minimum_required(VERSION 3.0)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Interactions)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Observables)

ADD_CUSTOM_TARGET(randisi
	DEPENDS DNA2ModInteraction
) 

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

# Observables

# Interactions
ADD_LIBRARY(DNA2ModInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/DNA2ModInteraction.cpp)
IF(CUDA)
	CUDA_ADD_LIBRARY(CUDADNA2ModInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDADNA2ModInteraction.cu src/Interactions/DNA2ModInteraction.cpp)
	ADD_DEPENDENCIES(randisi CUDADNA2ModInteraction)
ENDIF(CUDA)
