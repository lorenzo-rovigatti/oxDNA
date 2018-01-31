cmake_minimum_required(VERSION 2.8)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Interactions)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Observables)

ADD_CUSTOM_TARGET(romano
	DEPENDS DirkInteraction DirkInteraction2 DirkInteractionBias NematicS
) 

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

# Observables
ADD_LIBRARY(NematicS SHARED EXCLUDE_FROM_ALL src/Observables/NematicS.cpp)

# Interactions
ADD_LIBRARY(DirkInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction.cpp)
ADD_LIBRARY(DirkInteraction2 SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction2.cpp)
ADD_LIBRARY(DirkInteractionSin SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionSin.cpp)
ADD_LIBRARY(DirkInteractionBias SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionBias.cpp)

