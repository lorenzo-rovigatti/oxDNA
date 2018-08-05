cmake_minimum_required(VERSION 2.8)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Interactions)
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Observables)
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Particles)

ADD_CUSTOM_TARGET(romano
	#	DEPENDS PatchyShapeParticle PatchyShapeInteraction DirkInteraction DirkInteraction2 DirkInteractionBias NematicS Icosahedron HardIcoInteraction
	DEPENDS HardIcoInteraction  PatchyShapeParticle PatchyShapeInteraction MCMovePatchyShape
) 

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

# Observables
ADD_LIBRARY(NematicS SHARED EXCLUDE_FROM_ALL src/Observables/NematicS.cpp)

# Interactions
ADD_LIBRARY(DirkInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction.cpp)
ADD_LIBRARY(DirkInteraction2 SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction2.cpp)
ADD_LIBRARY(DirkInteractionSin SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionSin.cpp)
ADD_LIBRARY(DirkInteractionBias SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionBias.cpp)
ADD_LIBRARY(HardIcoInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/HardIcoInteraction.cpp src/Particles/Icosahedron.cpp)
ADD_LIBRARY(PatchyShapeInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PatchyShapeInteraction.cpp   src/Particles/PatchyShapeParticle.cpp)


#Backends
ADD_LIBRARY(MCMovePatchyShape SHARED EXCLUDE_FROM_ALL  src/Backends/MCMoves/MCMovePatchyShape.cpp   src/Interactions/PatchyShapeInteraction.cpp src/Particles/PatchyShapeParticle.cpp)
# Particles
ADD_LIBRARY(PatchyShapeParticle SHARED EXCLUDE_FROM_ALL src/Particles/PatchyShapeParticle.cpp)
ADD_LIBRARY(Icosahedron SHARED EXCLUDE_FROM_ALL src/Particles/Icosahedron.cpp)
