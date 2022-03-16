cmake_minimum_required(VERSION 3.0)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Interactions)
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Observables)
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/Particles)

# PREVIOUSLY
#ADD_CUSTOM_TARGET(romano
	#	DEPENDS PatchyShapeParticle PatchyShapeInteraction DirkInteraction DirkInteraction2 DirkInteractionBias NematicS Icosahedron HardIcoInteraction
	#	DEPENDS HardIcoInteraction PLCluster PatchyShapeParticle PatchyShapeInteraction MCMovePatchyShape ChiralRodInteraction NematicS Swim ChiralRodExplicit Reappear CutVolume Grow Exhaust FakePressure FreeVolume Depletion NDepletion DepletionVolume AVBDepletion Swap #MCMoveDesign
	#) 

ADD_CUSTOM_TARGET(romano
	#	DEPENDS PatchyShapeParticle PatchyShapeInteraction DirkInteraction DirkInteraction2 DirkInteractionBias NematicS Icosahedron HardIcoInteraction
	DEPENDS Reappear ChargedCylInteraction SideChange NematicS ClosestN
) 

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

# Observables
ADD_LIBRARY(ClosestN SHARED EXCLUDE_FROM_ALL src/Observables/ClosestN.cpp)
ADD_LIBRARY(NematicS SHARED EXCLUDE_FROM_ALL src/Observables/NematicS.cpp)
ADD_LIBRARY(FakePressure SHARED EXCLUDE_FROM_ALL src/Observables/FakePressure.cpp)
ADD_LIBRARY(FreeVolume SHARED EXCLUDE_FROM_ALL src/Observables/FreeVolume.cpp)
ADD_LIBRARY(PLCluster SHARED EXCLUDE_FROM_ALL  src/Observables/PLCluster.cpp  src/Interactions/PatchyShapeInteraction.cpp src/Particles/PatchyShapeParticle.cpp)

# Interactions
ADD_LIBRARY(DirkInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction.cpp)
ADD_LIBRARY(DirkInteraction2 SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteraction2.cpp)
ADD_LIBRARY(DirkInteractionSin SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionSin.cpp)
ADD_LIBRARY(DirkInteractionBias SHARED EXCLUDE_FROM_ALL src/Interactions/DirkInteractionBias.cpp)
ADD_LIBRARY(ChiralRodInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/ChiralRodInteraction.cpp ../../src/Particles/SpheroCylinder.cpp)
ADD_LIBRARY(ChiralRodExplicit SHARED EXCLUDE_FROM_ALL src/Interactions/ChiralRodExplicit.cpp ../../src/Particles/SpheroCylinder.cpp)
ADD_LIBRARY(HardIcoInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/HardIcoInteraction.cpp src/Particles/Icosahedron.cpp)
ADD_LIBRARY(PatchyShapeInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PatchyShapeInteraction.cpp   src/Particles/PatchyShapeParticle.cpp)
ADD_LIBRARY(ChargedCylInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/ChargedCylInteraction.cpp ../../src/Particles/SpheroCylinder.cpp)

#Backends
ADD_LIBRARY(MCMovePatchyShape SHARED EXCLUDE_FROM_ALL  src/Backends/MCMoves/MCMovePatchyShape.cpp   src/Interactions/PatchyShapeInteraction.cpp src/Particles/PatchyShapeParticle.cpp)
ADD_LIBRARY(Swim SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Swim.cpp)
ADD_LIBRARY(Reappear SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Reappear.cpp)
ADD_LIBRARY(SideChange SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/SideChange.cpp)
ADD_LIBRARY(Depletion SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Depletion.cpp)
ADD_LIBRARY(NDepletion SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/NDepletion.cpp)
ADD_LIBRARY(CutVolume SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/CutVolume.cpp)
ADD_LIBRARY(Grow SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Grow.cpp)
ADD_LIBRARY(Exhaust SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Exhaust.cpp)
ADD_LIBRARY(DepletionVolume SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/DepletionVolume.cpp)
#ADD_LIBRARY(VMMCDepletion SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/VMMCDepletion.cpp)
ADD_LIBRARY(AVBDepletion SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/AVBDepletion.cpp)
#ADD_LIBRARY(OverlapVolume SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/OverlapVolume.cpp)
ADD_LIBRARY(Swap SHARED EXCLUDE_FROM_ALL src/Backends/MCMoves/Swap.cpp)
#ADD_LIBRARY(MCMoveDesign SHARED EXCLUDE_FROM_ALL  src/Backends/MCMoves/MCMoveDesign.cpp   src/Interactions/PatchyShapeInteraction.cpp src/Particles/PatchyShapeParticle.cpp)

# Particles
ADD_LIBRARY(PatchyShapeParticle SHARED EXCLUDE_FROM_ALL src/Particles/PatchyShapeParticle.cpp)
ADD_LIBRARY(Icosahedron SHARED EXCLUDE_FROM_ALL src/Particles/Icosahedron.cpp)
