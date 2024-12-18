SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(rovigatti
	DEPENDS PolymerSwapInteraction MGAssemblyConf 
			DensityPressureProfile PolydisperseLTInteraction ConstructwisePressure 
			MicrogelElasticity MGAnalysis AOInteraction RadialDensityProfile MGInteraction 
			GenericGrByInsertion VoidPercolation Remoteness CPAnalysis Widom 
			CPMixtureInteraction StarrInteraction mWInteraction 
			FSInteraction StarrConf ConstructwiseBonds 
			EmptyVolume FSConf Gyradius PatchySwapInteraction DetailedPatchySwapInteraction 
			PatchyBonds Bonds DetailedPolymerSwapInteraction
)

# without this bit clang-based compilation fails with some "undefined symbols" errors during the linking stage
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	link_libraries(${ARGV0} "-undefined dynamic_lookup")
endif()

# this removes the prefix from each library. I didn't find any other way of doing it for all 
# libraries at once if the library type is MODULE
function(add_library_no_prefix target source)
	# ${ARGN} holds the list of arguments past the last expected argument, which is "source"
	add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
	set_target_properties(${target} PROPERTIES PREFIX "")
endfunction()

# Observables
add_library_no_prefix(Bonds src/Observables/Bonds.cpp src/Interactions/DetailedPolymerSwapInteraction.cpp)
add_library_no_prefix(PatchyBonds src/Observables/PatchyBonds.cpp)
add_library_no_prefix(ConstructwiseBonds src/Observables/ConstructwiseBonds.cpp)
add_library_no_prefix(EmptyVolume src/Observables/EmptyVolume.cpp)
add_library_no_prefix(FSConf src/Observables/FSConf.cpp src/Interactions/FSInteraction.cpp)
add_library_no_prefix(Gyradius src/Observables/Gyradius.cpp)
add_library_no_prefix(StarrConf src/Observables/StarrConf.cpp)
add_library_no_prefix(Widom src/Observables/Widom.cpp)
add_library_no_prefix(CPAnalysis src/Observables/CPAnalysis.cpp)
add_library_no_prefix(Remoteness src/Observables/Remoteness.cpp)
add_library_no_prefix(VoidPercolation src/Observables/VoidPercolation.cpp)
add_library_no_prefix(GenericGrByInsertion src/Observables/GenericGrByInsertion.cpp)
add_library_no_prefix(RadialDensityProfile src/Observables/RadialDensityProfile.cpp)
add_library_no_prefix(MGAnalysis src/Observables/MGAnalysis.cpp)
add_library_no_prefix(MicrogelElasticity src/Observables/MicrogelElasticity.cpp)
add_library_no_prefix(ConstructwisePressure src/Observables/ConstructwisePressure.cpp src/Interactions/PolymerSwapInteraction.cpp)
add_library_no_prefix(DensityPressureProfile src/Observables/DensityPressureProfile.cpp)
add_library_no_prefix(MGAssemblyConf src/Observables/MGAssemblyConf.cpp src/Interactions/FSInteraction.cpp)

# Interactions
add_library_no_prefix(FSInteraction src/Interactions/FSInteraction.cpp)
add_library_no_prefix(mWInteraction src/Interactions/mWInteraction.cpp)
add_library_no_prefix(StarrInteraction src/Interactions/StarrInteraction.cpp)
add_library_no_prefix(CPMixtureInteraction src/Interactions/CPMixtureInteraction.cpp)
add_library_no_prefix(MGInteraction src/Interactions/MGInteraction.cpp)
add_library_no_prefix(AOInteraction src/Interactions/AOInteraction.cpp)
add_library_no_prefix(PolydisperseLTInteraction src/Interactions/PolydisperseLTInteraction.cpp)
add_library_no_prefix(PolymerSwapInteraction src/Interactions/PolymerSwapInteraction.cpp)
add_library_no_prefix(PatchySwapInteraction src/Interactions/PatchySwapInteraction.cpp)
add_library_no_prefix(DetailedPatchySwapInteraction src/Interactions/DetailedPatchySwapInteraction.cpp)
add_library_no_prefix(DetailedPolymerSwapInteraction src/Interactions/DetailedPolymerSwapInteraction.cpp)

IF(CUDA)
	set(CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -rdynamic")
	find_package("CUDA")
	
	# same thing but for CUDA libs
	function(cuda_add_library_no_prefix target source)
		cuda_add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
		target_link_libraries(${target} ${CUDA_LIBRARIES})
		set_target_properties(${target} PROPERTIES PREFIX "")
	endfunction()

	cuda_add_library_no_prefix(CUDADetailedPatchySwapInteraction src/Interactions/CUDADetailedPatchySwapInteraction.cu src/Interactions/DetailedPatchySwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPatchySwapInteraction src/Interactions/CUDAPatchySwapInteraction.cu src/Interactions/PatchySwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDAFSInteraction src/Interactions/CUDAFSInteraction.cu src/Interactions/FSInteraction.cpp)
	cuda_add_library_no_prefix(CUDAmWInteraction src/Interactions/CUDAmWInteraction.cu src/Interactions/mWInteraction.cpp)
	cuda_add_library_no_prefix(CUDAStarrInteraction src/Interactions/CUDAStarrInteraction.cu src/Interactions/StarrInteraction.cpp)
	cuda_add_library_no_prefix(CUDACPMixtureInteraction src/Interactions/CUDACPMixtureInteraction.cu src/Interactions/CPMixtureInteraction.cpp)
	cuda_add_library_no_prefix(CUDAMGInteraction src/Interactions/CUDAMGInteraction.cu src/Interactions/MGInteraction.cpp)
	cuda_add_library_no_prefix(CUDAAOInteraction src/Interactions/CUDAAOInteraction.cu src/Interactions/AOInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPolymerSwapInteraction src/Interactions/CUDAPolymerSwapInteraction.cu src/Interactions/PolymerSwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDADetailedPolymerSwapInteraction src/Interactions/CUDADetailedPolymerSwapInteraction.cu src/Interactions/DetailedPolymerSwapInteraction.cpp)
	
	ADD_DEPENDENCIES(rovigatti CUDADetailedPolymerSwapInteraction CUDADetailedPatchySwapInteraction CUDAPatchySwapInteraction CUDAAOInteraction CUDAMGInteraction CUDACPMixtureInteraction CUDAFSInteraction CUDAmWInteraction CUDAStarrInteraction CUDAPolymerSwapInteraction)
ENDIF(CUDA)
