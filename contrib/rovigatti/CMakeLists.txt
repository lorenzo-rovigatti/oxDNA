SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(rovigatti
	DEPENDS PolymerSwapInteraction TSPInteraction MGAssemblyConf PolymerInteraction 
			DensityPressureProfile PolydisperseLTInteraction ConstructwisePressure 
			MicrogelElasticity MGAnalysis AOInteraction RadialDensityProfile MGInteraction 
			GenericGrByInsertion VoidPercolation Remoteness CPAnalysis Widom LevyDelta 
			LevyInteraction CPMixtureInteraction StarrInteraction mWInteraction 
			GraftedInteraction FSInteraction TSPAnalysis StarrConf ConstructwiseBonds 
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
add_library_no_prefix(TSPAnalysis src/Observables/TSPAnalysis.cpp)
add_library_no_prefix(Widom src/Observables/Widom.cpp)
add_library_no_prefix(LevyDelta src/Observables/LevyDelta.cpp src/Interactions/LevyInteraction.cpp)
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
add_library_no_prefix(TSPInteraction src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(PolymerInteraction src/Interactions/PolymerInteraction.cpp src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(FSInteraction src/Interactions/FSInteraction.cpp)
add_library_no_prefix(GraftedInteraction src/Interactions/GraftedInteraction.cpp src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(mWInteraction src/Interactions/mWInteraction.cpp)
add_library_no_prefix(StarrInteraction src/Interactions/StarrInteraction.cpp)
add_library_no_prefix(CPMixtureInteraction src/Interactions/CPMixtureInteraction.cpp)
add_library_no_prefix(LevyInteraction src/Interactions/LevyInteraction.cpp)
add_library_no_prefix(MGInteraction src/Interactions/MGInteraction.cpp)
add_library_no_prefix(AOInteraction src/Interactions/AOInteraction.cpp)
add_library_no_prefix(PolydisperseLTInteraction src/Interactions/PolydisperseLTInteraction.cpp)
add_library_no_prefix(PolymerSwapInteraction src/Interactions/PolymerSwapInteraction.cpp)
add_library_no_prefix(PatchySwapInteraction src/Interactions/PatchySwapInteraction.cpp)
add_library_no_prefix(DetailedPatchySwapInteraction src/Interactions/DetailedPatchySwapInteraction.cpp)
add_library_no_prefix(DetailedPolymerSwapInteraction src/Interactions/DetailedPolymerSwapInteraction.cpp)

FIND_PACKAGE(PkgConfig)
IF(PkgConfig_FOUND)
	PKG_CHECK_MODULES(GSL gsl)
	IF(GSL_FOUND)
		add_library_no_prefix(NathanNeighs src/Observables/NathanNeighs.cpp)
		target_link_libraries(NathanNeighs ${GSL_LIBRARIES})
	ENDIF(GSL_FOUND)
ENDIF(PkgConfig_FOUND)

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
	cuda_add_library_no_prefix(CUDALevyInteraction src/Interactions/CUDALevyInteraction.cu src/Interactions/LevyInteraction.cpp)
	cuda_add_library_no_prefix(CUDACPMixtureInteraction src/Interactions/CUDACPMixtureInteraction.cu src/Interactions/CPMixtureInteraction.cpp)
	cuda_add_library_no_prefix(CUDAMGInteraction src/Interactions/CUDAMGInteraction.cu src/Interactions/MGInteraction.cpp)
	cuda_add_library_no_prefix(CUDAAOInteraction src/Interactions/CUDAAOInteraction.cu src/Interactions/AOInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPolymerInteraction src/Interactions/CUDAPolymerInteraction.cu src/Interactions/PolymerInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPolymerSwapInteraction src/Interactions/CUDAPolymerSwapInteraction.cu src/Interactions/PolymerSwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDADetailedPolymerSwapInteraction src/Interactions/CUDADetailedPolymerSwapInteraction.cu src/Interactions/DetailedPolymerSwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDATSPInteraction src/Interactions/CUDATSPInteraction.cu src/Interactions/TSPInteraction.cpp)
	
	ADD_DEPENDENCIES(rovigatti CUDADetailedPolymerSwapInteraction CUDADetailedPatchySwapInteraction CUDAPatchySwapInteraction CUDAPolymerInteraction CUDAAOInteraction CUDAMGInteraction CUDACPMixtureInteraction CUDAFSInteraction CUDAmWInteraction CUDAStarrInteraction CUDALevyInteraction CUDAPolymerSwapInteraction CUDATSPInteraction)
ENDIF(CUDA)
