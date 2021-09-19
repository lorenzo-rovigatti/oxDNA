SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(rovigatti
	DEPENDS PolymerSwapInteraction TSPInteraction MGAssemblyConf PolymerInteraction 
			DensityPressureProfile PolydisperseLTInteraction ConstructwisePressure 
			MicrogelElasticity MGAnalysis AOInteraction RadialDensityProfile MGInteraction 
			GenericGrByInsertion VoidPercolation Remoteness CPAnalysis Widom LevyDelta 
			LevyInteraction CPMixtureInteraction StarrInteraction mWInteraction 
			GraftedInteraction FSInteraction TSPAnalysis StarrConf ConstructwiseBonds 
			EmptyVolume FSConf Gyradius PatchySwapInteraction DetailedPatchySwapInteraction 
			PatchyBonds Bonds ValerioInteraction
)

# without this bit clang-based compilation fails with some "undefined symbols" errors during the linking stage
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	link_libraries(${ARGV0} "-undefined dynamic_lookup")
endif()

# this removes the prefix from each library. I didn't find any other way of doing it for all 
# libraries at once if the library type is MODULE
function(add_library_no_prefix)
	add_library(${ARGN})
	set_target_properties(${ARGV0} PROPERTIES PREFIX "")
endfunction()

# Observables
add_library_no_prefix(Bonds MODULE EXCLUDE_FROM_ALL src/Observables/Bonds.cpp)
add_library_no_prefix(PatchyBonds MODULE EXCLUDE_FROM_ALL src/Observables/PatchyBonds.cpp)
add_library_no_prefix(ConstructwiseBonds MODULE EXCLUDE_FROM_ALL src/Observables/ConstructwiseBonds.cpp)
add_library_no_prefix(EmptyVolume MODULE EXCLUDE_FROM_ALL src/Observables/EmptyVolume.cpp)
add_library_no_prefix(FSConf MODULE EXCLUDE_FROM_ALL src/Observables/FSConf.cpp src/Interactions/FSInteraction.cpp)
add_library_no_prefix(Gyradius MODULE EXCLUDE_FROM_ALL src/Observables/Gyradius.cpp)
add_library_no_prefix(StarrConf MODULE EXCLUDE_FROM_ALL src/Observables/StarrConf.cpp)
add_library_no_prefix(TSPAnalysis MODULE EXCLUDE_FROM_ALL src/Observables/TSPAnalysis.cpp)
add_library_no_prefix(Widom MODULE EXCLUDE_FROM_ALL src/Observables/Widom.cpp)
add_library_no_prefix(LevyDelta MODULE EXCLUDE_FROM_ALL src/Observables/LevyDelta.cpp src/Interactions/LevyInteraction.cpp)
add_library_no_prefix(CPAnalysis MODULE EXCLUDE_FROM_ALL src/Observables/CPAnalysis.cpp)
add_library_no_prefix(Remoteness MODULE EXCLUDE_FROM_ALL src/Observables/Remoteness.cpp)
add_library_no_prefix(VoidPercolation MODULE EXCLUDE_FROM_ALL src/Observables/VoidPercolation.cpp)
add_library_no_prefix(GenericGrByInsertion MODULE EXCLUDE_FROM_ALL src/Observables/GenericGrByInsertion.cpp)
add_library_no_prefix(RadialDensityProfile MODULE EXCLUDE_FROM_ALL src/Observables/RadialDensityProfile.cpp)
add_library_no_prefix(MGAnalysis MODULE EXCLUDE_FROM_ALL src/Observables/MGAnalysis.cpp)
add_library_no_prefix(MicrogelElasticity MODULE EXCLUDE_FROM_ALL src/Observables/MicrogelElasticity.cpp)
add_library_no_prefix(ConstructwisePressure MODULE EXCLUDE_FROM_ALL src/Observables/ConstructwisePressure.cpp src/Interactions/PolymerSwapInteraction.cpp)
add_library_no_prefix(DensityPressureProfile MODULE EXCLUDE_FROM_ALL src/Observables/DensityPressureProfile.cpp)
add_library_no_prefix(MGAssemblyConf MODULE EXCLUDE_FROM_ALL src/Observables/MGAssemblyConf.cpp src/Interactions/FSInteraction.cpp)

# Interactions
add_library_no_prefix(TSPInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(PolymerInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/PolymerInteraction.cpp src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(FSInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/FSInteraction.cpp)
add_library_no_prefix(GraftedInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/GraftedInteraction.cpp src/Interactions/TSPInteraction.cpp)
add_library_no_prefix(mWInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/mWInteraction.cpp)
add_library_no_prefix(StarrInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/StarrInteraction.cpp)
add_library_no_prefix(CPMixtureInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CPMixtureInteraction.cpp)
add_library_no_prefix(LevyInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/LevyInteraction.cpp)
add_library_no_prefix(MGInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/MGInteraction.cpp)
add_library_no_prefix(AOInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/AOInteraction.cpp)
add_library_no_prefix(PolydisperseLTInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/PolydisperseLTInteraction.cpp)
add_library_no_prefix(PolymerSwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/PolymerSwapInteraction.cpp)
add_library_no_prefix(PatchySwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/PatchySwapInteraction.cpp)
add_library_no_prefix(DetailedPatchySwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/DetailedPatchySwapInteraction.cpp)
add_library_no_prefix(ValerioInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/ValerioInteraction.cpp)

FIND_PACKAGE(PkgConfig)
IF(PkgConfig_FOUND)
	PKG_CHECK_MODULES(GSL gsl)
	IF(GSL_FOUND)
		add_library_no_prefix(NathanNeighs MODULE EXCLUDE_FROM_ALL src/Observables/NathanNeighs.cpp)
		TARGET_LINK_LIBRARIES(NathanNeighs ${GSL_LIBRARIES})
	ENDIF(GSL_FOUND)
ENDIF(PkgConfig_FOUND)

IF(CUDA)
	SET(CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -rdynamic")
	
	# same thing but for CUDA libs
	function(cuda_add_library_no_prefix)
		cuda_add_library(${ARGN})
		set_target_properties(${ARGV0} PROPERTIES PREFIX "")
	endfunction()

	cuda_add_library_no_prefix(CUDADetailedPatchySwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDADetailedPatchySwapInteraction.cu src/Interactions/DetailedPatchySwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPatchySwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAPatchySwapInteraction.cu src/Interactions/PatchySwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDAFSInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAFSInteraction.cu src/Interactions/FSInteraction.cpp)
	cuda_add_library_no_prefix(CUDAmWInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAmWInteraction.cu src/Interactions/mWInteraction.cpp)
	cuda_add_library_no_prefix(CUDAStarrInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAStarrInteraction.cu src/Interactions/StarrInteraction.cpp)
	cuda_add_library_no_prefix(CUDALevyInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDALevyInteraction.cu src/Interactions/LevyInteraction.cpp)
	cuda_add_library_no_prefix(CUDACPMixtureInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDACPMixtureInteraction.cu src/Interactions/CPMixtureInteraction.cpp)
	cuda_add_library_no_prefix(CUDAMGInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAMGInteraction.cu src/Interactions/MGInteraction.cpp)
	cuda_add_library_no_prefix(CUDAAOInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAAOInteraction.cu src/Interactions/AOInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPolymerInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAPolymerInteraction.cu src/Interactions/PolymerInteraction.cpp)
	cuda_add_library_no_prefix(CUDAPolymerSwapInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAPolymerSwapInteraction.cu src/Interactions/PolymerSwapInteraction.cpp)
	cuda_add_library_no_prefix(CUDATSPInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDATSPInteraction.cu src/Interactions/TSPInteraction.cpp)
	
	ADD_DEPENDENCIES(rovigatti CUDADetailedPatchySwapInteraction CUDAPatchySwapInteraction CUDAPolymerInteraction CUDAAOInteraction CUDAMGInteraction CUDACPMixtureInteraction CUDAFSInteraction CUDAmWInteraction CUDAStarrInteraction CUDALevyInteraction CUDAPolymerSwapInteraction CUDATSPInteraction)
ENDIF(CUDA)

set_target_properties(Bonds PROPERTIES PREFIX "")
