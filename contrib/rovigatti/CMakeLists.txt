SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(rovigatti
	DEPENDS PolymerSwapInteraction TSPInteraction MGAssemblyConf PolymerInteraction 
			DensityPressureProfile PolydisperseLTInteraction ConstructwisePressure 
			MicrogelElasticity MGAnalysis AOInteraction RadialDensityProfile MGInteraction 
			GenericGrByInsertion VoidPercolation Remoteness CPAnalysis Widom LevyDelta 
			LevyInteraction CPMixtureInteraction StarrInteraction mWInteraction 
			GraftedInteraction FSInteraction TSPAnalysis StarrConf ConstructwiseBonds 
			EmptyVolume FSConf Gyradius PatchySwapInteraction DetailedPatchySwapInteraction 
			PatchyBonds
)

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

# Observables
ADD_LIBRARY(PatchyBonds SHARED EXCLUDE_FROM_ALL src/Observables/PatchyBonds.cpp)
ADD_LIBRARY(ConstructwiseBonds SHARED EXCLUDE_FROM_ALL src/Observables/ConstructwiseBonds.cpp)
ADD_LIBRARY(EmptyVolume SHARED EXCLUDE_FROM_ALL src/Observables/EmptyVolume.cpp)
ADD_LIBRARY(FSConf SHARED EXCLUDE_FROM_ALL src/Observables/FSConf.cpp)
TARGET_LINK_LIBRARIES(FSConf FSInteraction)
ADD_LIBRARY(Gyradius SHARED EXCLUDE_FROM_ALL src/Observables/Gyradius.cpp)
ADD_LIBRARY(StarrConf SHARED EXCLUDE_FROM_ALL src/Observables/StarrConf.cpp)
ADD_LIBRARY(TSPAnalysis SHARED EXCLUDE_FROM_ALL src/Observables/TSPAnalysis.cpp)
ADD_LIBRARY(Widom SHARED EXCLUDE_FROM_ALL src/Observables/Widom.cpp)
ADD_LIBRARY(LevyDelta SHARED EXCLUDE_FROM_ALL src/Observables/LevyDelta.cpp)
ADD_LIBRARY(CPAnalysis SHARED EXCLUDE_FROM_ALL src/Observables/CPAnalysis.cpp)
ADD_LIBRARY(Remoteness SHARED EXCLUDE_FROM_ALL src/Observables/Remoteness.cpp)
ADD_LIBRARY(VoidPercolation SHARED EXCLUDE_FROM_ALL src/Observables/VoidPercolation.cpp)
ADD_LIBRARY(GenericGrByInsertion SHARED EXCLUDE_FROM_ALL src/Observables/GenericGrByInsertion.cpp)
ADD_LIBRARY(RadialDensityProfile SHARED EXCLUDE_FROM_ALL src/Observables/RadialDensityProfile.cpp)
ADD_LIBRARY(MGAnalysis SHARED EXCLUDE_FROM_ALL src/Observables/MGAnalysis.cpp)
ADD_LIBRARY(MicrogelElasticity SHARED EXCLUDE_FROM_ALL src/Observables/MicrogelElasticity.cpp)
ADD_LIBRARY(ConstructwisePressure SHARED EXCLUDE_FROM_ALL src/Observables/ConstructwisePressure.cpp)
TARGET_LINK_LIBRARIES(ConstructwisePressure PolymerSwapInteraction)
ADD_LIBRARY(DensityPressureProfile SHARED EXCLUDE_FROM_ALL src/Observables/DensityPressureProfile.cpp)
ADD_LIBRARY(MGAssemblyConf SHARED EXCLUDE_FROM_ALL src/Observables/MGAssemblyConf.cpp)
TARGET_LINK_LIBRARIES(MGAssemblyConf FSInteraction)

# Interactions
ADD_LIBRARY(TSPInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/TSPInteraction.cpp)
ADD_LIBRARY(PolymerInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PolymerInteraction.cpp)
ADD_LIBRARY(FSInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/FSInteraction.cpp)
ADD_LIBRARY(GraftedInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/GraftedInteraction.cpp src/Interactions/TSPInteraction.cpp)
ADD_LIBRARY(mWInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/mWInteraction.cpp)
ADD_LIBRARY(StarrInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/StarrInteraction.cpp)
ADD_LIBRARY(CPMixtureInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CPMixtureInteraction.cpp)
ADD_LIBRARY(LevyInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/LevyInteraction.cpp)
ADD_LIBRARY(MGInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/MGInteraction.cpp)
ADD_LIBRARY(AOInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/AOInteraction.cpp)
ADD_LIBRARY(PolydisperseLTInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PolydisperseLTInteraction.cpp)
ADD_LIBRARY(PolymerSwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PolymerSwapInteraction.cpp)
ADD_LIBRARY(PatchySwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/PatchySwapInteraction.cpp)
ADD_LIBRARY(DetailedPatchySwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/DetailedPatchySwapInteraction.cpp)

FIND_PACKAGE(PkgConfig)
IF(PkgConfig_FOUND)
	PKG_CHECK_MODULES(GSL gsl)
	IF(GSL_FOUND)
		ADD_LIBRARY(NathanNeighs SHARED EXCLUDE_FROM_ALL src/Observables/NathanNeighs.cpp)
		TARGET_LINK_LIBRARIES(NathanNeighs ${GSL_LIBRARIES})
	ENDIF(GSL_FOUND)
ENDIF(PkgConfig_FOUND)

IF(CUDA)
	SET(CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -rdynamic")

	CUDA_ADD_LIBRARY(CUDADetailedPatchySwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDADetailedPatchySwapInteraction.cu src/Interactions/DetailedPatchySwapInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAPatchySwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAPatchySwapInteraction.cu src/Interactions/PatchySwapInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAFSInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAFSInteraction.cu src/Interactions/FSInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAmWInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAmWInteraction.cu src/Interactions/mWInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAStarrInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAStarrInteraction.cu src/Interactions/StarrInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDALevyInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDALevyInteraction.cu src/Interactions/LevyInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDACPMixtureInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDACPMixtureInteraction.cu src/Interactions/CPMixtureInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAMGInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAMGInteraction.cu src/Interactions/MGInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAAOInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAAOInteraction.cu src/Interactions/AOInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAPolymerInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAPolymerInteraction.cu src/Interactions/PolymerInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDAPolymerSwapInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDAPolymerSwapInteraction.cu src/Interactions/PolymerSwapInteraction.cpp)
	CUDA_ADD_LIBRARY(CUDATSPInteraction SHARED EXCLUDE_FROM_ALL src/Interactions/CUDATSPInteraction.cu src/Interactions/TSPInteraction.cpp)
	
	ADD_DEPENDENCIES(rovigatti CUDADetailedPatchySwapInteraction CUDAPatchySwapInteraction CUDAPolymerInteraction CUDAAOInteraction CUDAMGInteraction CUDACPMixtureInteraction CUDAFSInteraction CUDAmWInteraction CUDAStarrInteraction CUDALevyInteraction CUDAPolymerSwapInteraction CUDATSPInteraction)
ENDIF(CUDA)
