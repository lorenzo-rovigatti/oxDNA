set(CMAKE_POSITION_INDEPENDENT_CODE ON)

SET(lib_name oxdna_common)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
	SET(exe_name oxDNA_debug)
ELSEIF(MPI)
	SET(exe_name oxDNA_mpi)
ELSE()
	SET(exe_name oxDNA)
ENDIF()

SET(observables_SOURCES
	Observables/BaseObservable.cpp
	Observables/ObservableFactory.cpp
	Observables/ObservableOutput.cpp
	Observables/Step.cpp
	Observables/PotentialEnergy.cpp
	Observables/KineticEnergy.cpp
	Observables/TotalEnergy.cpp
	Observables/BackendInfo.cpp	
	Observables/HBEnergy.cpp
	Observables/OrderParameterValues.cpp
	Observables/HBList.cpp
	Observables/StrandwiseBonds.cpp
	Observables/PairEnergy.cpp
	Observables/PairForce.cpp
	Observables/ForceEnergy.cpp
	Observables/Pressure.cpp
	Observables/StressAutocorrelation.cpp
	Observables/Density.cpp
	Observables/Rdf.cpp
	Observables/Distance.cpp
	Observables/DensityProfile.cpp
	Observables/ParticlePosition.cpp
	Observables/CoaxVariables.cpp
	Observables/Pitch.cpp
	Observables/SaltExtrapolation.cpp
	Observables/ExternalTorque.cpp
	Observables/MeanVectorCosine.cpp
	Observables/VectorAngle.cpp
	Observables/Checkpoint.cpp
	Observables/Contacts.cpp
	Observables/Writhe.cpp
	Observables/UnstackedList.cpp
	Observables/PlectonemePosition.cpp
	Observables/StretchedBonds.cpp
	Observables/StructureFactor.cpp
	Observables/FormFactor.cpp
	Observables/UnstackedList.cpp
	Observables/PlectonemePosition.cpp
	Observables/TEPPlectonemePosition.cpp
	Observables/AverageEnergy.cpp
	Observables/ContactMap.cpp
	Observables/AllVectors.cpp
	Observables/ExternalForce.cpp
	Observables/Configurations/Configuration.cpp
	Observables/Configurations/BinaryConfiguration.cpp
	Observables/Configurations/TclOutput.cpp
	Observables/Configurations/PdbOutput.cpp
	Observables/Configurations/ChimeraOutput.cpp
	Observables/Configurations/TEPtclOutput.cpp
	Observables/Configurations/TEPxyzOutput.cpp
	Observables/Configurations/JordanOutput.cpp
)

SET(backends_SOURCES
	Backends/BackendFactory.cpp
	Backends/SimBackend.cpp
	Backends/AnalysisBackend.cpp 
	Backends/MDBackend.cpp
	Backends/MCBackend.cpp
	Backends/MD_CPUBackend.cpp 
	Backends/MC_CPUBackend.cpp
	Backends/MC_CPUBackend2.cpp
	Backends/VMMC_CPUBackend.cpp
	Backends/MinBackend.cpp
	Backends/FIREBackend.cpp
	Backends/FIREBackend.h
	Backends/FFS_MD_CPUBackend.cpp
	Backends/Thermostats/SRDThermostat.cpp
	Backends/Thermostats/DPDThermostat.cpp
	Backends/Thermostats/ThermostatFactory.cpp
	Backends/Thermostats/BaseThermostat.cpp
	Backends/Thermostats/NoThermostat.cpp
	Backends/Thermostats/BrownianThermostat.cpp
	Backends/Thermostats/RefreshThermostat.cpp
	Backends/Thermostats/BussiThermostat.cpp
	Backends/Thermostats/LangevinThermostat.cpp
)

SET(lists_SOURCES
	Lists/BaseList.cpp
	Lists/NoList.cpp
	Lists/Cells.cpp
	Lists/RodCells.cpp
	Lists/VerletList.cpp
	Lists/BinVerletList.cpp
	Lists/ListFactory.cpp
)

SET(interactions_SOURCES
	Interactions/InteractionFactory.cpp
	Interactions/InteractionUtils.cpp
	Interactions/Mesh.cpp
	Interactions/BaseInteraction.cpp
	Interactions/LJInteraction.cpp
	Interactions/DNAInteraction.cpp
	Interactions/DNA2Interaction.cpp
	Interactions/DNAInteraction_nomesh.cpp
	Interactions/DNAInteraction_relax.cpp
	Interactions/RNAInteraction.cpp
	Interactions/RNAInteraction2.cpp
	Interactions/RNAInteraction_relax.cpp
	Interactions/JordanInteraction.cpp
	Interactions/TEPInteraction.cpp
	Interactions/HSInteraction.cpp
	Interactions/DHSInteraction.cpp
	Interactions/HardCylinderInteraction.cpp
	Interactions/HardSpheroCylinderInteraction.cpp
	Interactions/BoxInteraction.cpp
	Interactions/PatchyInteraction.cpp
	Interactions/PatchyInteractionDan.cpp
	Interactions/KFInteraction.cpp
	Interactions/CustomInteraction.cpp
	Interactions/DRHInteraction.cpp
	Interactions/DRHInteraction_relax.cpp
)

SET(move_SOURCES
	Backends/MCMoves/BaseMove.cpp
	Backends/MCMoves/MoveFactory.cpp
	Backends/MCMoves/MCTras.cpp
	Backends/MCMoves/MCRot.cpp
	Backends/MCMoves/Pivot.cpp
	Backends/MCMoves/VMMC.cpp
	Backends/MCMoves/VolumeMove.cpp
	Backends/MCMoves/MoleculeVolumeMove.cpp
	Backends/MCMoves/ShapeMove.cpp
	Backends/MCMoves/RotateSite.cpp
)

SET(forces_SOURCES
	Forces/ForceFactory.cpp
	Forces/BaseForce.cpp
	Forces/COMForce.cpp
	Forces/ConstantRateForce.cpp
	Forces/SawtoothForce.cpp
	Forces/ConstantRateTorque.cpp
	Forces/ConstantTrap.cpp
	Forces/MovingTrap.cpp
	Forces/LowdimMovingTrap.cpp
	Forces/MutualTrap.cpp
	Forces/RepulsiveEllipsoid.cpp
	Forces/RepulsionPlane.cpp
	Forces/RepulsionPlaneMoving.cpp
	Forces/RepulsiveSphere.cpp
	Forces/RepulsiveSphereSmooth.cpp
	Forces/LJWall.cpp
	Forces/HardWall.cpp
	Forces/AlignmentField.cpp
	Forces/GenericCentralForce.cpp
	Forces/LJCone.cpp
	Forces/Metadynamics/LTCOMTrap.cpp
	Forces/Metadynamics/LTCOMAngleTrap.cpp
	Forces/Metadynamics/LT2DCOMTrap.cpp
	Forces/Metadynamics/LTAtanCOMTrap.cpp
	Forces/Metadynamics/meta_utils.cpp
)

SET(box_SOURCES
	Boxes/BaseBox.cpp
	Boxes/CubicBox.cpp
	Boxes/OrthogonalBox.cpp
	Boxes/LeesEdwardsCubicBox.cpp
	Boxes/BoxFactory.cpp
)

SET(common_SOURCES
	Particles/JordanParticle.cpp
	Particles/TEPParticle.cpp
	Particles/RNANucleotide.cpp
	Particles/BaseParticle.cpp
	Particles/DNANucleotide.cpp
	Particles/PatchyParticle.cpp
	Particles/PatchyParticleDan.cpp
	Particles/SpheroCylinder.cpp
	Particles/CustomParticle.cpp
	Particles/Molecule.cpp
	Managers/SimManager.cpp
	Utilities/OrderParameters.cpp
	Utilities/Weights.cpp
	Utilities/Histogram.cpp
	Utilities/LR_vector.cpp
	Utilities/LR_matrix.cpp
	Utilities/Utils.cpp
	Utilities/oxDNAException.cpp
	Utilities/Logger.cpp
	Utilities/parse_input/parse_input.cpp
	Utilities/time_scales/time_scales.cpp
	Utilities/SignalManager.cpp
	Utilities/ConfigInfo.cpp
	Utilities/FlattenedConfigInfo.cpp
	Utilities/TopologyParser.cpp
	PluginManagement/PluginManager.cpp
	${forces_SOURCES}
	${observables_SOURCES}
	${lists_SOURCES}
	${backends_SOURCES}
	${interactions_SOURCES}
	${move_SOURCES}
	${box_SOURCES}
)

SET(DNAnalysis_SOURCES
	DNAnalysis.cpp
	Managers/AnalysisManager.cpp
)

SET(confGenerator_SOURCES
	confGenerator.cpp
	Managers/GeneratorManager.cpp
)

IF(MPI)
	FIND_PACKAGE(MPI REQUIRED)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	
	LIST(APPEND common_SOURCES
		Backends/PT_VMMC_CPUBackend.cpp
		Managers/ParallelManager.cpp
	)
	
	SET(oxDNA_SOURCES
		oxDNA_mpi.cpp 
	)
ELSE()
	SET(oxDNA_SOURCES
		oxDNA.cpp 
	)
ENDIF()
	
IF(CUDA)
	CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
	FIND_PACKAGE("CUDA" REQUIRED)
	
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

	INCLUDE_DIRECTORIES(${CUDA_INSTALL_PREFIX}/include/)
	LINK_DIRECTORIES(${CUDA_INSTALL_PREFIX}/lib)
	
	INCLUDE(../select_compute_arch.cmake)
	
	IF(CUDA_COMMON_ARCH)
		CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Common)
	ELSE(CUDA_COMMON_ARCH)
		CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
	ENDIF(CUDA_COMMON_ARCH)
	
	MESSAGE(STATUS "Generating code for the following CUDA architectures: ${ARCH_FLAGS_readable}")
	LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
	MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")
   
	IF(G)
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lineinfo)
	ENDIF(G)

	IF(DEBUG)
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lineinfo -G)
	ENDIF(DEBUG)
	
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -use_fast_math --compiler-options "-ffast-math -O3" -Xcudafe "--diag_suppress=unrecognized_gcc_pragma")

	LIST(APPEND common_SOURCES
		CUDA/CUDA_sort.cu
		CUDA/CUDAUtils.cu
		CUDA/cuda_utils/cuda_device_utils.cu
		CUDA/Backends/CUDABaseBackend.cu
		CUDA/Backends/MD_CUDABackend.cu
		CUDA/Interactions/CUDABaseInteraction.cu
		CUDA/Interactions/CUDADNAInteraction.cu
		CUDA/Interactions/CUDALJInteraction.cu
		CUDA/Interactions/CUDARNAInteraction.cu
		CUDA/Interactions/CUDAPatchyInteraction.cu
		CUDA/Interactions/CUDATEPInteraction.cu
		CUDA/Interactions/CUDAInteractionFactory.cu
		CUDA/Lists/CUDAListFactory.cu
		CUDA/Lists/CUDANoList.cu
		CUDA/Lists/CUDASimpleVerletList.cu
		CUDA/Lists/CUDABinVerletList.cu
		CUDA/Thermostats/CUDABaseThermostat.cu
		CUDA/Thermostats/CUDANoThermostat.cu
		CUDA/Thermostats/CUDAThermostatFactory.cu
		CUDA/Thermostats/CUDABrownianThermostat.cu
		CUDA/Thermostats/CUDASRDThermostat.cu
		CUDA/Thermostats/CUDALangevinThermostat.cu
		CUDA/Thermostats/CUDABussiThermostat.cu
		Utilities/Timings.cpp
	)
	
	IF(CUDA_DOUBLE)
		ADD_DEFINITIONS(-DCUDA_DOUBLE_PRECISION)
	ELSE()
		LIST(APPEND common_SOURCES
			CUDA/Backends/MD_CUDAMixedBackend.cu 
			CUDA/Backends/FFS_MD_CUDAMixedBackend.cu
		)
	ENDIF()
	
	# we need the library to be shared or the oxpy compilation fails
	CUDA_ADD_LIBRARY(${lib_name} SHARED ${common_SOURCES})
	CUDA_ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES})
	CUDA_ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES} Utilities/Timings.cpp)
	
	SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
ELSE()
	ADD_DEFINITIONS(-DNOCUDA)
	SET(common_SOURCES
		${common_SOURCES}
		Utilities/Timings.cpp
	)

	ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES})
	ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES})
	ADD_LIBRARY(${lib_name} ${common_SOURCES})
ENDIF(CUDA)

ADD_EXECUTABLE(confGenerator ${confGenerator_SOURCES})

IF(MPI)
	TARGET_LINK_LIBRARIES(${lib_name} ${CMAKE_DL_LIBS} ${MPI_LIBRARIES})
ELSE()
	TARGET_LINK_LIBRARIES(${lib_name} ${CMAKE_DL_LIBS})
ENDIF()

TARGET_LINK_LIBRARIES(${exe_name} ${lib_name})
TARGET_LINK_LIBRARIES(DNAnalysis ${lib_name})
TARGET_LINK_LIBRARIES(confGenerator ${lib_name})

# we add these executable as dependencies for the test targets
ADD_DEPENDENCIES(test_run ${exe_name} DNAnalysis confGenerator)
ADD_DEPENDENCIES(test_quick ${exe_name} DNAnalysis confGenerator)
