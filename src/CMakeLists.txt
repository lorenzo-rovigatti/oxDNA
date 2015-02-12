IF(CMAKE_BUILD_TYPE MATCHES Debug)
	SET(exe_name oxDNA_debug)
	SET(lib_name common_debug)
ELSE()
	SET(exe_name oxDNA)
	SET(lib_name common)
ENDIF()

IF(MPI)
	FIND_PACKAGE(MPI REQUIRED)
	INCLUDE(CMakeForceCompiler)
	CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")

	IF(CMAKE_BUILD_TYPE MATCHES Debug)
		SET(exe_name oxDNA_MPI_debug)
		SET(lib_name common_MPI_debug)
	ELSE()
		SET(exe_name oxDNA_MPI)
		SET(lib_name common_MPI)
	ENDIF()
ENDIF(MPI)

SET(observables_SOURCES
	Observables/ObservableOutput.cpp
	Observables/Step.cpp
	Observables/PotentialEnergy.cpp
	Observables/KineticEnergy.cpp
	Observables/TotalEnergy.cpp
	Observables/HBEnergy.cpp
	Observables/OrderParameterValues.cpp
	Observables/HBList.cpp
	Observables/StrandwiseBonds.cpp
	Observables/PairEnergy.cpp
	Observables/PairForce.cpp
	Observables/BackendInfo.cpp
	Observables/ForceEnergy.cpp
	Observables/Pressure.cpp
	Observables/Density.cpp
	Observables/Rdf.cpp
	Observables/Distance.cpp
	Observables/DensityProfile.cpp
	Observables/ParticlePosition.cpp
	Observables/CoaxVariables.cpp
	Observables/Pitch.cpp
	Observables/SaltExtrapolation.cpp
	Observables/ExternalTorque.cpp
	Observables/Configurations/Configuration.cpp
	Observables/Configurations/BinaryConfiguration.cpp
	Observables/Configurations/TclOutput.cpp
	Observables/Configurations/PdbOutput.cpp
	Observables/Configurations/ChimeraOutput.cpp
	Observables/ObservableFactory.cpp
)

SET(backends_SOURCES
	Backends/BackendFactory.cpp
	Backends/MDBackend.cpp
	Backends/MCBackend.cpp 
	Backends/MD_CPUBackend.cpp 
	Backends/MC_CPUBackend.cpp
	Backends/MC_CPUBackend2.cpp
	Backends/SimBackend.cpp 
	Backends/FFS_MD_CPUBackend.cpp
	Backends/VMMC_CPUBackend.cpp
	Backends/Thermostats/ThermostatFactory.cpp
	Backends/Thermostats/BrownianThermostat.cpp
	Backends/Thermostats/NoThermostat.cpp
	Backends/Thermostats/RefreshThermostat.cpp
	Backends/Thermostats/LangevinThermostat.cpp
	Backends/Thermostats/SRDThermostat.cpp
)

SET(lists_SOURCES
	Lists/NoList.cpp
	Lists/Cells.cpp
	Lists/VerletList.cpp
	Lists/ListFactory.cpp
)

SET(interactions_SOURCES
	Interactions/TSPInteraction.cpp
	Interactions/LJInteraction.cpp
	Interactions/HSInteraction.cpp
	Interactions/DHSInteraction.cpp
	Interactions/DirkInteraction.cpp
	Interactions/DirkInteraction2.cpp
	Interactions/DirkInteractionBias.cpp
	Interactions/HardCylinderInteraction.cpp
	Interactions/HardSpheroCylinderInteraction.cpp
	Interactions/BoxInteraction.cpp
	Interactions/PatchyInteraction.cpp
	Interactions/RNAInteraction.cpp
	Interactions/RNAInteraction_relax.cpp
	Interactions/DNAInteraction.cpp
	Interactions/RNAInteraction2.cpp
	Interactions/DNA2Interaction.cpp
	Interactions/DNAInteraction_nomesh.cpp
	Interactions/DNAInteraction_relax.cpp
	Interactions/CustomInteraction.cpp
	Interactions/InteractionFactory.cpp
)

SET(move_SOURCES
	Backends/MCMoves/MoveFactory.cpp
	Backends/MCMoves/MCTras.cpp
	Backends/MCMoves/MCRot.cpp
)

SET(forces_SOURCES
	Forces/COMForce.cpp
	Forces/ConstantRateForce.cpp
	Forces/ConstantRateTorque.cpp
	Forces/ConstantTrap.cpp
	Forces/MovingTrap.cpp
	Forces/LowdimMovingTrap.cpp
	Forces/MutualTrap.cpp
	Forces/RepulsionPlane.cpp
	Forces/RepulsionPlaneMoving.cpp
	Forces/RepulsiveSphere.cpp
	Forces/ForceFactory.cpp
)

SET(common_SOURCES
	Particles/RNANucleotide.cpp
	Particles/BaseParticle.cpp
	Particles/DNANucleotide.cpp
	Particles/PatchyParticle.cpp
	Particles/TSPParticle.cpp
	Particles/CustomParticle.cpp
	Managers/SimManager.cpp
	Utilities/OrderParameters.cpp
	Utilities/Weights.cpp
	Utilities/Histogram.cpp
	Utilities/Utils.cpp
	Utilities/oxDNAException.cpp
	Utilities/Logger.cpp
	Utilities/parse_input/parse_input.cpp
	Utilities/timing/timing.c
	Utilities/time_scales/time_scales.c
	PluginManagement/PluginManager.cpp
	${forces_SOURCES}
	${observables_SOURCES}
	${lists_SOURCES}
	${backends_SOURCES}
	${interactions_SOURCES}
	${move_SOURCES}
)

ADD_LIBRARY(${lib_name} ${common_SOURCES})
TARGET_LINK_LIBRARIES(${lib_name} ${CMAKE_DL_LIBS})

SET(oxDNA_SOURCES
	oxDNA.cpp 
)
	
IF (MPI) 
SET(oxDNA_SOURCES ${oxDNA_SOURCES} 
		Backends/MD_MPIBackend.cpp
		Backends/PT_VMMC_CPUBackend.cpp
	)
ENDIF (MPI)


IF(CUDA)
	SET(oxDNA_CUDASOURCES
		CUDA/CUDA_rand.cu
		CUDA/CUDA_sort.cu
		CUDA/CUDAUtils.cu
		CUDA/cuda_utils/cuda_device_utils.cu
		CUDA/Backends/CUDABaseBackend.cu
		CUDA/Backends/MD_CUDABackend.cu
		CUDA/Backends/MD_CUDAMixedBackend.cu
		CUDA/Backends/FFS_MD_CUDAMixedBackend.cu
		CUDA/Interactions/CUDABaseInteraction.cu
		CUDA/Interactions/CUDADNAInteraction.cu
		CUDA/Interactions/CUDALJInteraction.cu
		CUDA/Interactions/CUDARNAInteraction.cu
		CUDA/Interactions/CUDAPatchyInteraction.cu
		CUDA/Interactions/CUDATSPInteraction.cu
		CUDA/Interactions/CUDAInteractionFactory.cu
		CUDA/Lists/CUDANoList.cu
		CUDA/Lists/CUDASimpleVerletList.cu
		CUDA/Lists/CUDASharedVerletList.cu
		CUDA/Lists/CUDAListFactory.cu
		CUDA/Thermostats/CUDABaseThermostat.cu
		CUDA/Thermostats/CUDANoThermostat.cu
		CUDA/Thermostats/CUDAThermostatFactory.cu
		CUDA/Thermostats/CUDABrownianThermostat.cu
		CUDA/Thermostats/CUDASRDThermostat.cu
		CUDA/Thermostats/CUDALangevinThermostat.cu
	)
	CUDA_ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES} ${oxDNA_CUDASOURCES})
ELSE()
	ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES})
ENDIF(CUDA)

TARGET_LINK_LIBRARIES(${exe_name} ${lib_name})

SET(DNAnalysis_SOURCES
	DNAnalysis.cpp
	Managers/AnalysisManager.cpp
	Backends/AnalysisBackend.cpp
)

ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES})
TARGET_LINK_LIBRARIES(DNAnalysis ${lib_name})

SET(confGenerator_SOURCES
	confGenerator.cpp
	Managers/GeneratorManager.cpp
)

ADD_EXECUTABLE(confGenerator ${confGenerator_SOURCES})
TARGET_LINK_LIBRARIES(confGenerator ${lib_name})
