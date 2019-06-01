IF(CMAKE_BUILD_TYPE MATCHES Debug)
	SET(exe_name oxDNA_debug)
	SET(lib_name common_debug)
ELSE()
	SET(exe_name oxDNA)
	SET(lib_name common)
ENDIF()

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
	Observables/Configurations/TEPtclOutput.cpp
	Observables/Configurations/TEPxyzOutput.cpp
	Observables/Configurations/JordanOutput.cpp
	Observables/ObservableFactory.cpp
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
)

SET(backends_SOURCES
	Backends/BackendFactory.cpp
	Backends/MDBackend.cpp
	Backends/MCBackend.cpp 
	Backends/MinBackend.cpp
	Backends/FIREBackend.cpp
	Backends/FIREBackend.h
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
	Backends/Thermostats/BussiThermostat.cpp
	Backends/Thermostats/LangevinThermostat.cpp
	Backends/Thermostats/SRDThermostat.cpp
	Backends/Thermostats/DPDThermostat.cpp
)

SET(lists_SOURCES
	Lists/NoList.cpp
	Lists/Cells.cpp
	Lists/RodCells.cpp
	Lists/VerletList.cpp
	Lists/BinVerletList.cpp
	Lists/ListFactory.cpp
)

SET(interactions_SOURCES
	Interactions/JordanInteraction.cpp
	Interactions/TEPInteraction.cpp
	Interactions/TSPInteraction.cpp
	Interactions/LJInteraction.cpp
	Interactions/HSInteraction.cpp
	Interactions/DHSInteraction.cpp
	Interactions/HardCylinderInteraction.cpp
	Interactions/HardSpheroCylinderInteraction.cpp
	Interactions/BoxInteraction.cpp
	Interactions/PatchyInteraction.cpp
	Interactions/PatchyInteractionDan.cpp
	Interactions/KFInteraction.cpp
	Interactions/RNAInteraction.cpp
	Interactions/RNAInteraction_relax.cpp
	Interactions/DNAInteraction.cpp
	Interactions/RNAInteraction2.cpp
	Interactions/DNA2Interaction.cpp
	Interactions/DNAInteraction_nomesh.cpp
	Interactions/DNAInteraction_relax.cpp
	Interactions/DNAInteraction_relax2.cpp
	Interactions/CustomInteraction.cpp
	Interactions/InteractionFactory.cpp
)

SET(move_SOURCES
	Backends/MCMoves/MoveFactory.cpp
	Backends/MCMoves/MCTras.cpp
	Backends/MCMoves/MCRot.cpp
	Backends/MCMoves/VMMC.cpp
	Backends/MCMoves/VolumeMove.cpp
	Backends/MCMoves/ShapeMove.cpp
	Backends/MCMoves/RotateSite.cpp
)

SET(forces_SOURCES
	Forces/COMForce.cpp
	Forces/ConstantRateForce.cpp
	Forces/SawtoothForce.cpp
	Forces/ConstantRateTorque.cpp
	Forces/ConstantTrap.cpp
	Forces/MovingTrap.cpp
	Forces/LowdimMovingTrap.cpp
	Forces/MutualTrap.cpp
	Forces/RepulsionPlane.cpp
	Forces/RepulsionPlaneMoving.cpp
	Forces/RepulsiveSphere.cpp
	Forces/RepulsiveSphereSmooth.cpp
	Forces/LJWall.cpp
	Forces/HardWall.cpp
	Forces/AlignmentField.cpp
	Forces/GenericCentralForce.cpp
	Forces/LJCone.cpp
	Forces/ForceFactory.cpp
	Forces/BaseForce.cpp
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
	Particles/TSPParticle.cpp
	Particles/SpheroCylinder.cpp
	Particles/CustomParticle.cpp
	Managers/SimManager.cpp
	Utilities/OrderParameters.cpp
	Utilities/Weights.cpp
	Utilities/Histogram.cpp
	Utilities/Utils.cpp
	Utilities/oxDNAException.cpp
	Utilities/Logger.cpp
	Utilities/parse_input/parse_input.cpp
	Utilities/time_scales/time_scales.c
	Utilities/SignalManager.cpp
	Utilities/ConfigInfo.cpp
	PluginManagement/PluginManager.cpp
	${forces_SOURCES}
	${observables_SOURCES}
	${lists_SOURCES}
	${backends_SOURCES}
	${interactions_SOURCES}
	${move_SOURCES}
	${box_SOURCES}
)

SET(oxDNA_SOURCES
	oxDNA.cpp 
)

SET(DNAnalysis_SOURCES
	DNAnalysis.cpp
	Managers/AnalysisManager.cpp
	Backends/AnalysisBackend.cpp
)

SET(confGenerator_SOURCES
	confGenerator.cpp
	Managers/GeneratorManager.cpp
)
	
IF(CUDA)
	SET(oxDNA_CUDASOURCES
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
		CUDA/Interactions/CUDATEPInteraction.cu
		CUDA/Interactions/CUDAInteractionFactory.cu
		CUDA/Lists/CUDANoList.cu
		CUDA/Lists/CUDASimpleVerletList.cu
		CUDA/Lists/CUDABinVerletList.cu
		CUDA/Lists/CUDAListFactory.cu
		CUDA/Thermostats/CUDABaseThermostat.cu
		CUDA/Thermostats/CUDANoThermostat.cu
		CUDA/Thermostats/CUDAThermostatFactory.cu
		CUDA/Thermostats/CUDABrownianThermostat.cu
		CUDA/Thermostats/CUDASRDThermostat.cu
		CUDA/Thermostats/CUDALangevinThermostat.cu
		CUDA/Thermostats/CUDABussiThermostat.cu
		Utilities/Timings.cpp
	)
	
	CUDA_ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES} ${oxDNA_CUDASOURCES})
	CUDA_ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES} Utilities/Timings.cpp)
ELSE()
	SET(common_SOURCES
		${common_SOURCES}
		Utilities/Timings.cpp
	)

	ADD_EXECUTABLE(${exe_name} ${oxDNA_SOURCES})
	ADD_EXECUTABLE(DNAnalysis ${DNAnalysis_SOURCES})
ENDIF(CUDA)

ADD_EXECUTABLE(confGenerator ${confGenerator_SOURCES})

ADD_LIBRARY(${lib_name} ${common_SOURCES})
TARGET_LINK_LIBRARIES(${lib_name} ${CMAKE_DL_LIBS})
TARGET_LINK_LIBRARIES(${exe_name} ${lib_name})
TARGET_LINK_LIBRARIES(DNAnalysis ${lib_name})
TARGET_LINK_LIBRARIES(confGenerator ${lib_name})

# we add these executable as dependencies for the test targets
ADD_DEPENDENCIES(test_run ${exe_name} DNAnalysis confGenerator)
ADD_DEPENDENCIES(test_quick ${exe_name} DNAnalysis confGenerator)
ADD_DEPENDENCIES(test_scientific ${exe_name} DNAnalysis confGenerator)

IF(MPI)
	FIND_PACKAGE(MPI REQUIRED)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

	SET(mpi_SOURCES
		oxDNA_mpi
		Managers/ParallelManager.cpp
                Backends/PT_VMMC_CPUBackend.cpp
	)
	ADD_EXECUTABLE(oxDNA_mpi ${mpi_SOURCES})
        SET(exe_name oxDNA_mpi)
	TARGET_LINK_LIBRARIES(oxDNA_mpi ${lib_name} ${MPI_LIBRARIES})
ENDIF()
