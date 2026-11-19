SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(tostiguerra
	DEPENDS ANNaMoInteraction
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

# Interactions
add_library_no_prefix(ANNaMoInteraction MODULE EXCLUDE_FROM_ALL src/ANNaMoInteraction.cpp)

if(CUDA)
	set(CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -rdynamic")

	# same thing but for CUDA libs. With CMake >= 3.18 the main build compiles CUDA as a first-class language, and the
	# plugin inherits its architectures and C++ standard (C++17 with CUDA 13, whose Thrust/CUB require it); the legacy
	# FindCUDA path (cuda_add_library) runs nvcc on its own and gets neither. Same as contrib/rovigatti
	function(cuda_add_library_no_prefix target source)
		if(CMAKE_VERSION VERSION_LESS 3.18)
			cuda_add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
			target_link_libraries(${target} ${CUDA_LIBRARIES})
		else()
			add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
			target_link_libraries(${target} CUDA::cudart)
		endif()
		set_target_properties(${target} PROPERTIES PREFIX "")
	endfunction()

	# the CUDA factory loads interaction_type = X from the plugin CUDAX
	cuda_add_library_no_prefix(CUDAANNaMoInteraction src/CUDANNaMoInteraction.cu src/ANNaMoInteraction.cpp)
	ADD_DEPENDENCIES(tostiguerra CUDAANNaMoInteraction)
endif(CUDA)

# smoke test: make test_annamo builds oxDNA, confGenerator and the plugin(s) and runs test/run_smoke_test.py
# (a short simulation on the CPU and, in CUDA builds with a GPU, on the GPU). As for test_run and test_quick, the
# script is run directly (it needs python3 with its standard library only), so CMake does not look for Python
IF(CMAKE_BUILD_TYPE MATCHES Debug)
	SET(annamo_oxdna_target oxDNA_debug)
ELSEIF(MPI)
	SET(annamo_oxdna_target oxDNA_mpi)
ELSE()
	SET(annamo_oxdna_target oxDNA)
ENDIF()
set(annamo_test_args
	--oxdna $<TARGET_FILE:${annamo_oxdna_target}>
	--confgen $<TARGET_FILE:confGenerator>
	--plugins ${CMAKE_CURRENT_SOURCE_DIR}
	--work ${CMAKE_CURRENT_BINARY_DIR}/test_annamo
)
if(CUDA)
	if(CUDA_DOUBLE)
		list(APPEND annamo_test_args --cuda double)
	else()
		list(APPEND annamo_test_args --cuda mixed)
	endif()
endif()
add_custom_target(test_annamo
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test/run_smoke_test.py ${annamo_test_args}
	COMMENT "Running the ANNaMo smoke test" VERBATIM
)
add_dependencies(test_annamo ${annamo_oxdna_target} confGenerator ANNaMoInteraction)
if(CUDA)
	add_dependencies(test_annamo CUDAANNaMoInteraction)
endif()
