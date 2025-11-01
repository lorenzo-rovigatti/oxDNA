SET(oxpy_SOURCES
	OxpyContext.cpp
	OxpyManager.cpp
)

IF(Python)
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
	add_subdirectory(pybind11 EXCLUDE_FROM_ALL) # the EXCLUDE_FROM_ALL forces cmake to use the local pybind11 files
	ADD_DEFINITIONS(-DPYTHON_BINDINGS)
	
	include_directories( ${PROJECT_SOURCE_DIR}/oxpy/pybind11/include )
	include_directories( ${PROJECT_SOURCE_DIR}/src )
	
	MESSAGE(STATUS "Enabling Python bindings")
	
	SET(CMAKE_SHARED_LIBRARY_PREFIX "")
	
	# we need the oxDNA "common" library to be compiled so that the resulting code is position-independent and the library can be linked dynamically
	set_target_properties(oxdna_common PROPERTIES POSITION_INDEPENDENT_CODE ON)
	
	Python_add_library(core MODULE bindings.cpp ${oxpy_SOURCES} WITH_SOABI)
	set_target_properties(core PROPERTIES POSITION_INDEPENDENT_CODE ON)
	target_include_directories(core PRIVATE ${Python_INCLUDE_DIRS})
	target_link_libraries(core PUBLIC oxdna_common pybind11::headers ${Python_LIBRARIES})
	
	set(PYTHON_SITE_PACKAGES ${Python_SITEARCH})
	install(TARGETS core
			LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/oxpy
			ARCHIVE DESTINATION ${PYTHON_SITE_PACKAGES}/oxpy)
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/oxpy/
			DESTINATION ${PYTHON_SITE_PACKAGES}/oxpy
			FILES_MATCHING PATTERN "*.py")
	
	# add a target to test oxpy
	add_custom_target(test_oxpy
	    ${PROJECT_SOURCE_DIR}/test/TestSuite.py test_folder_list.txt ${PYTHON_EXECUTABLE} oxpy
	    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
	    COMMENT "Running oxpy tests" VERBATIM
	)
	# and one to test oat
	add_custom_target(test_oat
	    bash test.sh ${PYTHON_EXECUTABLE} && bash cleanup.sh
	    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/analysis/tests
	    COMMENT "Running oxpy tests" VERBATIM
	)
	add_dependencies(test test_oxpy test_oat)
	
ENDIF(Python)
