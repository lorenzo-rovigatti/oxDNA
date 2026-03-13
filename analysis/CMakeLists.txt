IF(Python)
    # First, check that the user's default Python meets the minimum version requirement.
	# Without this check, find_package may silently find a different Python installation
	# that meets the version requirement, causing oxpy to be installed for a Python
	# interpreter other than the one the user expects.
	find_program(DEFAULT_PYTHON NAMES python3 python)
	if(NOT DEFAULT_PYTHON)
		message(FATAL_ERROR "No python3 or python executable found in PATH.")
	endif()

	# Pin find_package to the user's default Python, preventing it from picking a different one.
	set(Python_EXECUTABLE "${DEFAULT_PYTHON}" CACHE FILEPATH "Path to Python interpreter" FORCE)

	if(OxpyVersionOverride)
		message(STATUS "OxpyVersionOverride is ON: skipping oxDNA_analysis_tools (requires Python >= 3.11)")
	else()
		# Install oxDNA_analysis_tools using the user's default Python (already validated >= 3.11 by src/oxpy/)
		if(OxpySystemInstall)
			install(CODE "message(STATUS \"Installing oxDNA_analysis_tools system-wide via pip...\")
						execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR})")
		else()
			install(CODE "message(STATUS \"Installing oxDNA_analysis_tools for current user via pip...\")
						execute_process(COMMAND ${Python_EXECUTABLE} -m pip install --user ${CMAKE_CURRENT_SOURCE_DIR})")
		endif()
	endif()
ENDIF(Python)
