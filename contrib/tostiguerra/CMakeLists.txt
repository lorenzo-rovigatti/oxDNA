SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(tostiguerra
	DEPENDS CGNucleicAcidsInteraction
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
add_library_no_prefix(CGNucleicAcidsInteraction MODULE EXCLUDE_FROM_ALL src/CGNucleicAcidsInteraction.cpp)
