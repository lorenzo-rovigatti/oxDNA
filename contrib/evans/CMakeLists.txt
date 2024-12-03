# cmakelists.txt for Josh Evans oxDNA contributions.
# file created: 3 Dec. 2024, based on rovigatti/CMakeLists.txt
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(evans DEPENDS RaspberryPatchyInteraction)

# copying this next bit directly from Lorenzo's CMakeLists.txt
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

# interactions
add_library_no_prefix(RaspberryPatchyInteraction src/Interactions/RaspberryPatchyInteraction.cpp)