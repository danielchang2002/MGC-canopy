#Project and build process settings
cmake_minimum_required(VERSION 3.6)
project(mgs-canopy-clustering)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Release) #Otherwise it's debug, Release uses -O3. If you want a one time Debug build use -D CMAKE_BUILD_TYPE=Debug flag when running cmake

OPTION(SINGLEPRECISION "Use floats instead of doubles. Saves memory but reduces precision" ON)
IF(SINGLEPRECISION)
    message("Using reduced precision type (float instead of double)")
    ADD_DEFINITIONS(-DSINGLEPRECISION)
ENDIF(SINGLEPRECISION)

#General compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS_RELEASE "-O3 -msse4.2") #O3 is just here to make sure, msse4.2 might feel a bit stringent but it should be ok in practice

#Source files 
set(SOURCE_FILES src/CanopyClustering.cpp src/Canopy.cpp src/Log.cpp src/Point.cpp src/signal_handlers.cpp src/Stats.cpp src/TimeProfile.cpp)
set(HEADER_FILES src/CanopyClustering.hpp src/Canopy.hpp src/Log.hpp src/Point.hpp src/signal_handlers.hpp src/Stats.hpp src/TimeProfile.hpp src/prog_bar_misc.hpp src/program_options_misc.hpp src/precision_type.hpp)

#BOOST
FIND_PACKAGE( Boost 1.59 COMPONENTS program_options iostreams unit_test_framework REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

#Zlib
FIND_PACKAGE( ZLIB REQUIRED )
INCLUDE_DIRECTORIES( ${ZLIB_INCLUDE_DIRS} )

#OpenMP
FIND_PACKAGE( OpenMP REQUIRED)

#Main targets
add_library(libcanopy ${HEADER_FILES} ${SOURCE_FILES})
add_executable(canopy src/main.cpp)
TARGET_LINK_LIBRARIES(canopy ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} OpenMP::OpenMP_CXX libcanopy)

#Test targets
add_definitions(-DBOOST_TEST_DYN_LINK)
add_executable(canopy_unit_test test/unit_test.cpp)
TARGET_LINK_LIBRARIES(canopy_unit_test ${Boost_LIBRARIES} OpenMP::OpenMP_CXX libcanopy)


#Copy the canopy executable after build to the root directoryadd_executable(MyExe ...)
add_custom_command(TARGET canopy POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:canopy> ../)
add_custom_command(TARGET canopy_unit_test POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:canopy_unit_test> ../)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ../canopy)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ../canopy_unit_test)
