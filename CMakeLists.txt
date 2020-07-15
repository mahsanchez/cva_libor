cmake_minimum_required(VERSION 3.15)
project(untitled)

set(CMAKE_FIND_DEBUG_MODE 10)

set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES main.cpp)
set(SOURCE_FILES calibration.h)
set(SOURCE_FILES calibration_libor.h)
set(SOURCE_FILES simulation.h)
set(SOURCE_FILES product.h)
set(SOURCE_FILES q_numerics.h)
set(SOURCE_FILES curvepoint.h)
set(SOURCE_FILES leastsquares.h)

file(GLOB SOURCES *.h *.cpp *.hpp)

add_executable(untitled ${SOURCE_FILES} ${SOURCES})

#Boost Libraries
#set(BOOST_ROOT /opt/boost/boost_1_72_0)
#set(BOOST_INCLUDE_DIR /opt/boost/boost_1_72_0/include)
#set(BOOST_LIBRARY_DIR /opt/boost/boost_1_72_0/lib)

set(BOOST_ROOT /usr/local/Cellar/boost/1.72.0_3)
set(BOOST_INCLUDE_DIR /usr/local/Cellar/boost/1.72.0_3/include)
set(BOOST_LIBRARY_DIR /usr/local/Cellar/boost/1.72.0_3/lib)

include_directories(${BOOST_ROOT})
include_directories(${BOOST_INCLUDE_DIR})
link_directories(${BOOST_LIBRARY_DIR})

#QuantLib Libraries
set(QUANTLIB_ROOT /usr/local/Cellar/quantlib/1.18)
set(QUANTLIB_INCLUDE_DIR /usr/local/Cellar/quantlib/1.18/include)
set(QUANTLIB_LIBRARY_DIR /usr/local/Cellar/quantlib/1.18/lib)

include_directories(${QUANTLIB_ROOT})
include_directories(${QUANTLIB_INCLUDE_DIR})
link_directories(${QUANTLIB_LIBRARY_DIR})

# Find MKL ROOT
set(MKLROOT /opt/intel/mkl)
set(MKL_ROOT_DIR /opt/intel/mkl)
set(MKL_LIBRARIES /opt/intel/mkl/lib)
set(MKL_LIBRARY_DIR /opt/intel/mkl/lib)

set(MKL_INCLUDE_DIR ${MKL_ROOT_DIR}/include)
set(MKL_LINK_TOOL ${MKL_ROOT_DIR}/tools/mkl_link_tool)

set(MKL_LINK_TOOL_COMMAND ${MKL_LINK_TOOL} "-libs")
list(APPEND MKL_LINK_TOOL_COMMAND "--compiler=clang")
list(APPEND MKL_LINK_TOOL_COMMAND "--os=mac")
list(APPEND MKL_LINK_TOOL_COMMAND "--arch=intel64")
list(APPEND MKL_LINK_TOOL_COMMAND "--parallel=no")

execute_process(COMMAND ${MKL_LINK_TOOL_COMMAND}
        OUTPUT_VARIABLE MKL_LIBS
        RESULT_VARIABLE COMMAND_WORKED
        TIMEOUT 2 ERROR_QUIET)

set(MKL_LIBRARIES)
set(MKL_LIBRARY_DIR)

if (NOT ${COMMAND_WORKED} EQUAL 0)
    message(FATAL_ERROR "Cannot find the MKL libraries correctly. Please check your MKL input variables and mkl_link_tool. The command executed was:\n ${MKL_LINK_TOOL_COMMAND}.")
endif()

string(REGEX REPLACE "\n" "" MKL_LIBS ${MKL_LIBS})
string(REPLACE "$(MKLROOT)" "${MKL_ROOT_DIR}" MKL_LIBRARIES ${MKL_LIBS})
# hack for lin with libiomp5.a
if (APPLE)
    string(REPLACE "-liomp5" "${MKL_ROOT_DIR}/../compiler/lib/libiomp5.a" MKL_LIBRARIES ${MKL_LIBRARIES})
endif()
separate_arguments(MKL_LIBRARIES)

string(REPLACE "-libs" "-opts" MKL_LINK_TOOL_COMMAND "${MKL_LINK_TOOL_COMMAND}")
execute_process(COMMAND ${MKL_LINK_TOOL_COMMAND} OUTPUT_VARIABLE RESULT_OPTS TIMEOUT 2 ERROR_QUIET)
string(REGEX MATCHALL "[-/]D[^\ ]*" MKL_DEFINITIONS ${RESULT_OPTS})

if (CMAKE_FIND_DEBUG_MODE)
    message(STATUS "Exectuted command: \n${MKL_LINK_TOOL_COMMAND}")
    message(STATUS "Found MKL_LIBRARIES:\n${MKL_LIBRARIES} ")
    message(STATUS "Found MKL_DEFINITIONS:\n${MKL_DEFINITIONS} ")
    message(STATUS "Found MKL_LIBRARY_DIR:\n${MKL_LIBRARY_DIR} ")
    message(STATUS "Found MKL_INCLUDE_DIR:\n${MKL_INCLUDE_DIR} ")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)
mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES MKL_DEFINITIONS MKL_ROOT_DIR)

include_directories(
        ${MKL_INCLUDE_DIR}
        ${QUANTLIB_INCLUDE_DIR}
)
link_directories(
        ${MKL_LIBRARIES}
        ${QUANTLIB_LIBRARY_DIR}
)

target_link_libraries(untitled  "-L${MKLROOT}/lib -Wl,-rpath,$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl" ${MKL_LIBRARIES} "-L${QUANTLIB_LIBRARY_DIR} -lQuantLib"  ${BOOST_LIBRARY_DIR})


