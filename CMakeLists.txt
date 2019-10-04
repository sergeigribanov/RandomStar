cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(RndStar)
option(MEUTILS_DIR "Matrix Element Utils instalation path")
if (NOT MEUTILS_DIR)
  if (DEFINED ENV{MEUTILS_DIR})
    set(MEUTILS_DIR $ENV{MEUTILS_DIR})
  else()
    message(FATAL_ERROR "[!] MEUTILS_DIR option is not set!
You need to set this option or to define MEUTILS_DIR environment variable.")
  endif()
endif()
find_package(GSL REQUIRED)
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
link_directories(${MEUTILS_DIR}/lib)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB SOURCES_NOMAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SOURCES_NOMAIN main.cpp)
set(PUBLIC_HEADERS ${HEADERS})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rndstar-setup.sh.in rndstar-setup.sh)
add_executable(rstar-test ${SOURCES})
target_include_directories(rstar-test
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${MEUTILS_DIR}/inc
  ${GSL_INCLUDE_DIRS})
target_link_libraries(rstar-test
  -lMEUtils
  GSL::gsl
  GSL::gslcblas
  ${ROOT_LIBRARIES})
add_library(rstar SHARED ${SOURCES_NOMAIN})
target_include_directories(rstar
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
  ${MEUTILS_DIR}/inc
  ${GSL_INCLUDE_DIRS})
target_link_libraries(rstar
  -lMEUtils
  GSL::gsl
  GSL::gslcblas
  ${ROOT_LIBRARIES})
set_target_properties(rstar
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS rstar-test DESTINATION bin)
install(TARGETS rstar
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/rndstar-setup.sh DESTINATION env)
