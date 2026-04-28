
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was HepMC3Config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

SET(HEPMC3_VERSION 3.03.00)
SET(HEPMC3_VERSION_MAJOR  3)
SET(HEPMC3_VERSION_MINOR  3)
SET(HEPMC3_VERSION_PATCH  0)


set_and_check(HEPMC3_INCLUDE_DIR ${PACKAGE_PREFIX_DIR}/include)
set(HEPMC3_CXX_STANDARD 20)
set(HEPMC3_FEATURES search rootIO examples python)
set(HEPMC3_COMPONENTS search;rootIO)

if(EXISTS ${PACKAGE_PREFIX_DIR}/share/HepMC3/interfaces)
  set(HEPMC3_INTERFACES_DIR ${PACKAGE_PREFIX_DIR}/share/HepMC3/interfaces)
endif()

find_library(HEPMC3_LIB NAMES HepMC3 HINTS ${PACKAGE_PREFIX_DIR}/lib64)
find_library(HEPMC3_SEARCH_LIB NAMES HepMC3search HINTS ${PACKAGE_PREFIX_DIR}/lib64)
find_library(HEPMC3_ROOTIO_LIB NAMES HepMC3rootIO HINTS ${PACKAGE_PREFIX_DIR}/lib64)
find_library(HEPMC3_PROTOBUFIO_LIB NAMES HepMC3protobufIO HINTS ${PACKAGE_PREFIX_DIR}/lib64)

set(HEPMC3_LIBRARIES ${HEPMC3_LIB})
if(EXISTS ${HEPMC3_SEARCH_LIB})
  list( APPEND  HEPMC3_LIBRARIES ${HEPMC3_SEARCH_LIB})
endif()
if(EXISTS ${HEPMC3_ROOTIO_LIB})
  list( APPEND  HEPMC3_LIBRARIES ${HEPMC3_ROOTIO_LIB})
endif()
if(EXISTS ${HEPMC3_PROTOBUFIO_LIB})
  list( APPEND  HEPMC3_LIBRARIES ${HEPMC3_PROTOBUFIO_LIB})
endif()

include(${CMAKE_CURRENT_LIST_DIR}/HepMC3Targets.cmake)

if (TARGET HepMC3All)
  message(STATUS "WARNING: Please note that multiple calls to find_package(HepMC3) are not recommended!")
else()
add_library(HepMC3All INTERFACE)
target_link_libraries(HepMC3All INTERFACE HepMC3::HepMC3)

if(TARGET HepMC3::HepMC3_static)
  add_library(HepMC3All_static INTERFACE)
  target_link_libraries(HepMC3All_static INTERFACE HepMC3::HepMC3_static)
endif()

foreach(_comp ${HEPMC3_COMPONENTS})
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/HepMC3${_comp}Targets.cmake)
    include(${CMAKE_CURRENT_LIST_DIR}/HepMC3${_comp}Targets.cmake)
    target_link_libraries(HepMC3All INTERFACE HepMC3::${_comp})
    if(TARGET HepMC3::${_comp}_static)
      target_link_libraries(HepMC3All_static INTERFACE HepMC3::${_comp}_static)
    endif()
  endif()
endforeach()

add_library(HepMC3::All ALIAS HepMC3All)
if(TARGET HepMC3All_static)
  add_library(HepMC3::All_static ALIAS HepMC3All_static)
endif()

foreach(_comp ${HepMC3_FIND_COMPONENTS})
  set(_comps "search;rootIO")
  if(NOT _comp IN_LIST _comps)
    message(WARNING "Unsupported component ${_comp}")
    set(HepMC3_${_comp}_FOUND False)
  else()
    if(TARGET HepMC3::${_comp})
      set(HepMC3_${_comp}_FOUND True)
    else()
      message(WARNING "Uninstalled component ${_comp}")
      set(HepMC3_${_comp}_FOUND False)
    endif()
  endif()
endforeach()
endif()

check_required_components(HepMC3)
