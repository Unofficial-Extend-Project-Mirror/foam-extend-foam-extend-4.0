
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was FOAMConfig.cmake.in                            ########

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

include("${CMAKE_CURRENT_LIST_DIR}/FOAMConfigVersion.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/FOAMMacros.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/FOAMTargets.cmake")

message(STATUS "Importing FOAM ${FOAM_VERSION} from ${CMAKE_CURRENT_LIST_DIR}")

include(CMakeFindDependencyMacro)

find_dependency(MPI)
add_library(mpi SHARED IMPORTED)
set_property(TARGET mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MPI_C_INCLUDE_PATH})
set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${MPI_LIBRARY})

find_dependency(ZLIB)

#set_and_check(FOAM_INCLUDE_DIR "")
#set_and_check(FOAM_LIBRARY foam)

# Once we split into components
#set(_supported_components Plot Table)

#foreach(_comp ${FOAM_FIND_COMPONENTS})
#  if (NOT ";${_supported_components};" MATCHES _comp)
#    set(FOAM_FOUND False)
#    set(FOAM_NOTFOUND_MESSAGE "Unsupported component: ${_comp}")
#  endif()
#  include("${CMAKE_CURRENT_LIST_DIR}/FOAM${_comp}Targets.cmake")
#endforeach()

check_required_components(FOAM)
