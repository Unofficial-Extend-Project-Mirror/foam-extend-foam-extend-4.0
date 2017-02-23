add_custom_target(ImportMessage
  ${CMAKE_COMMAND} -E cmake_echo_color --red --bold
  "Using FOAM ${FOAM_VERSION} in ${CMAKE_CURRENT_LIST_DIR}"
)

function(add_foam_library lib)
  # Create target for lnInclude and use it
  include_directories(lnInclude)
  add_custom_command(
    OUTPUT lnInclude/uptodate
    DEPENDS CMakeLists.txt
    COMMAND rm -rf lnInclude
    COMMAND wmakeLnInclude .
    COMMAND touch lnInclude/uptodate
  )
  add_custom_target(${lib}_lnInclude DEPENDS lnInclude/uptodate)

  # Add the library to the targets and set include paths
  add_library(${lib} ${ARGN})
  add_dependencies(${lib} ${lib}_lnInclude)
  target_include_directories(${lib} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lnInclude>
    $<INSTALL_INTERFACE:include/${lib}>
  )

  # Install library and headers to install location
  install(TARGETS ${lib} DESTINATION lib EXPORT FOAMTargets)
  file(GLOB _files FOLLOW_SYMLINKS "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude/*")
  set (_resolvedFiles "")
  foreach (_file ${_files})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
  endforeach()
  install(FILES ${_resolvedFiles} DESTINATION include/${lib})

  # Export target to include them from externals builds to build location
  export(TARGETS ${lib} APPEND FILE ${CMAKE_BINARY_DIR}/cmake/FOAMTargets.cmake)

  #generate_export_header(${lib})
  set_property(TARGET ${lib} PROPERTY VERSION ${FOAM_VERSION})
  set_property(TARGET ${lib} PROPERTY SOVERSION 3)
  set_property(TARGET ${lib} PROPERTY INTERFACE_FOAM_MAJOR_VERSION 4)
  set_property(TARGET ${lib} APPEND PROPERTY
    COMPATIBLE_INTERFACE_STRING FOAM_MAJOR_VERSION
  )

  if(NOT DEFINED HAVE_FOAM)
    add_dependencies(${lib} ImportMessage)
  endif()
endfunction()

function(add_foam_executable exe)
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs DEPENDS SOURCES)
  cmake_parse_arguments(AFE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  add_executable(${exe} ${AFE_SOURCES})
  target_link_libraries(${exe} ${AFE_DEPENDS})
  install (TARGETS ${exe} DESTINATION bin)

  if(NOT DEFINED HAVE_FOAM)
    add_dependencies(${exe} ImportMessage)
  endif()
endfunction()


