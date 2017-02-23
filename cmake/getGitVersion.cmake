
if(NOT GIT_FOUND)
  find_package(Git QUIET)
endif()

if(GIT_FOUND)
  # Try to get version from from git
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty=-dirty
	OUTPUT_VARIABLE GIT_VERSION
	ERROR_VARIABLE dummy
	RESULT_VARIABLE res
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  if(res EQUAL 0)
    string(REPLACE "g" "" GIT_VERSION "${GIT_VERSION}")
    string(REGEX REPLACE "^v([0-9]+\\.?[0-9]*\\.?([0-9]*)).*" "\\1" FOAM_VERSION "${GIT_VERSION}")
  endif()
endif()

if(NOT GIT_VERSION)
  if(EXISTS cmake/storedVersion.cmake)
    # Fall-back to stored version
    include(cmake/storedVersion.cmake)
  else()
    # We should not be here. Set some defaults!
    set(GIT_VERSION "4.0-unknown")
    set(FOAM_VERSION "4.0")
  endif()
endif()

# Save version information so that it can be exported with the project
# and included as a fall-back - see above
configure_file(cmake/storedVersion.cmake.in
  cmake/storedVersion.cmake
)

# Configure global.C
configure_file(src/foam/global/global.C.in
  src/foam/global/global.C
)

