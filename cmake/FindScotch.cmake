
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(SCOTCH_PKGCONF scotch)

# Include dir
find_path(SCOTCH_INCLUDE_DIR
  NAMES scotch.h
  HINTS
  ThirdParty/packages/scotch-6.0.4/platforms/linux64GccDPOpt/include
  ${SCOTCH_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(SCOTCH_LIBRARY
  NAMES scotch
  HINTS
  ThirdParty/packages/scotch-6.0.4/platforms/linux64GccDPOpt/lib
  ${SCOTCH_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

message(STATUS ${SCOTCH_INCLUDE_DIR})
message(STATUS ${SCOTCH_LIBRARY})

set(SCOTCH_PROCESS_INCLUDES SCOTCH_INCLUDE_DIR)
set(SCOTCH_PROCESS_LIBS SCOTCH_LIBRARY)
libfind_process(SCOTCH)
