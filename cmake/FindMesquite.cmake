
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(MESQUITE_PKGCONF mesquite)

# Include dir
find_path(MESQUITE_INCLUDE_DIR
  NAMES Mesquite_all_headers.hpp
  HINTS
  ThirdParty/packages/mesquite-2.1.2/platforms/linux64GccDPOpt/include
  ${MESQUITE_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(MESQUITE_LIBRARY
  NAMES mesquite
  HINTS
  ThirdParty/packages/mesquite-2.1.2/platforms/linux64GccDPOpt/lib
  ${MESQUITE_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

set(MESQUITE_PROCESS_INCLUDE MESQUITE_INCLUDE_DIR)
set(MESQUITE_PROCESS_LIB MESQUITE_LIBRARY)
libfind_process(MESQUITE)
