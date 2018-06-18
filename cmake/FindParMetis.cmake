
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(PARMETIS_PKGCONF parmetis)

# Include dir
find_path(PARMETIS_INCLUDE_DIR
  NAMES parmetis.h
  HINTS
  ThirdParty/packages/parmetis-4.0.3/platforms/linux64GccDPOpt/include
  ${PARMETIS_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(PARMETIS_LIBRARY
  NAMES parmetis
  HINTS
  ThirdParty/packages/parmetis-4.0.3/platforms/linux64GccDPOpt/lib
  ${PARMETIS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

set(PARMETIS_PROCESS_INCLUDES PARMETIS_INCLUDE_DIR)
set(PARMETIS_PROCESS_LIBS PARMETIS_LIBRARY)
libfind_process(PARMETIS)
