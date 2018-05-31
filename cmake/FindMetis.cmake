
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(METIS_PKGCONF metis)

# Include dir
find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  HINTS
  ThirdParty/packages/metis-5.1.0/platforms/linux64GccDPOpt/include
  ${METIS_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(METIS_LIBRARY
  NAMES metis
  HINTS
  ThirdParty/packages/metis-5.1.0/platforms/linux64GccDPOpt/lib
  ${METIS_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

set(METIS_PROCESS_INCLUDES METIS_INCLUDE_DIR)
set(METIS_PROCESS_LIBS METIS_LIBRARY)
libfind_process(METIS)
