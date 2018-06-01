
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CCMIO_PKGCONF mesquite)

# Include dir
find_path(CCMIO_INCLUDE_DIR
  NAMES libccmio/ccmio.h
  HINTS
  ThirdParty/packages/libccmio-2.6.1/platforms/linux64GccDPOpt/include
  ${MESQUITE_PKGCONF_INCLUDE_DIRS}
)

find_path(CCMIO_INCLUDE_DIR2
  NAMES ccmio.h
  HINTS
  ThirdParty/packages/libccmio-2.6.1/platforms/linux64GccDPOpt/include/libccmio
  ${MESQUITE_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(CCMIO_LIBRARY
  NAMES ccmio
  HINTS
  ThirdParty/packages/libccmio-2.6.1/platforms/linux64GccDPOpt/lib
  ${CCMIO_PKGCONF_LIBRARY_DIRS}
)

find_library(CCMIO_LIBRARY2
  NAMES adf_ccmio
  HINTS
  ThirdParty/packages/libccmio-2.6.1/platforms/linux64GccDPOpt/lib
  ${CCMIO_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

set(CCMIO_PROCESS_INCLUDES CCMIO_INCLUDE_DIR CCMIO_INCLUDE_DIR2)
set(CCMIO_PROCESS_LIBS CCMIO_LIBRARY CCMIO_LIBRARY2)
libfind_process(CCMIO)
