
include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(PARMGRIDGEN_PKGCONF scotch)

# Include dir
find_path(PARMGRIDGEN_INCLUDE_DIR
  NAMES mgridgen.h
  HINTS
  ThirdParty/packages/ParMGridGen-1.0/platforms/linux64GccDPOpt/include/Lib
  ${PARMGRIDGEN_PKGCONF_INCLUDE_DIRS}
)

find_path(PARMGRIDGEN_INCLUDE_DIR2
  NAMES IMlib.h
  HINTS
  ThirdParty/packages/ParMGridGen-1.0/platforms/linux64GccDPOpt/include/IMlib
  ${PARMGRIDGEN_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(PARMGRIDGEN_LIBRARY
  NAMES MGridGen
  HINTS
  ThirdParty/packages/ParMGridGen-1.0/platforms/linux64GccDPOpt/lib
  ${PARMGRIDGEN_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.

set(PARMGRIDGEN_PROCESS_INCLUDE PARMGRIDGEN_INCLUDE_DIR)
set(PARMGRIDGEN_PROCESS_INCLUDES PARMGRIDGEN_INCLUDE_DIR2)
set(PARMGRIDGEN_PROCESS_LIB PARMGRIDGEN_LIBRARY)
libfind_process(PARMGRIDGEN)
