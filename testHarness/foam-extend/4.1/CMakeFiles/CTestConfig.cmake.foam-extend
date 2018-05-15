## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

# These settings will allow you to publish your dashboards results
# on the foam-extend CDash service hosted on SourceForge.Net
# See here: http://foam-extend.sourceforge.net/CDash/index.php
#
set(CTEST_PROJECT_NAME "foam-extend-4.1")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "foam-extend.sourceforge.net")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=foam-extend-4.1")
set(CTEST_DROP_SITE_CDASH TRUE)

# We can override those variables for local sites so you can use
# your own site CDash service
# This optional file will be located here:
# $FOAM_SITE_DIR/etc/CTestConfig.site.cmake
#
include($ENV{FOAM_SITE_DIR}/etc/CTestConfig.site.cmake OPTIONAL)

# We can override those variables from user space so you can use
# your own personal CDash service
# This optional file will be located here:
# $WM_PROJECT_USER_DIR/etc/CTestConfig.user.cmake
#
include($ENV{WM_PROJECT_USER_DIR}/etc/CTestConfig.user.cmake OPTIONAL)
