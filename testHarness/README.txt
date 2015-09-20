Description
===========
This part of the repository is for FOAM test harnesses.


Directory Structure
===================

foam-extend-3.2    : Test harness for foam-extend version 3.2. See the file foam-extend-3.2/README.txt for more information

OSIG                : CMake/CTest scripts for FOAM Special Interest Group (OSIG) test harness
OSIG/TurboMachinery : Test harness for the TurboMachinery OSIG. See the file OSIG/Turbomachinery/README.txt for more information.


Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved.


1: Select your git branch of choice: ie:

   git checkout master   # For Hrv master branch


2: Make sure your FOAM environment is properly configured to run FOAM.


3: The minimal cmake version number for running the test loop is 2.8.0. Make sure you are there.


4: Make sure you have the test harness scripts available under your git branch;
   otherwise, you will need to fetch this from Hrv's master branch, and merge it into yours

   ls $WM_PROJECT_DIR/testHarness  # Checking availability of testHarness under this branch


5: move to the runDir directory for the FOAM test harness

   cd $WM_PROJECT_DIR/testHarness/foam-extend/3.2/runDir


6: Normally, if using the master branch, everything should already be setup for you to run the test harness.
   Still, I recommand always checking that these two important files are up-to-date:

   cp ../CMakeFiles/CMakeLists.txt $WM_PROJECT_DIR
   cp ../CMakeFiles/CTestConfig.cmake.foam-extend $WM_PROJECT_DIR/CTestConfig.cmake


7:   Next, running the test loop is pretty simple:

   cd $WM_PROJECT_DIR/testHarness/foam-extend/3.2/runDir # you should already be there...
   ./Allclean
   ./Allrun_Experimental


8: The results will be published on the CDash dashboard on foam-extend.

   To see your results:
   URL      : http://foam-extend.sourceforge.net/CDash/index.php?project=foam-extend-3.2


9: You can customize your system identifier on the dashboard using the environment variable $CDASH_SUBMIT_LOCAL_HOST_ID.
   Otherwise, the fully qualified name of your system will be used.

   A good customization idea would be to add the name of your git branch in your system ID.
   I will probably modify my scripts to add this information automagically.

   NB: Please no "forward slash" or "/" in the system ID; it looks like CDash will choke on this.


10: In general, see the file $WM_PROJECT_DIR/testHarness/foam-extend/3.2/README.txt for the necessary information about running the
    test loop.


11: Please do not hesitate to report any problems, comments, suggestions about the test loop.
