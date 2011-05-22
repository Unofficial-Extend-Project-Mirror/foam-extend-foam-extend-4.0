Description
===========
This part of the repository is for OpenFOAM test harnesses.


Directory Structure
===================

OpenFOAM            : CMake/CTest scripts for compilation and execution test harness for OpenFOAM.
OpenFOAM/1.6-ext    : Test harness for OpenFOAM version 1.6-ext. See the file OpenFOAM/1.6-ext/README.txt for more information

OSIG                : CMake/CTest scripts for OpenFOAM Special Interest Group (OSIG) test harness
OSIG/TurboMachinery : Test harness for the TurboMachinery OSIG. See the file OSIG/Turbomachinery/README.txt for more information.


Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved.

 
1: Select your git branch of choice: ie:
   
   git checkout master   # For Hrv master branch


2: Make sure your OpenFOAM environment is properly configured to run OpenFOAM.


3: The minimal cmake version number for running the test loop is 2.8.0. Make sure you are there. 


4: Make sure you have the test harness scripts available under your git branch; 
   otherwise, you will need to fetch this from Hrv's master branch, and merge it into yours

   ls $WM_PROJECT_DIR/testHarness  # Checking availability of testHarness under this branch


5: move to the runDir directory for the OpenFOAM test harness
   
   cd $WM_PROJECT_DIR/testHarness/OpenFOAM/1.6-ext/runDir


6: Normally, if using the master branch, everything should already be setup for you to run the test harness.
   Still, I recommand always checking that these two important files are up-to-date:

   cp ../CMakeFiles/CMakeLists.txt $WM_PROJECT_DIR
   cp ../CMakeFiles/CTestConfig.cmake.openfoam-extend $WM_PROJECT_DIR/CTestConfig.cmake 


7:   Next, running the test loop is pretty simple:
 
   cd $WM_PROJECT_DIR/testHarness/OpenFOAM/1.6-ext/runDir # you should already be there...
   ./Allclean
   ./Allrun_Experimental


8: The results will be published on the CDash dashboard on openfoam-extend. 

   To see your results:
   URL      : http://openfoam-extend.sourceforge.net/CDash/index.php?project=OpenFOAM-1.6-ext


9: You can customize your system identifier on the dashboard using the environment variable $CDASH_SUBMIT_LOCAL_HOST_ID. 
   Otherwise, the fully qualified name of your system will be used. 

   A good customization idea would be to add the name of your git branch in your system ID. 
   I will probably modify my scripts to add this information automagically.

   NB: Please no "forward slash" or "/" in the system ID; it looks like CDash will choke on this.


10: In general, see the file $WM_PROJECT_DIR/testHarness/OpenFOAM/1.6-ext/README.txt for the necessary information about running the 
    test loop. 


11: Please do not hesitate to report any problems, comments, suggestions about the test loop. 
