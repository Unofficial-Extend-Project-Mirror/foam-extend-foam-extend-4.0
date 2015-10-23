"""
Application-class that implements pyFoamConvertMixingPlaneToNewSyntax.py

Adjust the mixingPlane interface definition in the boundary
file to the latest supported syntax.

Author:
  Martin Beaudoin, Hydro-Quebec, 2012.  All rights reserved

"""

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.ThirdParty.six import print_
from os import path
import sys

# ------> Start of Python code snippet copied from an external source
#
# Author : Brian Beck
# http://code.activestate.com/recipes/410692-readable-switch-construction-without-lambdas-or-di/
#
# License for this Python code snippet:
#
# This code that was deposited before July 15, 2008 on aspn.activestate.com.
# It is governed by the Python license (http://www.python.org/psf/license/)
# in accordance with the Python Cookbook agreement.
#
# Description:
# Python's lack of a 'switch' statement has garnered much discussion and even
# a PEP. The most popular substitute uses dictionaries to map cases to
# functions, which requires lots of defs or lambdas. While the approach shown
# here may be O(n) for cases, it aims to duplicate C's original 'switch'
# functionality and structure with reasonable accuracy.
#
# This class provides the functionality we want. You only need to look at
# this if you want to know how this works. It only needs to be defined
# once, no need to muck around with its internals.
#

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

# ------> End of Python code snippet copied from an external source

####################################################################
#
# The rest of this source code was written by me.
# Martin Beaudoin, June 2012

class ConvertMixingPlaneBoundaryToNewSyntax(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Change MixingPlane boundary condition parameters
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory>",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=1)

    def addOptions(self):

        self.parser.add_option("--test",
                               action="store_true",
                               default=False,
                               dest="test",
                               help="Only print the new boundary file")

    def run(self):
        fName=self.parser.getArgs()[0]

        boundary=ParsedParameterFile(path.join(".",fName,"constant","polyMesh","boundary"),debug=False,boundaryDict=True)

        bnd=boundary.content

        if type(bnd)!=list:
            print_("Problem with boundary file (not a list)")
            sys.exit(-1)

        found=False

        for index in range(0, len(bnd), 2):

            indexDefPatch=index+1

            oldAssembly=""
            oldOrientation=""

            if bnd[indexDefPatch]["type"]=="mixingPlane":
                if bnd[indexDefPatch].has_key("assembly"):
                    print_("    Replacing the parameter 'assembly' for patch", bnd[index])
                    oldAssembly=bnd[indexDefPatch]["assembly"]
                    del bnd[indexDefPatch]["assembly"]

                if bnd[indexDefPatch].has_key("orientation"):
                    print_("    Replacing the parameter 'orientation' for patch", bnd[index])
                    oldOrientation=bnd[indexDefPatch]["orientation"]
                    del bnd[indexDefPatch]["orientation"]

                if bnd[indexDefPatch].has_key("ribbonPatch")==False:
                    bnd[indexDefPatch]["ribbonPatch"]={}

                if bnd[indexDefPatch].has_key("zone")==False:
                    bnd[indexDefPatch]["zone"]=bnd[index] + "Zone"

                if oldAssembly != "":
                    # Converting "assembly" to ribbonPatch/discretisation
                    for case in switch(oldAssembly):
                        if case('master'):
                            bnd[indexDefPatch]["ribbonPatch"]["discretisation"]="masterPatch"
                            break
                        if case('slave'):
                            bnd[indexDefPatch]["ribbonPatch"]["discretisation"]="slavePatch"
                            break
                        if case('both'):
                            bnd[indexDefPatch]["ribbonPatch"]["discretisation"]="bothPatches"
                            break
                        if case('userdefined'):
                            bnd[indexDefPatch]["ribbonPatch"]["discretisation"]="userDefined"
                            break
                        if case(): # default
                            print_("Unsupported assembly type: ", oldAssembly)

                if oldOrientation != "":
                    # Converting "orientation" to ribbonPatch/ribbonPatchSweepAxis and
                    # ribbonPatch/ribbonPatchStackAxis
                    for case in switch(oldOrientation):

                        if case('dirX_spanY'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="X"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Y"
                            break
                        if case('dirX_spanZ'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="X"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Z"
                            break
                        if case('dirY_spanX'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Y"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="X"
                            break
                        if case('dirY_spanZ'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Y"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Z"
                            break
                        if case('dirZ_spanX'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Z"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="X"
                            break
                        if case('dirZ_spanY'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Z"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Y"
                            break
                        if case('dirR_spanTheta'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="R"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Theta"
                            break
                        if case('dirR_spanZ'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="R"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Z"
                            break
                        if case('dirTheta_spanZ'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Theta"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Z"
                            break
                        if case('dirTheta_spanR'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Theta"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="R"
                            break
                        if case('dirZ_spanTheta'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Z"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="Theta"
                            break
                        if case('dirZ_spanR'):
                            bnd[indexDefPatch]["ribbonPatch"]["stackAxis"]="Z"
                            bnd[indexDefPatch]["ribbonPatch"]["sweepAxis"]="R"
                            break
                        if case(): # default
                            print_("Unsupported orientation type: ", oldOrientation)

        if self.parser.getOptions().test:
            print_(boundary)
        else:
            boundary.writeFile()
