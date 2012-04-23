"""
Application-class that implements pyFoamChangeMixingPlaneBoundary.py

Change various mixingPlane interface attributes in
constant/polymesh/boundary file.

Author:
  Martin Beaudoin, Hydro-Quebec, 2012.  All rights reserved

"""

from PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from os import path
import sys

class ChangeMixingPlaneBoundary(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Change MixingPlane boundary condition parameters
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory> mixingPlanePatchName",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=2)

    def addOptions(self):
        self.parser.add_option("--shadowPatch",
                               action="store",
                               dest="shadowPatch",
                               default=None,
                               help='Name of the shadowPatch')
        self.parser.add_option("--coordinateSystemName",
                               action="store",
                               dest="coordinateSystemName",
                               default="mixingCS",
                               help='coordinateSystemName (mixingCS)')
        self.parser.add_option("--coordinateSystemType",
                               action="store",
                               dest="coordinateSystemType",
                               default=None,
                               help='coordinateSystemType (cyindrical/spherical)')
        self.parser.add_option("--coordinateSystemOrigin",
                               action="store",
                               dest="coordinateSystemOrigin",
                               default=None,
                               help='origin for coordinate system of mixingPlane')
        self.parser.add_option("--coordinateSystemE1",
                               action="store",
                               dest="coordinateSystemE1",
                               default=None,
                               help='axis E1 for coordinate system of mixingPlane')
        self.parser.add_option("--coordinateSystemE3",
                               action="store",
                               dest="coordinateSystemE3",
                               default=None,
                               help='axis E3 for coordinate system of mixingPlane')
        self.parser.add_option("--assembly",
                               action="store",
                               dest="assembly",
                               default=None,
                               help='Assembly (master|slave|both|userdefined')
        self.parser.add_option("--orientation",
                               action="store",
                               dest="orientation",
                               default=None,
                               help='Orientation of profile (\
dirX_spanY| \
dirX_spanZ|dirY_spanX|dirY_spanZ|dirZ_spanX|dirZ_spanY|dirR_spanTheta|dirR_spanZ|dirTheta_spanZ|dirTheta_spanR|dirZ_spanTheta|dirZ_spanR|unknown)')
        
        self.parser.add_option("--test",
                               action="store_true",
                               default=False,
                               dest="test",
                               help="Only print the new boundary file")

    def run(self):
        fName=self.parser.getArgs()[0]
        bName=self.parser.getArgs()[1]
 
        boundary=ParsedParameterFile(path.join(".",fName,"constant","polyMesh","boundary"),debug=False,boundaryDict=True)

        bnd=boundary.content

        if type(bnd)!=list:
            print "Problem with boundary file (not a list)"
            sys.exit(-1)

        found=False

        for val in bnd:
            if val==bName:
                found=True
            elif found:
                if val["type"]=="mixingPlane":
                    if self.parser.getOptions().shadowPatch!=None:
                        val["shadowPatch"]=self.parser.getOptions().shadowPatch

                    if self.parser.getOptions().orientation!=None:
                        val["orientation"]=self.parser.getOptions().orientation

                    if val.has_key("coordinateSystem")==False:
                        val["coordinateSystem"]={}

                    if self.parser.getOptions().coordinateSystemName!=None:
                        val["coordinateSystem"]["name"]=self.parser.getOptions().coordinateSystemName

                    if self.parser.getOptions().coordinateSystemType!=None:
                        val["coordinateSystem"]["type"]=self.parser.getOptions().coordinateSystemType

                    if self.parser.getOptions().coordinateSystemOrigin!=None:
                        val["coordinateSystem"]["origin"]=self.parser.getOptions().coordinateSystemOrigin

                    if self.parser.getOptions().coordinateSystemE1!=None:
                        val["coordinateSystem"]["e1"]=self.parser.getOptions().coordinateSystemE1

                    if self.parser.getOptions().coordinateSystemE3!=None:
                        val["coordinateSystem"]["e3"]=self.parser.getOptions().coordinateSystemE3

                    if self.parser.getOptions().assembly!=None:
                        val["assembly"]=self.parser.getOptions().assembly

                    if self.parser.getOptions().orientation!=None:
                        val["orientation"]=self.parser.getOptions().orientation

                break

        if not found:
            print "Boundary",bName,"not found in",bnd[::2]
            sys.exit(-1)

        if self.parser.getOptions().test:
            print boundary
        else:
            boundary.writeFile()

