"""
Application-class that implements pyFoamChangeMixingPlaneBoundary.py

Change various mixingPlane interface attributes in
constant/polymesh/boundary file.

Author:
  Martin Beaudoin, Hydro-Quebec, 2012.  All rights reserved

"""

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.ThirdParty.six import print_
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
        self.parser.add_option("--zone",
                               action="store",
                               dest="zone",
                               default=None,
                               help='Name of the zone for mixingPlanePatch')
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
        self.parser.add_option("--ribbonPatchSweepAxis",
                               action="store",
                               dest="ribbonPatchSweepAxis",
                               default=None,
                               help='ribbonPatch sweepAxis (X|Y|Z|R|Theta')
        self.parser.add_option("--ribbonPatchStackAxis",
                               action="store",
                               dest="ribbonPatchStackAxis",
                               default=None,
                               help='ribbonPatch stackAxis (X|Y|Z|R|Theta')
        self.parser.add_option("--ribbonPatchDiscretisation",
                               action="store",
                               dest="ribbonPatchDiscretisation",
                               default=None,
                               help='ribbonPatch discretisation (masterPatch|slavePatch|bothPatches|uniform|userDefined)')

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
            print_("Problem with boundary file (not a list)")
            sys.exit(-1)

        found=False

        for val in bnd:
            if val==bName:
                found=True
            elif found:
                if val["type"]=="mixingPlane":
                    if self.parser.getOptions().shadowPatch!=None:
                        val["shadowPatch"]=self.parser.getOptions().shadowPatch

                    if self.parser.getOptions().zone!=None:
                        val["zone"]=self.parser.getOptions().zone

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

                    if val.has_key("ribbonPatch")==False:
                        val["ribbonPatch"]={}

                    if self.parser.getOptions().ribbonPatchSweepAxis!=None:
                        val["ribbonPatch"]["sweepAxis"]=self.parser.getOptions().ribbonPatchSweepAxis

                    if self.parser.getOptions().ribbonPatchStackAxis!=None:
                        val["ribbonPatch"]["stackAxis"]=self.parser.getOptions().ribbonPatchStackAxis

                    if self.parser.getOptions().ribbonPatchDiscretisation!=None:
                        val["ribbonPatch"]["discretisation"]=self.parser.getOptions().ribbonPatchDiscretisation

                break

        if not found:
            print_("Boundary",bName,"not found in",bnd[::2])
            sys.exit(-1)

        if self.parser.getOptions().test:
            print_(boundary)
        else:
            boundary.writeFile()
