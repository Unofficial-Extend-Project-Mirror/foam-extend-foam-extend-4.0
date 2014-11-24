"""
Application-class that implements pyFoamChangeGGIBoundary.py

Modification of GGI and cyclicGGI interface parameters in
constant/polymesh/boundary file.

Author:
  Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

"""

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.ThirdParty.six import print_
from os import path
import sys
import re

class ChangeGGIBoundary(PyFoamApplication):
    def __init__(self,args=None):
        description="""\
Change GGI boundary condition parameters
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory> ggiPatchName",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=2)

    def addOptions(self):
        self.parser.add_option("--shadowPatch",
                               action="store",
                               dest="shadowPatch",
                               default=None,
                               help='Name of the shadowPatch')
        self.parser.add_option("--shadowName",
                               action="store",
                               dest="shadowName",
                               default=None,
                               help='Name of the shadowPatch. Deprecated. Use --shadowPatch instead')
        self.parser.add_option("--zone",
                               action="store",
                               dest="zone",
                               default=None,
                               help='Name of the zone for the GGI patch')
        self.parser.add_option("--patchZoneName",
                               action="store",
                               dest="patchZoneName",
                               default=None,
                               help='Name of the zone for the GGI patch. Deprecated. Use --zone instead')
        self.parser.add_option("--bridgeOverlap",
                               action="store",
                               dest="bridgeOverlap",
                               default=None,
                               help='bridgeOverlap flag (on/off)')
        self.parser.add_option("--bridgeOverlapFlag",
                               action="store",
                               dest="bridgeOverlapFlag",
                               default=None,
                               help='bridgeOverlap flag (on/off). Deprecated. Use --bridgeOverlap instead')
        self.parser.add_option("--rotationAxis",
                               action="store",
                               dest="rotationAxis",
                               default=None,
                               help='rotation axis for cyclicGgi')
        self.parser.add_option("--rotationAngle",
                               action="store",
                               dest="rotationAngle",
                               default=None,
                               help='rotation axis angle for cyclicGgi')
        self.parser.add_option("--separationOffset",
                               action="store",
                               dest="separationOffset",
                               default=None,
                               help='separation offset for cyclicGgi')

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
            self.error("Problem with boundary file (not a list)")

        found=False

        for val in bnd:
            if val==bName:
                found=True
            elif found:
                bcType=val["type"]
                if re.match("cyclicGgi", bcType)!= None or re.match("ggi", bcType)!= None:
                    if self.parser.getOptions().shadowPatch!=None:
                        shadowPatch=self.parser.getOptions().shadowPatch
                        val["shadowPatch"]=shadowPatch
                        if shadowPatch not in bnd:
                            self.error("\n    Option --shadowPatch for patch:",bName,": there is no patch called",shadowPatch,"\n")

                    if self.parser.getOptions().zone!=None:
                        val["zone"]=self.parser.getOptions().zone

                    if self.parser.getOptions().bridgeOverlap!=None:
                        val["bridgeOverlap"]=self.parser.getOptions().bridgeOverlap

                    if val["type"]=="cyclicGgi":
                        if self.parser.getOptions().rotationAxis!=None:
                            val["rotationAxis"]=self.parser.getOptions().rotationAxis

                        if self.parser.getOptions().rotationAngle!=None:
                            val["rotationAngle"]=self.parser.getOptions().rotationAngle

                        if self.parser.getOptions().separationOffset!=None:
                            val["separationOffset"]=self.parser.getOptions().separationOffset


                    # Deprecated
                    if self.parser.getOptions().shadowName!=None:
                        self.warning("\n    PatchName:",bName,":  Option --shadowName is deprecated. Use --shadowPatch instead\n")
                        shadowName=self.parser.getOptions().shadowName
                        val["shadowPatch"]=shadowName
                        if shadowName not in bnd:
                            self.error("\n    Option --shadowName for patch:",bName,": there is no patch called",shadowName,"\n")

                    # Deprecated
                    if self.parser.getOptions().patchZoneName!=None:
                        self.warning("\n    PatchName:",bName,":  Option --patchZoneName is deprecated. Use --zone instead\n")
                        val["zone"]=self.parser.getOptions().patchZoneName

                    # Deprecated
                    if self.parser.getOptions().bridgeOverlapFlag!=None:
                        self.warning("\n    PatchName:",bName,":  Option --bridgeOverlapFlag is deprecated. Use --bridgeOverlap instead\n")
                        val["bridgeOverlap"]=self.parser.getOptions().bridgeOverlapFlag


                else:
                    print_("Unsupported GGI type '",bcType,"' for patch",bName)
                break

        if not found:
            self.error("Boundary",bName,"not found in",bnd[::2])

        if self.parser.getOptions().test:
            print_(boundary)
        else:
            boundary.writeFile()
