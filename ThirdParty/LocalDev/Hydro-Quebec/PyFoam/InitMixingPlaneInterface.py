"""
Application-class that implements pyFoamInitMixingPlaneInterface.py

Initialize various mixingPlane interface attributes in the
constant/polymesh/boundary file, and in the time directories.

Backups of the modified files are created

Author:
  Martin Beaudoin, Hydro-Quebec, 2012.  All rights reserved

"""

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.TimeDirectory import TimeDirectory
from PyFoam.ThirdParty.six import print_
from os import path, listdir
import sys, fnmatch, re

class InitMixingPlaneInterface(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Init MixingPlane boundary condition parameters
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory> mixingPlane_MasterPatchName mixingPlane_ShadowPatchName",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=3)

    def addOptions(self):
        # No need for option --shadowPatch since the name of the shadowPatch is a mandatory parameter to the script
        #self.parser.add_option("--shadowPatch",
        #                       action="store",
        #                       dest="shadowPatch",
        #                       default=None,
        #                       help='Name of the shadowPatch')
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

        self.parser.add_option("--timeDirs",
                               action="store",
                               dest="timeDirs",
                               default=None,
                               help='time directories for the mixingPlane boundaryfields. Accept expressions like "[0-9]*", "0", etc.')

        self.parser.add_option("--test",
                               action="store_true",
                               default=False,
                               dest="test",
                               help="Only print the new boundary file")

    def createMixingPlanePatch(self, patch, patchName):
        description="""\
Create a default definition for a mixingPlane patch, and replace
the current definition
        """
        print_("Replacing definition of patch: ", patchName, ":", patch)
        newPatch={
            'type'        : "mixingPlane",
            'nFaces'      : patch["nFaces"],
            'startFace'   : patch["startFace"],
            'shadowPatch' : 'unknown',
            'coordinateSystem' : {
                'name'   : 'mixingCS',
                'type'   : 'cylindrical',
                'origin' : '(0 0 0)',
                'e1'     : '(1 0 0)',
                'e3'     : '(0 0 1)'
                },
            'ribbonPatch' : {
                'sweepAxis'      : 'Theta',
                'stackAxis'      : 'Z',
                'discretisation' : 'bothPatches',
                }
            }
        return newPatch

    def modifyMixinPlanePatchDefinition(self, patch, patchName, shadowName):
        description="""\
Modify the definition of a mixingPlane patch
        """
        print_("    Modifying mixingPlane boundary definition in constant/polyMesh/boundary for patch", patchName)

        patch["shadowPatch"]=shadowName

        if patch.has_key("coordinateSystem")==False:
            patch["coordinateSystem"]={}

        if self.parser.getOptions().coordinateSystemName!=None:
            patch["coordinateSystem"]["name"]=self.parser.getOptions().coordinateSystemName

        if self.parser.getOptions().coordinateSystemType!=None:
            patch["coordinateSystem"]["type"]=self.parser.getOptions().coordinateSystemType

        if self.parser.getOptions().coordinateSystemOrigin!=None:
            patch["coordinateSystem"]["origin"]=self.parser.getOptions().coordinateSystemOrigin

        if self.parser.getOptions().coordinateSystemE1!=None:
            patch["coordinateSystem"]["e1"]=self.parser.getOptions().coordinateSystemE1

        if self.parser.getOptions().coordinateSystemE3!=None:
            patch["coordinateSystem"]["e3"]=self.parser.getOptions().coordinateSystemE3

        if patch.has_key("ribbonPatch")==False:
            patch["ribbonPatch"]={}

        if self.parser.getOptions().ribbonPatchSweepAxis!=None:
            patch["ribbonPatch"]["sweepAxis"]=self.parser.getOptions().ribbonPatchSweepAxis

        if self.parser.getOptions().ribbonPatchStackAxis!=None:
            patch["ribbonPatch"]["stackAxis"]=self.parser.getOptions().ribbonPatchStackAxis

        if self.parser.getOptions().ribbonPatchDiscretisation!=None:
            patch["ribbonPatch"]["discretisation"]=self.parser.getOptions().ribbonPatchDiscretisation



    def modifyMixinPlanePatchDefinitionInTimeDirs(self, caseDir, patchName, timeDirs):
        description="""\
Modify the definition of a mixingPlane patch in the time directories
        """
        regex = fnmatch.translate(timeDirs)

        reobj = re.compile(regex)

        for timeDir in listdir(caseDir):
            if reobj.match(timeDir):
                print_("    Modifying mixingPlane boundaryFields in timeDir", timeDir, "for patch", patchName)

                td=TimeDirectory(caseDir, timeDir, yieldParsedFiles=True)

                for f in td:
                    print_("        Modifying field", f.name)
                    f["boundaryField"][patchName]["type"]='mixingPlane'
                    f.writeFile()

    def run(self):
        fName=self.parser.getArgs()[0]
        masterbName=self.parser.getArgs()[1]
        shadowbName=self.parser.getArgs()[2]

        boundary=ParsedParameterFile(path.join(".",fName,"constant","polyMesh","boundary"),debug=False,boundaryDict=True,backup=True)

        bnd=boundary.content

        if type(bnd)!=list:
            print_("Problem with boundary file (not a list)")
            sys.exit(-1)

        masterFound=False
        shadowFound=False
        updateTimeDirs=False

        timeDirs="0"
        if self.parser.getOptions().timeDirs!=None:
            timeDirs=self.parser.getOptions().timeDirs
            updateTimeDirs=True

        print_("UpdateTimeDirs: ", updateTimeDirs)

        for index in range(len(bnd)):

            indexDefPatch=index+1

            if bnd[index]==masterbName:
                masterFound=True
                if bnd[indexDefPatch]["type"]!="mixingPlane":
                    bnd[indexDefPatch] = self.createMixingPlanePatch(bnd[indexDefPatch], masterbName)

                self.modifyMixinPlanePatchDefinition(bnd[indexDefPatch], masterbName, shadowbName)

                if updateTimeDirs:
                    self.modifyMixinPlanePatchDefinitionInTimeDirs(fName, masterbName, timeDirs)

            elif bnd[index]==shadowbName:
                shadowFound=True
                if bnd[indexDefPatch]["type"]!="mixingPlane":
                    bnd[indexDefPatch] = self.createMixingPlanePatch(bnd[indexDefPatch], shadowbName)

                self.modifyMixinPlanePatchDefinition(bnd[indexDefPatch], shadowbName, masterbName)

                if updateTimeDirs:
                    self.modifyMixinPlanePatchDefinitionInTimeDirs(fName, shadowbName, timeDirs)

            if masterFound and shadowFound:
                break;

        if not masterFound:
            self.error("Boundary patch",masterbName,"not found in",bnd[::2])

        if not shadowFound:
            self.error("Boundary patch",shadowbName,"not found in",bnd[::2])

        if self.parser.getOptions().test:
            print_(boundary)
        else:
            boundary.writeFile()
