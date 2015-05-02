"""
Application-class that implements pyFoamInitGgiInterface.py

Inititialize various ggi interface attributes in the
constant/polymesh/boundary file, and in the time directories.

Backups of the boundary file is created.

Generate companion scripts for initializing the ggi zone faceSets.

Modify the decomposeParDict file for the new ggi zones names.

Author:
  Martin Beaudoin, Hydro-Quebec, 2012.  All rights reserved

"""

import sys, fnmatch, re
from os import path, listdir, chmod
from stat import *

from PyFoam.Applications.PyFoamApplication import PyFoamApplication
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.ThirdParty.six import print_
from PyFoam.RunDictionary.TimeDirectory import TimeDirectory
from PyFoam.Basics.BasicFile import BasicFile

class InitGgiInterface(PyFoamApplication):
    def __init__(self,args=None):
        description="""
Init GGI boundary condition parameters in boundary file.
Init GGI boundary fields in time directories.
Generate faceSet scripts for ggi zones.
Modify GGI zones information in decomposeParDict file.
"""
        PyFoamApplication.__init__(self,
                                   args=args,
                                   description=description,
                                   usage="%prog <caseDirectory> ggi_MasterPatchName ggi_ShadowPatchName",
                                   interspersed=True,
                                   changeVersion=False,
                                   nr=3)

    def addOptions(self):
        self.parser.add_option("--type",
                               action="store",
                               dest="ggiType",
                               default='ggi',
                               help='ggi type: ggi | cyclicGgi | overlapGgi')
        self.parser.add_option("--patchZoneName",
                               action="store",
                               dest="patchZoneName",
                               default=None,
                               help='Name of the zone for the GGI patch')
        self.parser.add_option("--bridgeOverlapFlag",
                               action="store",
                               dest="bridgeOverlapFlag",
                               default=None,
                               help='bridgeOverlap flag (on/off)')
        self.parser.add_option("--rotationAxis",
                               action="store",
                               dest="rotationAxis",
                               default=None,
                               help='rotation axis for cyclicGgi or overlapGgi')
        self.parser.add_option("--rotationAngle",
                               action="store",
                               dest="rotationAngle",
                               default=None,
                               help='rotation axis angle for cyclicGgi. Accept Python math expressions like 360.0/17.0')
        self.parser.add_option("--separationOffset",
                               action="store",
                               dest="separationOffset",
                               default=None,
                               help='separation offset for cyclicGgi')
        self.parser.add_option("--nCopies",
                               action="store",
                               dest="nCopies",
                               default=None,
                               help='number of copies for overlapGgi')
        self.parser.add_option("--timeDirs",
                               action="store",
                               dest="timeDirs",
                               default=None,
                               help='time directories for modifying the ggi boundaryfields. Accept expressions like "[0-9]*", "0", etc.')

        self.parser.add_option("--genFaceSetForGgiZonesScriptName",
                               action="store",
                               dest="genFaceSetForGgiZonesScriptName",
                               default="genFaceSetForGgiZones.setSet",
                               help='setSet batch file for generating faceSets for GGI zones. Default: genFaceSetForGgiZones.setSet')

        self.parser.add_option("--initGgiZonesScriptName",
                               action="store",
                               dest="initGgiZonesScriptName",
                               default="initGgiZones.sh",
                               help='script name for initializing the GGI zone faceSets. Default: initGgiZones.sh')

        self.parser.add_option("--test",
                               action="store_true",
                               default=False,
                               dest="test",
                               help="Only print the new boundary file")

    def createGGIPatch(self, patch, patchName, ggiType):
        description="""\
Create a default definition for a ggi patch, and replace
the current definition
        """
        print_("Replacing definition of patch: ", patchName, ":", patch)
        newPatch={
            'type'             : ggiType,
            'nFaces'           : patch["nFaces"],
            'startFace'        : patch["startFace"],
            'shadowPatch'      : 'unknown',
            'zone'             : patchName+'Zone',
            'bridgeOverlap'    : 'true',
            'rotationAxis'     : '(0 0 1)',
            'rotationAngle'    : '0.0',
            'separationOffset' : '(0 0 0)'
            }
        return newPatch

    def modifyGGIPatchDefinition(self, patch, patchName, shadowName, ggiType, rotationAngle):
        description="""\
Modify the definition of a ggi patch
        """
        print_("    Modifying ggi boundary definition in constant/polyMesh/boundary for patch", patchName)

        patch["type"]=ggiType

        patch["shadowPatch"]=shadowName

        if self.parser.getOptions().patchZoneName!=None:
            patch["zone"]=self.parser.getOptions().patchZoneName
        else:
            patch["zone"]=patchName+'Zone'

        if self.parser.getOptions().bridgeOverlapFlag!=None:
            patch["bridgeOverlap"]=self.parser.getOptions().bridgeOverlapFlag

        if ggiType=="cyclicGgi":
            if self.parser.getOptions().rotationAxis!=None:
                patch["rotationAxis"]=self.parser.getOptions().rotationAxis

            # Use the rotationAngle value passed as a parameter.
            # The calling function will change the sign for the slave ggi patch
            if self.parser.getOptions().rotationAngle!=None:
                patch["rotationAngle"]=rotationAngle

            if self.parser.getOptions().separationOffset!=None:
                patch["separationOffset"]=self.parser.getOptions().separationOffset

        if ggiType=="overlapGgi":
            if self.parser.getOptions().rotationAxis!=None:
                patch["rotationAxis"]=self.parser.getOptions().rotationAxis

            if self.parser.getOptions().nCopies!=None:
                patch["nCopies"]=self.parser.getOptions().nCopies


    def modifyGGIPatchDefinitionInTimeDirs(self, caseDir, patchName, ggiType, timeDirs):
        description="""\
Modify the definition of a ggi patch in the time directories
        """
        regex = fnmatch.translate(timeDirs)

        reobj = re.compile(regex)

        for timeDir in listdir(caseDir):
            if reobj.match(timeDir):
                print_("    Modifying ggi boundaryFields in timeDir", timeDir, "for patch", patchName)

                td=TimeDirectory(caseDir, timeDir, yieldParsedFiles=True)

                for f in td:
                    print_("        Modifying field", f.name)
                    f["boundaryField"][patchName]["type"]=ggiType
                    f.writeFile()


    def generateCompanionFiles(self, caseDir, boundary):
        description="""\
Generate a setSet batch file based on the zone info specified in the ggi interfaces definition.
Generate a bash file for invoking setSet and setsToZones
Update GGI zone infoprmation in decomposeParDict
        """
        # Default file: genFaceSetForGgiZones.setSet
        bfGenFaceSets = BasicFile(path.join(caseDir, self.parser.getOptions().genFaceSetForGgiZonesScriptName))

        print_("    Updating file ", bfGenFaceSets.name, " for generating GGI zones faceSet using the setSet command")

        bnd=boundary.content

        if type(bnd)!=list:
            self.error("Problem with boundary file (not a list)")

        # Memorize list of GGI zones for later processing
        listOfGgiZones = []

        for index in range(0, len(bnd), 2):
            patchName = bnd[index]
            indexDefPatch=index+1
            if bnd[indexDefPatch]["type"]=='ggi' or bnd[indexDefPatch]["type"]=='cyclicGgi' or bnd[indexDefPatch]["type"]=='overlapGgi':
                bfGenFaceSets.writeLine([ "faceSet " + bnd[indexDefPatch]["zone"] + " new patchToFace "+ patchName ])
                listOfGgiZones.append(bnd[indexDefPatch]["zone"])

        bfGenFaceSets.writeLine([ "quit" ])
        bfGenFaceSets.close()

        # Default file: initGgiZones.sh
        bfInitGgiZones = BasicFile(path.join(caseDir, self.parser.getOptions().initGgiZonesScriptName))

        print_("    Updating file ", bfInitGgiZones.name, " for inititalizing GGI zones")

        bfInitGgiZones.writeLine([ "#!/bin/bash" ])
        bfInitGgiZones.writeLine([ "setSet -batch " + self.parser.getOptions().genFaceSetForGgiZonesScriptName ])
        bfInitGgiZones.writeLine([ "setsToZones -noFlipMap" ])
        bfInitGgiZones.close()

        # Set execution permissions for this script (755)
        chmod(bfInitGgiZones.name, S_IRWXU|S_IRGRP|S_IXGRP|S_IXOTH|S_IROTH)

        # DecomposeParDict
        decomposeParDictPath=path.join(caseDir,"system","decomposeParDict")
        if path.exists(decomposeParDictPath):
            print_("    Updating file ", decomposeParDictPath, " for GGI zones")
            decomposeParDict=ParsedParameterFile(decomposeParDictPath,debug=False,backup=True)
            dcp=decomposeParDict.content
            dcp["globalFaceZones"]="(\n    " + '\n    '.join(list(listOfGgiZones)) + "\n)"
            decomposeParDict.writeFile()

    def run(self):
        caseDir=self.parser.getArgs()[0]
        masterbName=self.parser.getArgs()[1]
        shadowbName=self.parser.getArgs()[2]

        boundary=ParsedParameterFile(path.join(".",caseDir,"constant","polyMesh","boundary"),debug=False,boundaryDict=True,backup=True)

        bnd=boundary.content

        if type(bnd)!=list:
            self.error("Problem with boundary file (not a list)")

        masterFound=False
        shadowFound=False
        updateTimeDirs=False

        timeDirs="0"
        if self.parser.getOptions().timeDirs!=None:
            timeDirs=self.parser.getOptions().timeDirs
            updateTimeDirs=True

        ggiType=self.parser.getOptions().ggiType

        rotationAngle=0.0
        if self.parser.getOptions().rotationAngle!=None:
            rotationAngle=float(eval(self.parser.getOptions().rotationAngle))

        for index in range(len(bnd)):
            indexDefPatch=index+1

            if bnd[index]==masterbName:
                masterFound=True
                if bnd[indexDefPatch]["type"]!=ggiType:
                    bnd[indexDefPatch] = self.createGGIPatch(bnd[indexDefPatch], masterbName, ggiType)

                self.modifyGGIPatchDefinition(bnd[indexDefPatch], masterbName, shadowbName, ggiType, rotationAngle)

                if updateTimeDirs:
                    self.modifyGGIPatchDefinitionInTimeDirs(caseDir, masterbName, ggiType, timeDirs)

            elif bnd[index]==shadowbName:
                shadowFound=True
                if bnd[indexDefPatch]["type"]!=ggiType:
                    bnd[indexDefPatch] = self.createGGIPatch(bnd[indexDefPatch], shadowbName, ggiType)
                self.modifyGGIPatchDefinition(bnd[indexDefPatch], shadowbName, masterbName, ggiType, -rotationAngle)

                if updateTimeDirs:
                    self.modifyGGIPatchDefinitionInTimeDirs(caseDir, shadowbName, ggiType, timeDirs)

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

        # Write companion files
        self.generateCompanionFiles(caseDir, boundary)
