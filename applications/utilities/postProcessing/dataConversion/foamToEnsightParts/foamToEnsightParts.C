/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    foamToEnsightParts

Description
    Translates OpenFOAM data to Ensight format.
    An Ensight part is created for each cellZone and patch.

Usage
    - foamToEnsightParts [OPTION] \n
    Translates OpenFOAM data to Ensight format

    @param -ascii \n
    Write Ensight data in ASCII format instead of "C Binary"

    @param -zeroTime \n
    Include the often incomplete initial conditions.

Note
    - no parallel data.
    - writes to @a Ensight directory to avoid collisions with foamToEnsight.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "IOobjectList.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ensightParts.H"
#include "ensightOutputFunctions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    // with -constant and -zeroTime
    timeSelector::addOptions(true, false);
    argList::noParallel();
    argList::validOptions.insert("ascii", "");

    const label nTypes = 2;
    const word fieldTypes[] =
    {
        volScalarField::typeName,
        volVectorField::typeName
    };

    const label nSprayFieldTypes = 2;
    const word sprayFieldTypes[] =
    {
        scalarIOField::typeName,
        vectorIOField::typeName
    };

#   include "setRootCase.H"
#   include "createTime.H"

    // get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    // default to binary output, unless otherwise specified
    IOstream::streamFormat format = IOstream::BINARY;
    if (args.options().found("ascii"))
    {
        format = IOstream::ASCII;
    }

    fileName ensightDir = args.rootPath()/args.globalCaseName()/"Ensight";
    fileName dataDir = ensightDir/"data";
    fileName caseFileName = "Ensight.case";
    fileName dataMask = fileName("data")/ensightFile::mask();

    // Ensight and Ensight/data directories must exist
    if (dir(ensightDir))
    {
        rmDir(ensightDir);
    }
    mkDir(ensightDir);
    mkDir(dataDir);

#   include "createMesh.H"
    // Construct the list of ensight parts for the entire mesh
    ensightParts partsList(mesh);

    // write summary information
    {
        OFstream partsInfoFile(ensightDir/"partsInfo");

        partsInfoFile
            << "// summary of ensight parts" << nl << nl;
        partsList.writeSummary(partsInfoFile);
    }

#   include "checkHasMovingMesh.H"
#   include "checkHasLagrangian.H"

    // only take the objects that exists at the end of the calculation
    IOobjectList objects(mesh, timeDirs[timeDirs.size()-1].name());
    IOobjectList sprayObjects(mesh, timeDirs[timeDirs.size()-1].name(), "lagrangian");

    // write single geometry or one per time step
    fileName geometryFileName("geometry");
    if (hasMovingMesh)
    {
        geometryFileName = dataMask/geometryFileName;
    }

    // the case file is always ASCII
    Info << "write case: " << caseFileName.c_str() << endl;

    OFstream caseFile
    (
        ensightDir/caseFileName,
        ios_base::out|ios_base::trunc,
        IOstream::ASCII
    );
    caseFile.setf(ios_base::left);
    caseFile
        << "FORMAT" << nl
        << setw(16) << "type:" << "ensight gold" << nl << nl
        << "GEOMETRY" << nl
        << setw(16) << "model: 1" << geometryFileName.c_str() << nl;

    if (hasLagrangian)
    {
        caseFile
            << setw(16) << "measured: 2"
            << fileName(dataMask/"lagrangian"/"positions").c_str() << nl;
    }
    caseFile
        << nl << "VARIABLE" << nl;

    label nFieldTime = timeDirs.size();
    if (nFieldTime < 0)
    {
        nFieldTime = 0;
    }

    List<label> fieldFileNumbers(nFieldTime);
    List<label> sprayFileNumbers(nFieldTime);

    // map used times used
    Map<scalar>  timeIndices;

    nFieldTime = 0;
    label nSprayTime = 0;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

#       include "getTimeIndex.H"

        fieldFileNumbers[nFieldTime++] = timeIndex;

        // the data/ITER subdirectory must exist
        fileName subDir = ensightFile::subDir(timeIndex);
        mkDir(dataDir/subDir);

        // place a timestamp in the directory for future reference
        {
            OFstream timeStamp(dataDir/subDir/"time");
            timeStamp
                << "#   timestep time" << nl
                << subDir.c_str() << " " << runTime.timeName() << nl;
        }

#       include "moveMesh.H"

        if (nFieldTime == 1 || mesh.moving())
        {
            if (hasMovingMesh)
            {
                geometryFileName = dataDir/subDir/"geometry";
            }
            if (mesh.moving())
            {
                partsList.recalculate(mesh);
            }

            ensightGeoFile geoFile(ensightDir/geometryFileName, format);
            partsList.writeGeometry(geoFile);
            Info << nl;
        }

        Info<< "write volume field: " << flush;

        for (label i=0; i < nTypes; i++)
        {
            wordList fieldNames = objects.names(fieldTypes[i]);

            forAll (fieldNames, fieldI)
            {
                word fieldName = fieldNames[fieldI];

#               include "checkHasValidField.H"

                if (!hasValidField)
                {
                    continue;
                }

                IOobject fieldObject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (fieldTypes[i] == volScalarField::typeName)
                {
                    if (nFieldTime == 1)
                    {
                        ensightCaseEntry<scalar>
                        (
                            caseFile,
                            fieldObject,
                            dataMask
                        );
                    }

                    ensightVolField<scalar>
                    (
                        partsList,
                        fieldObject,
                        mesh,
                        dataDir,
                        subDir,
                        format
                    );

                }
                else if (fieldTypes[i] == volVectorField::typeName)
                {
                    if (nFieldTime == 1)
                    {
                        ensightCaseEntry<vector>
                        (
                            caseFile,
                            fieldObject,
                            dataMask
                        );
                    }

                    ensightVolField<vector>
                    (
                        partsList,
                        fieldObject,
                        mesh,
                        dataDir,
                        subDir,
                        format
                    );

                }
                else if (fieldTypes[i] == volSphericalTensorField::typeName)
                {
                    if (nFieldTime == 1)
                    {
                        ensightCaseEntry<sphericalTensor>
                        (
                            caseFile,
                            fieldObject,
                            dataMask
                        );
                    }

                    ensightVolField<sphericalTensor>
                    (
                        partsList,
                        fieldObject,
                        mesh,
                        dataDir,
                        subDir,
                        format
                    );

                }
                else if (fieldTypes[i] == volSymmTensorField::typeName)
                {
                    if (nFieldTime == 1)
                    {
                        ensightCaseEntry<symmTensor>
                        (
                            caseFile,
                            fieldObject,
                            dataMask
                        );
                    }

                    ensightVolField<symmTensor>
                    (
                        partsList,
                        fieldObject,
                        mesh,
                        dataDir,
                        subDir,
                        format
                    );

                }
                else if (fieldTypes[i] == volTensorField::typeName)
                {
                    if (nFieldTime == 1)
                    {
                        ensightCaseEntry<tensor>
                        (
                            caseFile,
                            fieldObject,
                            dataMask
                        );
                    }

                    ensightVolField<tensor>
                    (
                        partsList,
                        fieldObject,
                        mesh,
                        dataDir,
                        subDir,
                        format
                    );

                }
            }
        }
        Info<< endl;


        if (hasLagrangian)
        {
            // check that the positions field is present for this time
            {
                IOobject ioHeader
                (
                    "positions",
                    mesh.time().timeName(),
                    "lagrangian",
                    mesh,
                    IOobject::NO_READ
                );

                if (ioHeader.headerOk())
                {
                    sprayFileNumbers[nSprayTime++] = timeIndex;
                }
            }

            Info<< "write  spray field: " << flush;

            ensightParticlePositions
            (
                mesh,
                dataDir,
                subDir,
                format
            );

            for (label i=0; i < nSprayFieldTypes; i++)
            {
                wordList fieldNames = sprayObjects.names(sprayFieldTypes[i]);

                forAll (fieldNames, fieldI)
                {
                    word fieldName = fieldNames[fieldI];

#                   include "checkHasSprayField.H"

                    if (!hasSprayField)
                    {
                        continue;
                    }

                    IOobject fieldObject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        "lagrangian",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (sprayFieldTypes[i] == scalarIOField::typeName)
                    {
                        if (nSprayTime == 1)
                        {
                            ensightCaseEntry<scalar>
                            (
                                caseFile,
                                fieldObject,
                                dataMask,
                                true
                            );
                        }

                        ensightSprayField<scalar>
                        (
                            fieldObject,
                            dataDir,
                            subDir,
                            format
                        );

                    }
                    else if (sprayFieldTypes[i] == vectorIOField::typeName)
                    {
                        if (nSprayTime == 1)
                        {
                            ensightCaseEntry<vector>
                            (
                                caseFile,
                                fieldObject,
                                dataMask,
                                true
                            );
                        }

                        ensightSprayField<vector>
                        (
                            fieldObject,
                            dataDir,
                            subDir,
                            format
                        );

                    }
                    else if (sprayFieldTypes[i] == tensorIOField::typeName)
                    {
                        if (nSprayTime == 1)
                        {
                            ensightCaseEntry<tensor>
                            (
                                caseFile,
                                fieldObject,
                                dataMask,
                                true
                            );
                        }

                        ensightSprayField<tensor>
                        (
                            fieldObject,
                            dataDir,
                            subDir,
                            format
                        );

                    }
                }
            }
            Info<< endl;
        }
    }

    fieldFileNumbers.setSize(nFieldTime);
    sprayFileNumbers.setSize(nSprayTime);

    // add time values
    caseFile << nl << "TIME" << nl;
#   include "ensightCaseTimes.H"

    Info<< "\nEnd\n"<< endl;

    return 0;
}


// ************************************************************************* //
