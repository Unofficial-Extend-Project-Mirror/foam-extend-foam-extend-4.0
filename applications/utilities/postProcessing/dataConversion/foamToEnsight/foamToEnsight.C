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

Description
    Translates FOAM data to EnSight format

Usage
    - foamToEnsight [OPTION] \n
    Translates OpenFOAM data to Ensight format

    @param -ascii \n
    Write Ensight data in ASCII format instead of "C Binary"

Note
    Parallel support for cloud data is not supported

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"

#include "volFields.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ensightMesh.H"
#include "ensightField.H"

#include "ensightParticlePositions.H"
#include "ensightCloudField.H"

#include "fvc.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool inFileNameList
(
    const fileNameList& nameList,
    const word& name
)
{
    forAll(nameList, i)
    {
        if (nameList[i] == name)
        {
            return true;
        }
    }

    return false;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("patches", "patch list");
    argList::validOptions.insert("ascii", "" );
#   include "addTimeOptions.H"

#   include "setRootCase.H"

    // Check options
    bool binary = true;
    if (args.options().found("ascii"))
    {
        binary = false;
    }

#   include "createTime.H"

    // get the available time-steps
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createNamedMesh.H"

    // Mesh instance (region0 gets filtered out)
    fileName regionPrefix = "";

    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    const label nTypes = 2;
    const word fieldTypes[] =
    {
        volScalarField::typeName,
        volVectorField::typeName
    };

    // Create the output folder
    const word postProcDir = "EnSight";

    // Path to EnSight folder at case level only
    // - For parallel cases, data only written from master
//    fileName postProcPath = runTime.path()/postProcDir;
    fileName postProcPath = args.rootPath()/args.globalCaseName()/postProcDir;

    if (Pstream::master())
    {
        if (dir(postProcPath))
        {
            rmDir(postProcPath);
        }

        mkDir(postProcPath);
    }

    // Start of case file header output
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const word prepend = args.globalCaseName() + '.';

    OFstream *ensightCaseFilePtr = NULL;
    if (Pstream::master())
    {
        fileName ensightCaseFileName = prepend + "case";

        if (!binary)
        {
            ensightCaseFilePtr = new OFstream
            (
                postProcPath/ensightCaseFileName,
                ios_base::out|ios_base::trunc,
                runTime.writeFormat(),
                runTime.writeVersion(),
                runTime.writeCompression()
            );
        }
        else
        {
            ensightCaseFilePtr = new OFstream
            (
                postProcPath/ensightCaseFileName,
                ios_base::out|ios_base::trunc,
                runTime.writeFormat(),
                runTime.writeVersion(),
                IOstream::UNCOMPRESSED
            );
        }

        Info<< nl << "Case file is " << ensightCaseFileName << endl;
    }

    OFstream& ensightCaseFile = *ensightCaseFilePtr;

#   include "ensightCaseHeader.H"

    // Construct the EnSight mesh
    ensightMesh eMesh(mesh, args, binary);

    // Set Time to the last time before looking for the lagrangian objects
    runTime.setTime(Times[Times.size()-1], Times.size()-1);

    IOobjectList objects(mesh, runTime.timeName());

#   include "checkMeshMoving.H"

    wordHashSet allCloudNames;
    word geomCaseFileName = prepend + "000";
    if (Pstream::master())
    {
        // test pre check variable if there is a moving mesh
        if (meshMoving == true)
        {
            geomCaseFileName = prepend + "***";
        }

        ensightCaseFile
            << "GEOMETRY" << nl
            << "model:        1     "
            << (geomCaseFileName + ".mesh").c_str() << nl;
    }

    // Identify if lagrangian data exists at each time, and add clouds
    // to the 'allCloudNames' hash set
    for (label n=startTime; n<endTime; n++)
    {
        runTime.setTime(Times[n], n);

        fileNameList cloudDirs = readDir
        (
            runTime.timePath()/regionPrefix/"lagrangian",
            fileName::DIRECTORY
        );

        forAll(cloudDirs, cloudI)
        {
            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                "lagrangian"/cloudDirs[cloudI]
            );

            IOobject* positionsPtr = cloudObjs.lookup("positions");

            if (positionsPtr)
            {
                allCloudNames.insert(cloudDirs[cloudI]);
            }
        }
    }

    HashTable<HashTable<word> > allCloudFields;
    forAllConstIter(wordHashSet, allCloudNames, cloudIter)
    {
        // Add the name of the cloud(s) to the case file header
        if (Pstream::master())
        {
            ensightCaseFile
            <<  (
                    "measured:     1     "
                  + prepend
                  + "***."
                  + cloudIter.key()
                ).c_str()
            << nl;
        }

        // Create a new hash table for each cloud
        allCloudFields.insert(cloudIter.key(), HashTable<word>());

        // Identify the new cloud in the hash table
        HashTable<HashTable<word> >::iterator newCloudIter =
            allCloudFields.find(cloudIter.key());

        // Loop over all times to build list of fields and field types
        // for each cloud
        for (label n=startTime; n<endTime; n++)
        {
            runTime.setTime(Times[n], n);

            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                "lagrangian"/cloudIter.key()
            );

            forAllConstIter(IOobjectList, cloudObjs, fieldIter)
            {
                const IOobject obj = *fieldIter();

                if (obj.name() != "positions")
                {
                    // Add field and field type
                    newCloudIter().insert
                    (
                        obj.name(),
                        obj.headerClassName()
                    );
                }
            }
        }
    }

    label nTimeSteps = 0;
    for (label n=startTime; n<endTime; n++)
    {
        nTimeSteps++;
        runTime.setTime(Times[n], n);
        label timeIndex = n - startTime;

        word timeName = itoa(timeIndex);
        word timeFile = prepend + timeName;

        Info<< "Translating time = " << runTime.timeName() << nl;

#       include "moveMesh.H"

        if (timeIndex == 0 || mesh.moving())
        {
            eMesh.write
            (
                postProcPath,
                prepend,
                timeIndex,
                ensightCaseFile
            );
        }


        // Start of field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (timeIndex == 0 && Pstream::master())
        {
            ensightCaseFile<< nl << "VARIABLE" << nl;
        }


        // Cell field data output
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label i=0; i<nTypes; i++)
        {
            wordList fieldNames = objects.names(fieldTypes[i]);

            for (label j=0; j<fieldNames.size(); j++)
            {
                word fieldName = fieldNames[j];

                bool variableGood = true;

#               include "checkData.H"

                if (variableGood)
                {
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
                        ensightField<scalar>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volVectorField::typeName)
                    {
                        ensightField<vector>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volSphericalTensorField::typeName)
                    {
                        ensightField<sphericalTensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volSymmTensorField::typeName)
                    {
                        ensightField<symmTensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volTensorField::typeName)
                    {
                        ensightField<tensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                }
            }
        }


        // Cloud field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~

        forAllConstIter(HashTable<HashTable<word> >, allCloudFields, cloudIter)
        {
            const word& cloudName = cloudIter.key();

            fileNameList currentCloudDirs = readDir
            (
                runTime.timePath()/regionPrefix/"lagrangian",
                fileName::DIRECTORY
            );

            bool cloudExists = inFileNameList(currentCloudDirs, cloudName);
            ensightParticlePositions
            (
                mesh,
                postProcPath,
                timeFile,
                cloudName,
                cloudExists
            );

            forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
            {
                const word& fieldName = fieldIter.key();
                const word& fieldType = fieldIter();

                IOobject fieldObject
                (
                    fieldName,
                    mesh.time().timeName(),
                    "lagrangian"/cloudName,
                    mesh,
                    IOobject::MUST_READ
                );

                bool fieldExists = fieldObject.headerOk();
                if (fieldType == scalarIOField::typeName)
                {
                    ensightCloudField<scalar>
                    (
                        fieldObject,
                        postProcPath,
                        prepend,
                        timeIndex,
                        cloudName,
                        ensightCaseFile,
                        fieldExists
                    );
                }
                else if (fieldType == vectorIOField::typeName)
                {
                    ensightCloudField<vector>
                    (
                        fieldObject,
                        postProcPath,
                        prepend,
                        timeIndex,
                        cloudName,
                        ensightCaseFile,
                        fieldExists
                    );
                }
                else
                {
                    Info<< "Unable to convert field type " << fieldType
                        << " for field " << fieldName << endl;
                }
            }
        }
    }

#   include "ensightCaseTail.H"

    if (Pstream::master())
    {
        delete ensightCaseFilePtr;
    }

    return 0;
}


// ************************************************************************* //
