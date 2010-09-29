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
    PODecomposition

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Calculates proper orthogonal decomposition of a given field set

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PODOrthoNormalBases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    const label startTime = 1;
    const label endTime = Times.size();
    const label nSnapshots = Times.size() - 1;

    Info << "Number of snapshots: " << nSnapshots << endl;

    // Create a list of snapshots
    PtrList<volScalarField> fields(nSnapshots);

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    IOdictionary PODsolverDict
    (
        IOobject
        (
            "PODsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar accuracy =
        readScalar
        (
            PODsolverDict.subDict("scalarTransportCoeffs").lookup("accuracy")
        );

    Info << "Seeking accuracy: " << accuracy << endl;

    word fieldName
    (
        PODsolverDict.subDict("scalarTransportCoeffs").lookup("field")
    );

    label snapI = 0;

    labelList timeIndices(nSnapshots);

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        Info<< "    Reading " << fieldName << endl;
        fields.set
        (
            snapI,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
            )
        );

        // Rename the field
        fields[snapI].rename(fieldName + name(i));
        timeIndices[snapI] = i;
        snapI++;

        Info<< endl;
    }

    timeIndices.setSize(snapI);

    // Read accurary
    Info<< "Reading \n" << endl;

    scalarPODOrthoNormalBase eb(fields, accuracy);

    const scalarRectangularMatrix& coeffs = eb.interpolationCoeffs();

    // Check all snapshots
    forAll (fields, fieldI)
    {
        runTime.setTime(Times[timeIndices[fieldI]], timeIndices[fieldI]);

        volScalarField pReconstruct
        (
            IOobject
            (
                fieldName + "PODreconstruct",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("zero", fields[fieldI].dimensions(), 0)
        );

        for (label baseI = 0; baseI < eb.baseSize(); baseI++)
        {
            pReconstruct +=
                coeffs[fieldI][baseI]*eb.orthoField(baseI);
        }

        scalar sumFieldError =
            Foam::sqrt
            (
                sumSqr
                (
                    pReconstruct.internalField()
                  - fields[fieldI].internalField()
                )
            );

        scalar measure =
            Foam::sqrt(sumSqr(fields[fieldI].internalField())) + SMALL;

        scalar sumFieldRelError = sumFieldError/measure;

        Info<< "Field error: absolute = " << sumFieldError
            << " relative = " << sumFieldRelError 
            << " measure = " << measure << endl;

        pReconstruct.write();
    }

    // Write out all fields
    runTime.setTime(Times[startTime], startTime);
    Info<< "Writing POD base for Time = " << runTime.timeName() << endl;

    for (label i = 0; i < eb.baseSize(); i++)
    {
        eb.orthoField(i).write();
    }

    Info << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
