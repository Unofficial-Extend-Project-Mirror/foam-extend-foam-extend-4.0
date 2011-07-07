/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM

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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "turboFvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turboFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, turboFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turboFvMesh::calcMovingPoints() const
{
    if (debug)
    {
        Info<< "void turboFvMesh::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsPtr_)
    {
        FatalErrorIn("void turboFvMesh::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Retrieve the cell zone Names
    const wordList cellZoneNames = cellZones().names();

    // Set the points
    movingPointsPtr_ = new vectorField(allPoints().size(), vector::zero);

    vectorField& movingPoints = *movingPointsPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    scalar rpm;

    forAll (cellZoneNames,cellZoneI)
    {
        const labelList& cellAddr =
            cellZones()[cellZones().findZoneID(cellZoneNames[cellZoneI])];

        if (dict_.subDict("rpm").found(cellZoneNames[cellZoneI]))
        {
            rpm = readScalar
            (
                dict_.subDict("rpm").lookup(cellZoneNames[cellZoneI])
            );

            Info<< "Moving Cell Zone Name: " << cellZoneNames[cellZoneI]
                << " rpm: " << rpm << endl;

            forAll (cellAddr, cellI)
            {
                const cell& curCell = c[cellAddr[cellI]];

                forAll (curCell, faceI)
                {
                    const face& curFace = f[curCell[faceI]];

                    forAll (curFace, pointI)
                    {
                        // The rotation data is saved within the cell data. For
                        // non-rotating regions rpm is zero, so mesh movement
                        // is also zero. The conversion of rotational speed

                        // Note: deltaT changes during the run: moved to
                        // turboFvMesh::update().  HJ, 14/Oct/2010
                        movingPoints[curFace[pointI]] =
                            vector(0, rpm/60.0*360.0, 0);
                    }
                }
            }
        }
    }

    // Retrieve the face zone Names
    const wordList faceZoneNames = faceZones().names();

    // This is not bullet proof, as one could specify the wrong name of a
    // faceZone, which is then not selected. The solver will crash after the
    // first meshMotion iteration.
    forAll (faceZoneNames, faceZoneI)
    {
        if (dict_.subDict("slider").found(faceZoneNames[faceZoneI]))
        {
            rpm = readScalar
            (
                dict_.subDict("slider").lookup(faceZoneNames[faceZoneI])
            );

            Info<< "Moving Face Zone Name: " << faceZoneNames[faceZoneI]
                << " rpm: " << rpm << endl;

            faceZoneID zoneID(faceZoneNames[faceZoneI], faceZones());

            const labelList& movingSliderAddr = faceZones()[zoneID.index()];

            forAll (movingSliderAddr, faceI)
            {
                const face& curFace = f[movingSliderAddr[faceI]];

                forAll (curFace, pointI)
                {
                    movingPoints[curFace[pointI]] =
                        vector( 0, rpm/60.0*360.0, 0);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::turboFvMesh::turboFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            dict_.subDict("coordinateSystem")
        )
    ),
    movingPointsPtr_(NULL)
{
    Info<< "Turbomachine Mixer mesh:" << nl
        << "    origin: " << cs().origin() << nl
        << "    axis  : " << cs().axis() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turboFvMesh::~turboFvMesh()
{
    deleteDemandDrivenData(movingPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask
const Foam::vectorField& Foam::turboFvMesh::movingPoints() const
{
    if (!movingPointsPtr_)
    {
        calcMovingPoints();
    }

    return *movingPointsPtr_;
}


bool Foam::turboFvMesh::update()
{
    movePoints
    (
        csPtr_->globalPosition
        (
            csPtr_->localPosition(allPoints())
          + movingPoints()*time().deltaT().value()
        )
    );

    // The mesh is not morphing
    return false;
}


// ************************************************************************* //
