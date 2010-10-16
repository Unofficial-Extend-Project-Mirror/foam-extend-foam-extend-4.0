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

\*---------------------------------------------------------------------------*/

#include "movingBoxFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(movingBoxFvMesh, 0);

addToRunTimeSelectionTable(dynamicFvMesh, movingBoxFvMesh, IOobject);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingBoxFvMesh::movingBoxFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    movingMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    splitDirection_(movingMeshCoeffs_.lookup("splitDirection")),
    leftEdge_(movingMeshCoeffs_.lookup("leftEdge")),
    rightEdge_(movingMeshCoeffs_.lookup("rightEdge")),
    velocity_(movingMeshCoeffs_.lookup("velocity")),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    motionMarkup_(stationaryPoints_.size(), 0)
{
    if (mag(splitDirection_) < SMALL)
    {
        FatalErrorIn("movingBoxFvMesh::movingBoxFvMesh(const IOobject& io)")
            << "Incorrect definition of split direction"
            << abort(FatalError);
    }

    splitDirection_ /= mag(splitDirection_);

    Info<< "Performing a moving mesh calculation: " << nl
        << "splitDirection: " << splitDirection_ << nl
        << "leftEdge: " << leftEdge_ << nl
        << "rightEdge: " << rightEdge_ << nl
        << "velocity: " << velocity_ << endl;

    // Calculate vertex markup
    scalarField p = points() & splitDirection_;
    scalar leftP = (leftEdge_ & splitDirection_);
    scalar rightP = (rightEdge_ & splitDirection_);

    // Re-scale p to be between 0 and 1

    scalar minP = min(p);
    p -= minP;
    leftP -= minP;
    rightP -= minP;

    scalar maxP = max(p);
    p /= maxP - minP;
    leftP /= maxP - minP;
    rightP /= maxP - minP;

    Info << "leftP: " << leftP << " rightP: " << rightP << " min: " << min(p)
         << " max: " << max(p);

    // Left part
    motionMarkup_ += neg(p - leftP - SMALL)*p/leftP;

    // Centre part
    motionMarkup_ += neg(p - rightP + SMALL)*pos(p - leftP - SMALL);

    // Right part
    motionMarkup_ += pos(p - rightP + SMALL)*(1.0 - p)/(1.0 - rightP);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

movingBoxFvMesh::~movingBoxFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool movingBoxFvMesh::update()
{
    pointField newPoints =
        stationaryPoints_ + motionMarkup_*velocity_*time().value();

    fvMesh::movePoints(newPoints);

    // Mesh motion only - return false
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
