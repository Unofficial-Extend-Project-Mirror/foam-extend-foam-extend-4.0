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

#include "oversetSolidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetSolidBodyMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        oversetSolidBodyMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetSolidBodyMotionFvMesh::oversetSolidBodyMotionFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
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
    motionRegions_(),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    curTimeIndex_(-1)
{
    // Read motion function for all regions
    PtrList<entry> motionDicts(dynamicMeshCoeffs_.lookup("motionFunctions"));

    motionRegions_.setSize(motionDicts.size());

    forAll (motionDicts, mI)
    {
        motionRegions_.set
        (
            mI,
            new movingOversetRegion
            (
                motionDicts[mI].keyword(),
                *this,
                motionDicts[mI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetSolidBodyMotionFvMesh::~oversetSolidBodyMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::oversetSolidBodyMotionFvMesh::update()
{
    // Move the mesh only once per time step as this is prescribed motion
    if (curTimeIndex_ < this->time().timeIndex())
    {
        pointField newPoints = undisplacedPoints_;

        forAll (motionRegions_, regionI)
        {
            newPoints +=
                motionRegions_[regionI].motionIncrement(undisplacedPoints_);
        }

        fvMesh::movePoints(newPoints);

        curTimeIndex_ = this->time().timeIndex();
    }

    return false;
}


// ************************************************************************* //
