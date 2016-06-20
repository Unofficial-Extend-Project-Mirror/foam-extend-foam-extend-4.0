/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "immersedBoundarySolidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundarySolidBodyMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        immersedBoundarySolidBodyMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundarySolidBodyMotionFvMesh::
immersedBoundarySolidBodyMotionFvMesh
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
    ibMotions_()
{
    // Read motion function for all regions
    PtrList<entry> motionDicts(dynamicMeshCoeffs_.lookup("motionFunctions"));

    ibMotions_.setSize(motionDicts.size());

    forAll (motionDicts, mI)
    {
        ibMotions_.set
        (
            mI,
            new movingImmersedBoundary
            (
                motionDicts[mI].keyword(),
                *this,
                motionDicts[mI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundarySolidBodyMotionFvMesh::
~immersedBoundarySolidBodyMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::immersedBoundarySolidBodyMotionFvMesh::update()
{
    forAll (ibMotions_, ibI)
    {
        ibMotions_[ibI].movePoints();
    }

    // Force flux and addressing recalculation as in topo change
    pointField newAllPoints = allPoints();

    movePoints(newAllPoints);

    return true;
}


// ************************************************************************* //
