/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundarySolidBodyMotionFvMesh::
immersedBoundarySolidBodyMotionFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    ibMotions_(),
    curTimeIndex_(-1)
{
    // Read motion function for all regions
    dictionary dynamicMeshCoeffs
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    );

    PtrList<entry> motionDicts(dynamicMeshCoeffs.lookup("motionFunctions"));

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
    // Move points based on new motion
    if (curTimeIndex_ < this->time().timeIndex())
    {
        // Grab old volumes before moving the mesh
        setV0();

        curTimeIndex_ = this->time().timeIndex();
    }

    forAll (ibMotions_, ibI)
    {
        ibMotions_[ibI].movePoints();
    }

    // Force quasi-topological change to rebuild addressing on motion of the
    // immersed boundary
    // HJ, 12/Dec/2017
    fvMesh::syncUpdateMesh();

    // Execute dummy mesh motion for the background mesh
    const pointField oldPoints = allPoints();
    fvMesh::movePoints(oldPoints);

    return true;
}


// ************************************************************************* //
