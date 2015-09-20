/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "coupledFaceTetPolyPatch.H"
#include "tetPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledFaceTetPolyPatch, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledFaceTetPolyPatch::coupledFaceTetPolyPatch
(
    const polyPatch& patch,
    const tetPolyBoundaryMesh& bm
)
:
    faceTetPolyPatch(patch, bm),
    coupledPolyPatch_(refCast<const coupledPolyPatch>(patch)),
    nonGlobalPatchPointsPtr_(NULL),
    meshPointsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledFaceTetPolyPatch::~coupledFaceTetPolyPatch()
{
    deleteDemandDrivenData(nonGlobalPatchPointsPtr_);
    deleteDemandDrivenData(meshPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList&
Foam::coupledFaceTetPolyPatch::nonGlobalPatchPoints() const
{
    if (!nonGlobalPatchPointsPtr_)
    {
        calcMeshPoints();
    }

    return *nonGlobalPatchPointsPtr_;
}


const Foam::labelList&
Foam::coupledFaceTetPolyPatch::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshPoints();
    }

    return *meshPointsPtr_;
}


const Foam::pointField&
Foam::coupledFaceTetPolyPatch::localPoints() const
{
    notImplemented("coupledFaceTetPolyPatch::localPoints() const");
    return pointField::null();
}


const Foam::vectorField&
Foam::coupledFaceTetPolyPatch::pointNormals() const
{
    notImplemented("coupledFaceTetPolyPatch::pointNormals() const");
    return vectorField::null();
}


// ************************************************************************* //
