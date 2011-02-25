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

#include "coupledFaceTetPolyPatchCellDecomp.H"
#include "tetPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledFaceTetPolyPatchCellDecomp, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledFaceTetPolyPatchCellDecomp::coupledFaceTetPolyPatchCellDecomp
(
    const polyPatch& patch,
    const tetPolyBoundaryMeshCellDecomp& bm
)
:
    faceTetPolyPatchCellDecomp(patch, bm),
    coupledPolyPatch_(refCast<const coupledPolyPatch>(patch)),
    nonGlobalPatchPointsPtr_(NULL),
    meshPointsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledFaceTetPolyPatchCellDecomp::~coupledFaceTetPolyPatchCellDecomp()
{
    deleteDemandDrivenData(nonGlobalPatchPointsPtr_);
    deleteDemandDrivenData(meshPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList&
Foam::coupledFaceTetPolyPatchCellDecomp::nonGlobalPatchPoints() const
{
    if (!nonGlobalPatchPointsPtr_)
    {
        calcMeshPoints();
    }

    return *nonGlobalPatchPointsPtr_;
}


const Foam::labelList&
Foam::coupledFaceTetPolyPatchCellDecomp::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshPoints();
    }

    return *meshPointsPtr_;
}


const Foam::pointField&
Foam::coupledFaceTetPolyPatchCellDecomp::localPoints() const
{
    notImplemented("coupledFaceTetPolyPatchCellDecomp::localPoints() const");
    return pointField::null();
}


const Foam::vectorField&
Foam::coupledFaceTetPolyPatchCellDecomp::pointNormals() const
{
    notImplemented("coupledFaceTetPolyPatchCellDecomp::pointNormals() const");
    return vectorField::null();
}


// ************************************************************************* //
