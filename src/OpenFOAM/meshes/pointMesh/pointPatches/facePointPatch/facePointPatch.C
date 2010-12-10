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

#include "facePointPatch.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "demandDrivenData.H"
#include "boolList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(facePointPatch, 0);
defineRunTimeSelectionTable(facePointPatch, polyPatch);

addToRunTimeSelectionTable
(
    facePointPatch,
    facePointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void facePointPatch::initGeometry()
{}


void facePointPatch::calcGeometry()
{}


void facePointPatch::initMovePoints(const pointField&)
{}


void facePointPatch::movePoints(const pointField&)
{}


void facePointPatch::initUpdateMesh()
{
    facePointPatch::initGeometry();
}


void facePointPatch::updateMesh()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

facePointPatch::facePointPatch
(
    const polyPatch& p,
    const pointBoundaryMesh& bm
)
:
    pointPatch(bm),
    polyPatch_(p)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
