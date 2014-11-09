/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "boundBox.H"
#include "demandDrivenData.H"


// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshOctreeCreator::meshOctreeCreator(meshOctree& mo)
:
    octree_(mo),
    scalingFactor_(1.0),
    meshDictPtr_(NULL),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0)),
    surfRefThickness_(mo.surface().size(), 0.0)
{}

meshOctreeCreator::meshOctreeCreator
(
    meshOctree& mo,
    const IOdictionary& dict
)
:
    octree_(mo),
    scalingFactor_(1.0),
    meshDictPtr_(&dict),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0)),
    surfRefThickness_(mo.surface().size(), 0.0)
{}

/*
meshOctreeCreator::meshOctreeCreator
(
    meshOctree& mo,
    const IOdictionary& dict,
    const volScalarField& localCellSize
)
:
    octree_(mo),
    scalingFactor_(1.0),
    meshDict_(dict),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size(), direction(0)),
    surfRefThickness_(mo.surface().size(), 0.0)
{
    FatalErrorIn
    (
        "meshOctreeCreator::meshOctreeCreator"
        "("
        "meshOctree& mo,"
        "const IOdictionary& dict,"
        "const volScalarField& localCellSize"
        ")"
    ) << "Not implemented!" << exit(FatalError);

    Info << "Constructing octree" << endl;

    Info << "Finished constructing octree" << endl;
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctreeCreator::~meshOctreeCreator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeCreator::setScalingFactor(const scalar s)
{
    scalingFactor_ = s;
}

void meshOctreeCreator::activateHexRefinement()
{
    hexRefinement_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
