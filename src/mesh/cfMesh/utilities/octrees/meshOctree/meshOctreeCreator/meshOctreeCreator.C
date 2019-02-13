/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    meshDictPtr_(nullptr),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size())
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
    surfRefLevel_(mo.surface().size())
{}

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
