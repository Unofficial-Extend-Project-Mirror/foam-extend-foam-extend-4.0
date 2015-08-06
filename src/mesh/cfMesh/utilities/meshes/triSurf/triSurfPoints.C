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

#include "triSurfPoints.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfPoints::triSurfPoints()
:
    points_(),
    pointSubsets_()
{}

triSurfPoints::triSurfPoints(const pointField& points)
:
    points_(points),
    pointSubsets_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
triSurfPoints::~triSurfPoints()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label triSurfPoints::addPointSubset(const word& subsetName)
{
    label id = pointSubsetIndex(subsetName);
    if( id >= 0 )
    {
        Warning << "Point subset " << subsetName << " already exists!" << endl;
        return id;
    }

    id = 0;
    forAllConstIter(Map<meshSubset>, pointSubsets_, it)
        id = Foam::max(id, it.key()+1);

    pointSubsets_.insert
    (
        id,
        meshSubset(subsetName, meshSubset::POINTSUBSET)
    );

    return id;
}

void triSurfPoints::removePointSubset(const label subsetID)
{
    if( pointSubsets_.find(subsetID) == pointSubsets_.end() )
        return;

    pointSubsets_.erase(subsetID);
}

word triSurfPoints::pointSubsetName(const label subsetID) const
{
    Map<meshSubset>::const_iterator it = pointSubsets_.find(subsetID);
    if( it == pointSubsets_.end() )
    {
        Warning << "Subset " << subsetID << " is not a point subset" << endl;
        return word();
    }

    return it().name();
}

label triSurfPoints::pointSubsetIndex(const word& subsetName) const
{
    forAllConstIter(Map<meshSubset>, pointSubsets_, it)
    {
        if( it().name() == subsetName )
            return it.key();
    }

    return -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
