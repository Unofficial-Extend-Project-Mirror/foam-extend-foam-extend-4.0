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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfPoints.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label triSurfPoints::nPoints() const
{
    return points_.size();
}

inline const pointField& triSurfPoints::points() const
{
    return points_;
}

inline void triSurfPoints::appendVertex(const point& p)
{
    const label s = points_.size();
    points_.setSize(s+1);
    points_[s] = p;
}

inline void triSurfPoints::addPointToSubset
(
    const label setI,
    const label pointI
)
{
    Map<meshSubset>::iterator it = pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it().addElement(pointI);
}

inline void triSurfPoints::removePointFromSubset
(
    const label setI,
    const label pointI
)
{
    Map<meshSubset>::iterator it = pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it().removeElement(pointI);
}

inline void triSurfPoints::pointInSubsets
(
    const label pointI,
    DynList<label>& pointSubsets
) const
{
    pointSubsets.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        pointSubsets_,
        it
    )
    {
        if( it().contains(pointI) )
            pointSubsets.append(it.key());
    }
}

inline void triSurfPoints::pointSubsetIndices(DynList<label>& indices) const
{
    indices.clear();

    forAllConstIter
    (
        Map<meshSubset>,
        pointSubsets_,
        it
    )
        indices.append(it.key());
}

template<class ListType>
inline void triSurfPoints::pointsInSubset
(
    const label setI,
    ListType& pointLabels
) const
{
    pointLabels.clear();

    Map<meshSubset>::const_iterator it =
        pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it().containedElements(pointLabels);
}

template<class ListType>
inline void triSurfPoints::updatePointSubsets(const ListType& newNodeLabels)
{
    for
    (
        Map<meshSubset>::iterator it=pointSubsets_.begin();
        it!=pointSubsets_.end();
        ++it
    )
        it().updateSubset(newNodeLabels);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
