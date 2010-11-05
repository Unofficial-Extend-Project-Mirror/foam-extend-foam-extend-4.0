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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Type
inline label& changeMap::type()
{
    return type_;
}


inline label changeMap::type() const
{
    return type_;
}


// For sliver-type cells, opposite edges
// are identified for removal.
inline label& changeMap::firstEdge()
{
    return firstEdge_;
}


inline label& changeMap::secondEdge()
{
    return secondEdge_;
}


// For cap-type cells, the face requiring splitting
// is identified for removal.
inline label& changeMap::apexPoint()
{
    return apexPoint_;
}


inline label& changeMap::opposingFace()
{
    return opposingFace_;
}


// Added entities
inline void changeMap::addPoint
(
    const label pIndex,
    const label master
)
{
    label curSize = addedPoints_.size();

    addedPoints_.setSize(curSize + 1);

    addedPoints_[curSize][0] = pIndex;
    addedPoints_[curSize][1] = master;
}


inline void changeMap::addEdge
(
    const label eIndex,
    const label master
)
{
    label curSize = addedEdges_.size();

    addedEdges_.setSize(curSize + 1);

    addedEdges_[curSize][0] = eIndex;
    addedEdges_[curSize][1] = master;
}


inline void changeMap::addFace
(
    const label fIndex,
    const label master
)
{
    label curSize = addedFaces_.size();

    addedFaces_.setSize(curSize + 1);

    addedFaces_[curSize][0] = fIndex;
    addedFaces_[curSize][1] = master;
}


inline void changeMap::addCell
(
    const label cIndex,
    const label master
)
{
    label curSize = addedCells_.size();

    addedCells_.setSize(curSize + 1);

    addedCells_[curSize][0] = cIndex;
    addedCells_[curSize][1] = master;
}


// Return an added point
inline const List<FixedList<label,2> >&
changeMap::addedPointList() const
{
    return addedPoints_;
}


// Return the list of added entities
inline const List<FixedList<label,2> >&
changeMap::addedEdgeList() const
{
    return addedEdges_;
}


inline const List<FixedList<label,2> >&
changeMap::addedFaceList() const
{
    return addedFaces_;
}


inline const List<FixedList<label,2> >&
changeMap::addedCellList() const
{
    return addedCells_;
}


inline void changeMap::operator=(const changeMap& rhs)
{
    type_ = rhs.type_;

    firstEdge_    = rhs.firstEdge_;
    secondEdge_   = rhs.secondEdge_;
    apexPoint_    = rhs.apexPoint_;
    opposingFace_ = rhs.opposingFace_;

    addedPoints_  = rhs.addedPoints_;
    addedEdges_   = rhs.addedEdges_;
    addedFaces_   = rhs.addedFaces_;
    addedCells_   = rhs.addedCells_;
}

} // End namespace Foam

// ************************************************************************* //
