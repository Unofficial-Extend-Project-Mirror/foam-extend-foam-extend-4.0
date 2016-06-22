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

Class
    changeMap

Description
    Accumulate topology change statistics

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Entity index
inline label& changeMap::index()
{
    return index_;
}


inline label changeMap::index() const
{
    return index_;
}


// Coupled patch index
inline label& changeMap::patchIndex()
{
    return pIndex_;
}


inline label changeMap::patchIndex() const
{
    return pIndex_;
}


// Type
inline label& changeMap::type()
{
    return type_;
}


inline label changeMap::type() const
{
    return type_;
}


inline void changeMap::addPoint
(
    const label pIndex,
    const labelList& master
)
{
    addedPoints_.append(objectMap(pIndex, master));
}


inline void changeMap::addEdge
(
    const label eIndex,
    const labelList& master
)
{
    addedEdges_.append(objectMap(eIndex, master));
}


inline void changeMap::addFace
(
    const label fIndex,
    const labelList& master
)
{
    addedFaces_.append(objectMap(fIndex, master));
}


inline void changeMap::addCell
(
    const label cIndex,
    const labelList& master
)
{
    addedCells_.append(objectMap(cIndex, master));
}


inline const List<objectMap>&
changeMap::addedPointList() const
{
    return addedPoints_;
}


inline const List<objectMap>&
changeMap::addedEdgeList() const
{
    return addedEdges_;
}


inline const List<objectMap>&
changeMap::addedFaceList() const
{
    return addedFaces_;
}


inline const List<objectMap>&
changeMap::addedCellList() const
{
    return addedCells_;
}


inline void changeMap::removePoint(const label pIndex)
{
    removedPoints_.append(pIndex);
}


inline void changeMap::removeEdge(const label eIndex)
{
    removedEdges_.append(eIndex);
}


inline void changeMap::removeFace(const label fIndex)
{
    removedFaces_.append(fIndex);
}


inline void changeMap::removeCell(const label cIndex)
{
    removedCells_.append(cIndex);
}


inline const labelList&
changeMap::removedPointList() const
{
    return removedPoints_;
}


inline const labelList&
changeMap::removedEdgeList() const
{
    return removedEdges_;
}


inline const labelList&
changeMap::removedFaceList() const
{
    return removedFaces_;
}


inline const labelList&
changeMap::removedCellList() const
{
    return removedCells_;
}


// Clear existing lists
inline void changeMap::clear()
{
    // Clear base dictionary
    dictionary::clear();

    index_ = -1;
    pIndex_ = -1;

    type_ = -1;

    // Clear added entities
    addedPoints_.clear();
    addedEdges_.clear();
    addedFaces_.clear();
    addedCells_.clear();

    // Clear removed entities
    removedPoints_.clear();
    removedEdges_.clear();
    removedFaces_.clear();
    removedCells_.clear();
}


inline void changeMap::operator=(const changeMap& rhs)
{
    // Copy base dictionary
    dictionary::operator=(rhs);

    index_ = rhs.index_;
    pIndex_ = rhs.pIndex_;

    type_ = rhs.type_;

    // Copy maps
    addedPoints_.setSize(rhs.addedPoints_.size());

    forAll(addedPoints_, indexI)
    {
        addedPoints_[indexI].index() = rhs.addedPoints_[indexI].index();

        addedPoints_[indexI].masterObjects() =
        (
            rhs.addedPoints_[indexI].masterObjects()
        );
    }

    addedEdges_.setSize(rhs.addedEdges_.size());

    forAll(addedEdges_, indexI)
    {
        addedEdges_[indexI].index() = rhs.addedEdges_[indexI].index();

        addedEdges_[indexI].masterObjects() =
        (
            rhs.addedEdges_[indexI].masterObjects()
        );
    }

    addedFaces_.setSize(rhs.addedFaces_.size());

    forAll(addedFaces_, indexI)
    {
        addedFaces_[indexI].index() = rhs.addedFaces_[indexI].index();

        addedFaces_[indexI].masterObjects() =
        (
            rhs.addedFaces_[indexI].masterObjects()
        );
    }

    addedCells_.setSize(rhs.addedCells_.size());

    forAll(addedCells_, indexI)
    {
        addedCells_[indexI].index() = rhs.addedCells_[indexI].index();

        addedCells_[indexI].masterObjects() =
        (
            rhs.addedCells_[indexI].masterObjects()
        );
    }

    removedPoints_ = rhs.removedPoints_;
    removedEdges_ = rhs.removedEdges_;
    removedFaces_ = rhs.removedFaces_;
    removedCells_ = rhs.removedCells_;
}


inline Ostream& operator<<(Ostream& os, const changeMap& cm)
{
    // Write base dictionary
    os << static_cast<const dictionary&>(cm) << endl;

    os << " index: " << cm.index_ << nl
       << " patchIndex: " << cm.pIndex_ << nl
       << " type: " << cm.type_ << nl
       << " addedPoints: " << cm.addedPoints_ << nl
       << " addedEdges: " << cm.addedEdges_ << nl
       << " addedFaces: " << cm.addedFaces_ << nl
       << " addedCells: " << cm.addedCells_ << nl
       << " removedPoints: " << cm.removedPoints_ << nl
       << " removedEdges: " << cm.removedEdges_ << nl
       << " removedFaces: " << cm.removedFaces_ << nl
       << " removedCells: " << cm.removedCells_ << endl;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, changeMap&)");

    return os;
}


} // End namespace Foam

// ************************************************************************* //
