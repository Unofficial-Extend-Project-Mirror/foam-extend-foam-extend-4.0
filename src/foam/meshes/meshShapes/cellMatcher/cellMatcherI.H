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

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "cellModeller.H"
#include "cellModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::Map<Foam::label>& Foam::cellMatcher::localPoint() const
{
    return localPoint_;
}


inline const Foam::faceList& Foam::cellMatcher::localFaces() const
{
    return localFaces_;
}


inline const Foam::labelList& Foam::cellMatcher::faceSize() const
{
    return faceSize_;
}


inline const Foam::labelList& Foam::cellMatcher::pointMap() const
{
    return pointMap_;
}


inline const Foam::labelList& Foam::cellMatcher::faceMap() const
{
    return faceMap_;
}


inline const Foam::labelList& Foam::cellMatcher::edgeFaces() const
{
    return edgeFaces_;
}


inline const Foam::labelListList& Foam::cellMatcher::pointFaceIndex() const
{
    return pointFaceIndex_;
}


inline const Foam::labelList& Foam::cellMatcher::vertLabels() const
{
    return vertLabels_;
}


inline const Foam::labelList& Foam::cellMatcher::faceLabels() const
{
    return faceLabels_;
}


inline const Foam::cellModel& Foam::cellMatcher::model() const
{
    if (cellModelPtr_ == NULL)
    {
        cellModelPtr_ = cellModeller::lookup(cellModelName_);
    }
    return *cellModelPtr_;
}



// Key into edgeFaces_. key and key+1 are the entries for edge going from
// v0 to v1
inline Foam::label Foam::cellMatcher::edgeKey
(
    const label numVert,
    const label v0,
    const label v1
)
{
    return 2*(v0*numVert + v1);
}


// Walk along face consistent with face orientation
inline Foam::label Foam::cellMatcher::nextVert
(
    const label localVertI,
    const label size,
    const bool rightHand
)
{
    if (rightHand)
    {
        // face oriented acc. to righthand rule
        return (localVertI + 1) % size;
    }
    else
    {
        // face oriented acc. to lefthand rule
        return (size + localVertI - 1) % size;
    }
}


// ************************************************************************* //
