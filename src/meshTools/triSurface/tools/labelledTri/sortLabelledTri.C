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

#include "sortLabelledTri.H"
#include "labelledTri.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

inline bool surfAndLabel::less::operator()
(
    const surfAndLabel& one,
    const surfAndLabel& two
) const
{
    const triSurface& surf = *one.surfPtr_;
    return surf[one.index_].region() < surf[two.index_].region();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sortLabelledTri::sortLabelledTri(const triSurface& surf)
:
    List<surfAndLabel>(surf.size(), surfAndLabel(surf, -1))
{

    // Set the face label
    forAll(surf, faceI)
    {
        operator[](faceI).index_ = faceI;
    }

    // Sort according to region number.
    sort(*this, surfAndLabel::less());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sortLabelledTri::indices(labelList& newIndices) const
{
    newIndices.setSize(size());

    forAll(newIndices, i)
    {
        newIndices[i] = operator[](i).index_;
    }
}


labelList sortLabelledTri::indices() const
{
    labelList newIndices(size());
    indices(newIndices);
    return newIndices;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
