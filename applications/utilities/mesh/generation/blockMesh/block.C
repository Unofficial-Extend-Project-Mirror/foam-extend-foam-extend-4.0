/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "block.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


label block::vtxLabel(label a, label b, label c)
{
    return (a + b*(blockDef_.n().x() + 1)
            + c*(blockDef_.n().x() + 1)*(blockDef_.n().y() + 1));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from description
block::block(const blockDescriptor& definition)
:
    blockDef_(definition),
    vertices_
    (
        ((blockDef_.n().x() + 1)*(blockDef_.n().y() + 1)*(blockDef_.n().z() + 1))
    ),
    cells_
    (
        (blockDef_.n().x()*blockDef_.n().y()*blockDef_.n().z())
    ),
    boundaryPatches_(6)
{
    // create points
    blockPoints();

    // generate internal cells
    blockCells();

    // generate boundary patches
    blockBoundary();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const blockDescriptor& block::blockDef() const
{
    return blockDef_;
}

const pointField& block::points() const
{
    return vertices_;
}

const labelListList& block::cells() const
{
    return cells_;
}

const labelListListList& block::boundaryPatches() const
{
    return boundaryPatches_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const block& b)
{
    os << b.vertices_ << nl
       << b.cells_ << nl
       << b.boundaryPatches_ << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

