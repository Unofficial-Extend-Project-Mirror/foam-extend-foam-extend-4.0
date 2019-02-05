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
    An graph of faces which supports automated output.

\*---------------------------------------------------------------------------*/

#include "cellIOGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellIOGraph::cellIOGraph(const IOobject& io)
:
    regIOobject(io),
    VRWGraph()
{}

cellIOGraph::cellIOGraph
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io),
    VRWGraph(size)
{}

cellIOGraph::cellIOGraph
(
    const IOobject& io,
    const VRWGraph& g
)
:
    regIOobject(io),
    VRWGraph(g)
{}

void cellIOGraph::operator=(const cellIOGraph& rhs)
{
    VRWGraph::operator=(rhs);
}

void cellIOGraph::operator=(const VRWGraph& rhs)
{
    VRWGraph::operator=(rhs);
}

bool cellIOGraph::writeData(Ostream& os) const
{
    return (os << *this).good();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameWithName(cellIOGraph, "cellList");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
