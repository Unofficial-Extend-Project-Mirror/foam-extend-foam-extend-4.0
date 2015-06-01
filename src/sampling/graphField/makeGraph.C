/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Global
    makeGraph

Description
    Write a graph file for a field given the data point locations field,
    the field of interest and the name of the field to be used for the
    graph file name.

\*---------------------------------------------------------------------------*/

#include "makeGraph.H"
#include "volFields.H"
#include "fvMesh.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void makeGraph
(
    const scalarField& x,
    const volScalarField& vsf,
    const word& graphFormat
)
{
    makeGraph(x, vsf, vsf.name(), graphFormat);
}


void makeGraph
(
    const scalarField& x,
    const volScalarField& vsf,
    const word& name,
    const word& graphFormat
)
{
    makeGraph(x, vsf.internalField(), name, vsf.path(), graphFormat);
}


void makeGraph
(
    const scalarField& x,
    const scalarField& sf,
    const word& name,
    const fileName& path,
    const word& graphFormat
)
{
    graph
    (
        name,
        "x",
        name,
        x,
        sf
    ).write(path/name, graphFormat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
