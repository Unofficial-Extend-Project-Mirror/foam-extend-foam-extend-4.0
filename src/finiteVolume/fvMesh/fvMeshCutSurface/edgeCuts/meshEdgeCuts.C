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

Description

\*---------------------------------------------------------------------------*/

#include "meshEdgeCuts.H"

// * * * * * * * * * * * * * Private Static Functions  * * * * * * * * * * * //

void Foam::meshEdgeCuts::mark(const label elem, labelHashSet& markedElems)
{
    labelHashSet::const_iterator iter = markedElems.find(elem);

    if (iter == markedElems.end())
    {
        markedElems.insert(elem);
    }
}


void Foam::meshEdgeCuts::mark
(
    const labelList& elems,
    labelHashSet& markedElems
)
{
    forAll(elems, elemI)
    {
        mark(elems[elemI], markedElems);
    }
}


bool Foam::meshEdgeCuts::crosses
(
    const scalar isoVal,
    const scalar val0,
    const scalar val1,
    scalar& weight
)
{
    if
    ( 
        ((val0 <= isoVal) && (val1 >= isoVal))
     || ((val1 <= isoVal) && (val0 >= isoVal))
    )
    {
        if (val0 == val1)
        {
            return false;
        }
        else
        {
            weight = (isoVal - val0)/((val1 - val0) + VSMALL);

            return true;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshEdgeCuts::meshEdgeCuts
(
    const primitiveMesh& mesh,
    const labelList& cells,
    const labelList& meshVerts,
    const labelList& meshEdges,
    const scalarField& meshEdgeWeights
)
:
    mesh_(mesh),
    cells_(cells),
    meshVerts_(meshVerts),
    meshEdges_(meshEdges),
    meshEdgeWeights_(meshEdgeWeights)
{}


// ************************************************************************* //
