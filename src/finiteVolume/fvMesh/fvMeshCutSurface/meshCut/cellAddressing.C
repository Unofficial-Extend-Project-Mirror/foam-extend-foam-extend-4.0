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

#include "cellAddressing.H"
#include "labelHashSet.H"
#include "cellModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::labelList Foam::cellAddressing::makeIdentity(const label size)
{
    labelList result(size);
    forAll(result, resultI)
    {
        result[resultI] = resultI;
    }
    return result;
}


Foam::label Foam::cellAddressing::findEdge(const edge& e)
{
    forAll(edges_, edgeI)
    {
        if (edges_[edgeI] == e)
        {
            return edgeI;
        }
    }
    FatalErrorIn
    (
        "cellAddressing::findEdge(const edge&)"
    )   << "Problem: cannot find edge " << e << " in edges " << edges_
        << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellAddressing::cellAddressing(const Foam::cellModel& model)
:
    modelFaces_(model.modelFaces()),
    edges_(model.edges(makeIdentity(model.nPoints()))),
    edgeFaces_(model.nEdges()),
    faceEdges_(model.nFaces()),
    pointEdges_(model.nPoints())
{
    List<labelHashSet > dynEdgeFaces(edgeFaces_.size());

    forAll(faceEdges_, faceI)
    {
        DynamicList<label> fEdges;

        const face& f = modelFaces_[faceI];

        forAll(f, fp)
        {
            label edgeI = findEdge(edge(f[fp], f[(fp+1)%f.size()]));

            // Append to faceEdges
            fEdges.append(edgeI);

            // Append to edgefaces
            labelHashSet& eFaces = dynEdgeFaces[edgeI];

            if (!eFaces.found(faceI))
            {
                eFaces.insert(faceI);
            }
        }
        fEdges.shrink();

        // Copy into faceEdges_
        faceEdges_[faceI].transfer(fEdges.shrink());
    }

    // Convert dynEdgeFaces into edgeFaces_

    forAll(dynEdgeFaces, edgeI)
    {
        edgeFaces_[edgeI] = dynEdgeFaces[edgeI].toc();
    }


    // Convert edges_ into pointEdges
    List<DynamicList<label> > dynPointEdges(model.nPoints());

    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        dynPointEdges[e.start()].append(edgeI);
        dynPointEdges[e.end()].append(edgeI);
    }

    forAll(dynPointEdges, pointI)
    {
        pointEdges_[pointI].transfer(dynPointEdges[pointI].shrink());
    }
}


// ************************************************************************* //
