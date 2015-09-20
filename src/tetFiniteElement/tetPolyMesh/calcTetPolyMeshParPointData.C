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

Description

\*---------------------------------------------------------------------------*/

#include "tetPolyMesh.H"
#include "processorPolyPatch.H"
#include "demandDrivenData.H"

using namespace Foam;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetPolyMesh::calcParPointData() const
{
    if (!Pstream::parRun())
    {
        parPointsPtr_ = new labelList(0);
        parEdgesPtr_ = new edgeList(0);

        return;
    }
    else
    {
        // Make a global point markup list and find out which vertices
        // belong to more than one patch.

        // Create a local check list for all vertices

        labelList multiProcVertices(mesh_.nPoints(), 0);

        forAll (mesh_.boundaryMesh(), patchI)
        {
            if
            (
                isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI])
            )
            {
                const labelList& p = mesh_.boundaryMesh()[patchI].meshPoints();

                forAll (p, pointI)
                {
                    multiProcVertices[p[pointI]]++;
                }
            }
        }

        SLList<label> parPoints;

        // Find all vertices that have been addressed more than once
        forAll (multiProcVertices, pointI)
        {
            if (multiProcVertices[pointI] > 1)
            {
                parPoints.append(pointI);
            }
        }

        // Re-pack the list
        parPointsPtr_ = new labelList(parPoints);

        // At the same time calculate all edges located between two
        // parallel points.  The list will have duplicates, which need
        // to be eliminated

        SLList<edge> parEdges;

        forAll (mesh_.boundaryMesh(), patchI)
        {
            if
            (
                isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI])
            )
            {
                const labelList& p = mesh_.boundaryMesh()[patchI].meshPoints();

                const edgeList& e = mesh_.boundaryMesh()[patchI].edges();

                forAll (e, eI)
                {
                    if
                    (
                        multiProcVertices[p[e[eI].start()]] > 1
                     && multiProcVertices[p[e[eI].end()]] > 1
                    )
                    {
                        // Found a possible edge
                        edge newEdge = edge(p[e[eI].start()], p[e[eI].end()]);

                        // Check if the edge is already on the list
                        bool found = false;

                        for
                        (
                            SLList<edge>::iterator parEIter =
                                parEdges.begin();
                            parEIter != parEdges.end();
                            ++parEIter
                        )
                        {
                            if (parEIter() == newEdge)
                            {
                                found = true;
                                break;
                            }
                        }

                        if (!found)
                        {
                            parEdges.append(newEdge);
                        }
                    }
                }
            }
        }

        // Re-pack the list
        parEdgesPtr_ = new edgeList(parEdges);
    }
}


void Foam::tetPolyMesh::clearOutParPointData() const
{
    deleteDemandDrivenData(parPointsPtr_);
    deleteDemandDrivenData(parEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& Foam::tetPolyMesh::parallelPoints() const
{
    if (!parPointsPtr_)
    {
        calcParPointData();
    }

    return *parPointsPtr_;
}


const edgeList& Foam::tetPolyMesh::parallelEdges() const
{
    if (!parEdgesPtr_)
    {
        calcParPointData();
    }

    return *parEdgesPtr_;
}


// ************************************************************************* //
