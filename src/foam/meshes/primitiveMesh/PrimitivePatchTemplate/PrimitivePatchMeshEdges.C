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

#include "PrimitivePatchTemplate.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::labelList
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
meshEdges
(
    const edgeList& allEdges,
    const labelListList& cellEdges,
    const labelList& faceCells
) const
{
    if (debug)
    {
        Info<< "labelList PrimitivePatch<Face, FaceList, PointField, PointType>"
            << "::meshEdges() : "
            << "calculating labels of patch edges in mesh edge list"
            << endl;
    }

    // get reference to the list of edges on the patch
    const edgeList& PatchEdges = edges();

    const labelListList& EdgeFaces = edgeFaces();

    // create the storage
    labelList meshEdges(PatchEdges.size());

    register bool found = false;

    // get reference to the points on the patch
    const labelList& pp = meshPoints();

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll (PatchEdges, edgeI)
    {
        const edge curEdge
            (pp[PatchEdges[edgeI].start()], pp[PatchEdges[edgeI].end()]);

        found = false;

        // get the patch faces sharing the edge
        const labelList& curFaces = EdgeFaces[edgeI];

        forAll (curFaces, faceI)
        {
            // get the cell next to the face
            label curCell = faceCells[curFaces[faceI]];

            // get reference to edges on the cell
            const labelList& ce = cellEdges[curCell];

            forAll (ce, cellEdgeI)
            {
                if (allEdges[ce[cellEdgeI]] == curEdge)
                {
                    found = true;

                    meshEdges[edgeI] = ce[cellEdgeI];

                    break;
                }
            }

            if (found) break;
        }
    }

    return meshEdges;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::labelList
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
meshEdges
(
    const edgeList& allEdges,
    const labelListList& pointEdges
) const
{
    if (debug)
    {
        Info<< "labelList PrimitivePatch<Face, FaceList, PointField, PointType>"
            << "::meshEdges() : "
            << "calculating labels of patch edges in mesh edge list"
            << endl;
    }

    // get reference to the list of edges on the patch
    const edgeList& PatchEdges = edges();

    // create the storage
    labelList meshEdges(PatchEdges.size());

    // get reference to the points on the patch
    const labelList& pp = meshPoints();

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll (PatchEdges, edgeI)
    {
        const label globalPointI = pp[PatchEdges[edgeI].start()];
        const edge curEdge(globalPointI, pp[PatchEdges[edgeI].end()]);

        const labelList& pe = pointEdges[globalPointI];

        forAll (pe, i)
        {
            if (allEdges[pe[i]] == curEdge)
            {
                meshEdges[edgeI] = pe[i];
                break;
            }
        }
    }

    return meshEdges;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::label
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
whichEdge
(
    const edge& e
) const
{
    // Get pointEdges from the starting point and search all the candidates
    const edgeList& Edges = edges();

    if (e.start() > -1 && e.start() < nPoints())
    {
        const labelList& pe = pointEdges()[e.start()];

        forAll (pe, peI)
        {
            if (e == Edges[pe[peI]])
            {
                return pe[peI];
            }
        }
    }

    // Edge not found.  Return -1
    return -1;
}


// ************************************************************************* //
