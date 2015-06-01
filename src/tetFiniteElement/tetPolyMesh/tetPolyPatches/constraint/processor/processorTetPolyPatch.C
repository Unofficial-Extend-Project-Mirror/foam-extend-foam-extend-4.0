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

\*---------------------------------------------------------------------------*/

#include "processorTetPolyPatch.H"
#include "tetPolyMesh.H"
#include "globalTetPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class FaceList>
labelList processorTetPolyPatch::calcProcLocalEdgesIndices
(
    const PrimitivePatch<face, FaceList, const pointField&>& p
) const
{
    if (debug)
    {
        Info<< "labelList processorTetPolyPatch::"
            << "calcProcLocalEdgesIndices(const primitivePatch& p) const : "
            << "calculating local edge indices"
            << endl;
    }

    // Count number of edges in the patch

    // Get reference to the mesh
    const tetPolyMesh& mesh = boundaryMesh().mesh();

    // Get reference to edges of the primitive patch
    const edgeList& patchEdges = p.edges();

    // Get reference to faces of the primitive patch
    const faceList& patchFaces = p;

    // Get reference to mesh points of the primitive patch
    const labelList& patchMeshPoints = p.meshPoints();

    // Make a point rejection list to remove shared processor points

    // Get reference to shared processor points
    const labelList& sharedPoints =
        boundaryMesh().globalPatch().meshPoints();

    boolList acceptPoint(patchMeshPoints.size(), true);

    forAll (patchMeshPoints, pointI)
    {
        label curP = patchMeshPoints[pointI];

        forAll (sharedPoints, sharedI)
        {
            if (sharedPoints[sharedI] == curP)
            {
                acceptPoint[pointI] = false;
                break;
            }
        }
    }

    // edges of the polyPatch
    label nEdgesInPatch = patchEdges.size();

    // diagonal edges across faces
    forAll (patchFaces, faceI)
    {
        nEdgesInPatch += patchFaces[faceI].size();
    }

    labelList localEdgeInd(nEdgesInPatch, -1);

    label nEdges = 0;

    const lduAddressing& lduAddr = mesh.lduAddr();

    // First do the edges of the primitive patch
    forAll (patchEdges, edgeI)
    {
        if
        (
            acceptPoint[patchEdges[edgeI].start()]
         || acceptPoint[patchEdges[edgeI].end()]
        )
        {
            localEdgeInd[nEdges] =
                lduAddr.triIndex
                (
                    patchMeshPoints[patchEdges[edgeI].start()],
                    patchMeshPoints[patchEdges[edgeI].end()]
                );

            nEdges++;
        }
    }

    // Now do the face-internal edges
    // All face edges are local to the plane

    // Calculate the offset to the first face centre
    label offset = mesh.faceOffset() + procPolyPatch_.start();

    forAll (patchFaces, faceI)
    {
        const face& curFace = patchFaces[faceI];

        forAll (curFace, pointI)
        {
            localEdgeInd[nEdges] =
                lduAddr.triIndex(curFace[pointI], offset + faceI);

            nEdges++;
        }
    }

    // Reset the size of the list
    localEdgeInd.setSize(nEdges);

#   ifdef DEBUGtetFemMatrix
    if (min(localEdgeInd) < 0)
    {
        FatalErrorIn
        (
            "void processorTetPolyPatch::"
            "calcProcLocalEdgesIndices(const primitivePatch& p) const"
        )   << "Problem in local edge addressing"
            << abort(FatalError);
    }
#   endif

    if (debug)
    {
        Info<< "labelList processorTetPolyPatch::"
            << "calcProcLocalEdgesIndices(const primitivePatch& p ) const : "
            << "finished calculating local edge indices"
            << endl;
    }

    return localEdgeInd;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
