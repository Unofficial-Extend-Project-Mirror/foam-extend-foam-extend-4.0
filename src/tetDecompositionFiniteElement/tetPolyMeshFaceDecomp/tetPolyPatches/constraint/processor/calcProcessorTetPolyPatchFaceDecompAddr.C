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
    Calculate cut edge processor addressing, needed for vector-matrix
    multiply on processor boundaries.

\*---------------------------------------------------------------------------*/

#include "processorTetPolyPatchFaceDecomp.H"
#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "boolList.H"
#include "globalTetPolyPatchFaceDecomp.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void processorTetPolyPatchFaceDecomp::calcCutEdgeIndices() const
{
    if (debug)
    {
        Info<< "void processorTetPolyPatchFaceDecomp::"
            << "calcCutEdgeIndices() const : "
            << endl << "calculating cut edge indices"
            << endl;
    }

    if (cutEdgeIndicesPtr_)
    {
        FatalErrorIn
        (
            "void processorTetPolyPatchFaceDecomp::calcCutEdgesIndices() const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }

    // Make a list over all edges in the mesh.  Mark the ones that are local
    // to the patch and then collect the rest

    boolList isLocal(boundaryMesh().mesh().nEdges(), false);

    // get reference to local edge indices
    const labelList& localEdges = localEdgeIndices();

    forAll (localEdges, edgeI)
    {
        isLocal[localEdges[edgeI]] = true;
    }

    // Get reference to parallel shared edges
    const labelList& sharedEdges =
        boundaryMesh().globalPatch().localEdgeIndices();

    forAll (sharedEdges, sharedI)
    {
        isLocal[sharedEdges[sharedI]] = true;
    }

    const labelList& sharedCutEdges =
        boundaryMesh().globalPatch().cutEdgeIndices();

    forAll (sharedCutEdges, sharedCutI)
    {
        isLocal[sharedCutEdges[sharedCutI]] = true;
    }

    // Count the maximum number of edges coming from the patch
    label maxEdgesOnPatch = 0;

    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    const labelList& mp = meshPoints();

    forAll (mp, pointI)
    {
        maxEdgesOnPatch += mesh.nEdgesForPoint(mp[pointI]);
    }

    // Allocate the array
    cutEdgeIndicesPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& cutEdgeInd = *cutEdgeIndicesPtr_;

    label nCutEdgeInd = 0;

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as local;
    // if not, add it to the list of cut edges.
    forAll (mp, pointI)
    {
        labelList curEdges = mesh.edgesForPoint(mp[pointI]);

        forAll (curEdges, edgeI)
        {
            if (!isLocal[curEdges[edgeI]])
            {
                cutEdgeInd[nCutEdgeInd] = curEdges[edgeI];
                nCutEdgeInd++;
            }
        }
    }

    // Reset the size of the edge list
    cutEdgeInd.setSize(nCutEdgeInd);

    if (debug)
    {
        Info<< "void processorTetPolyPatchFaceDecomp::"
            << "calcCutEdgeIndices() const : "
            << endl << "finished calculating cut edge indices"
            << endl;
    }
}


void processorTetPolyPatchFaceDecomp::calcCutEdgeAddressing() const
{
    if
    (
        cutEdgeOwnerIndicesPtr_
     || cutEdgeOwnerStartPtr_
     || cutEdgeNeighbourIndicesPtr_
     || cutEdgeNeighbourStartPtr_
     || doubleCutEdgeIndicesPtr_
     || doubleCutOwnerPtr_
     || doubleCutNeighbourPtr_
    )
    {
        FatalErrorIn
        (
            "void processorTetPolyPatchFaceDecomp::calcCutEdgeAddressing() "
            "const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }


    // Make a list over all edges in the mesh.  Mark the ones that are local
    // to the patch and then collect the rest
    // For doubly cut edges, mark up the local points

    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    // Get reference to local edge indices
    const labelList& localEdges = localEdgeIndices();

    // Get reference to mesh point addressing for the patch
    const labelList& mp = meshPoints();

    // Mark up the local edges

    boolList isLocal(mesh.nEdges(), false);

    forAll (localEdges, edgeI)
    {
        isLocal[localEdges[edgeI]] = true;
    }

    // Get reference to parallel shared edges
    const labelList& sharedEdges =
        boundaryMesh().globalPatch().localEdgeIndices();

    forAll (sharedEdges, sharedI)
    {
        isLocal[sharedEdges[sharedI]] = true;
    }

    const labelList& sharedCutEdges =
        boundaryMesh().globalPatch().cutEdgeIndices();

    forAll (sharedCutEdges, sharedCutI)
    {
        isLocal[sharedCutEdges[sharedCutI]] = true;
    }

    labelList localPointLabel(mesh.nPoints(), -1);

    forAll (mp, pointI)
    {
        localPointLabel[mp[pointI]] = pointI;
    }

    // Count the maximum number of edges coming from the patch
    label maxEdgesOnPatch = 0;

    forAll (mp, pointI)
    {
        maxEdgesOnPatch += mesh.nEdgesForPoint(mp[pointI]);
    }

    // Get reference to addressing
    const lduAddressing& ldu = mesh.lduAddr();

    // Allocate doubly cut arrays
    doubleCutEdgeIndicesPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& doubleCutEdges = *doubleCutEdgeIndicesPtr_;

    doubleCutOwnerPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& doubleCutOwn = *doubleCutOwnerPtr_;

    doubleCutNeighbourPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& doubleCutNei = *doubleCutNeighbourPtr_;

    label nDoubleCut = 0;

    const labelList& globalOwner = ldu.lowerAddr();
    const labelList& globalNeighbour = ldu.upperAddr();

    // Owner side
    //~~~~~~~~~~~

    // Allocate the array
    cutEdgeOwnerIndicesPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& own = *cutEdgeOwnerIndicesPtr_;
    label nOwn = 0;

    cutEdgeOwnerStartPtr_ = new labelList(meshPoints().size() + 1, -1);
    labelList& ownStart = *cutEdgeOwnerStartPtr_;

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as local;
    // if not, add it to the list of cut edges.
    forAll (mp, pointI)
    {
        ownStart[pointI] = nOwn;

        const label curPointID = mp[pointI];

        // get owner edges indices
        const label startFaceOwn = ldu.ownerStartAddr()[curPointID];
        const label endFaceOwn = ldu.ownerStartAddr()[curPointID + 1];

        for
        (
            label edgeLabel = startFaceOwn;
            edgeLabel < endFaceOwn;
            edgeLabel++
        )
        {
            if (!isLocal[edgeLabel])
            {
                isLocal[edgeLabel] = true;

                if (localPointLabel[globalNeighbour[edgeLabel]] == -1)
                {
                    // Singly cut edge
                    own[nOwn] = edgeLabel;
                    nOwn++;
                }
                else
                {
                    // Doubly cut edge
                    doubleCutEdges[nDoubleCut] = edgeLabel;
                    doubleCutOwn[nDoubleCut] = pointI;
                    doubleCutNei[nDoubleCut] =
                        localPointLabel[globalNeighbour[edgeLabel]];

                    nDoubleCut++;
                }
            }
        }
    }

    // reset the size of owner edges
    own.setSize(nOwn);

    // set the last start label by hand
    ownStart[meshPoints().size()] = nOwn;


    // Neighbour side
    //~~~~~~~~~~~~~~~

    // Allocate the array
    cutEdgeNeighbourIndicesPtr_ = new labelList(maxEdgesOnPatch, -1);
    labelList& nei = *cutEdgeNeighbourIndicesPtr_;
    label nNei = 0;

    cutEdgeNeighbourStartPtr_ = new labelList(meshPoints().size() + 1, -1);
    labelList& neiStart = *cutEdgeNeighbourStartPtr_;

    const unallocLabelList& losort = ldu.losortAddr();

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as local;
    // if not, add it to the list of cut edges.
    forAll (mp, pointI)
    {
        neiStart[pointI] = nNei;

        const label curPointID = mp[pointI];

        // get neighbour edges indices
        const label startFaceNei = ldu.losortStartAddr()[curPointID];
        const label endFaceNei = ldu.losortStartAddr()[curPointID + 1];

        for
        (
            label edgeLabel = startFaceNei;
            edgeLabel < endFaceNei;
            edgeLabel++
        )
        {
            if (!isLocal[losort[edgeLabel]])
            {
                isLocal[losort[edgeLabel]] = true;

                if (localPointLabel[globalOwner[losort[edgeLabel]]] == -1)
                {
                    // Singly cut edge
                    nei[nNei] = losort[edgeLabel];
                    nNei++;
                }
                else
                {
                    // Doubly cut edge
                    doubleCutEdges[nDoubleCut] = losort[edgeLabel];
                    doubleCutOwn[nDoubleCut] =
                        localPointLabel[globalOwner[losort[edgeLabel]]];
                    doubleCutNei[nDoubleCut] = pointI;

                    nDoubleCut++;
                }
            }
        }
    }

    // reset the size of neighbour edges
    nei.setSize(nNei);

    // set the last start label by hand
    neiStart[meshPoints().size()] = nNei;

    // Reset the size of double cut edge data
    doubleCutEdges.setSize(nDoubleCut);
    doubleCutOwn.setSize(nDoubleCut);
    doubleCutNei.setSize(nDoubleCut);
}


void processorTetPolyPatchFaceDecomp::calcOwnNeiDoubleMask() const
{
    if (ownNeiDoubleMaskPtr_)
    {
        FatalErrorIn
        (
            "void processorTetPolyPatchFaceDecomp::calcOwnNeiDoubleMask() "
            "const"
        )   << "ownNeiDoubleMaskPtr_ already allocated"
            << abort(FatalError);
    }

    // Algorithm:
    // Make a point lookup field which marks all points belonging to
    // Other processor patches.
    // Go through all the cut edges and cheeck the points of the edge.
    // If the point belongs to another parallel patch, set the mask to zero

    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    boolList otherProcPoint(mesh.nPoints(), false);

    forAll (boundaryMesh(), patchI)
    {
        if
        (
            isA<processorTetPolyPatchFaceDecomp>(boundaryMesh()[patchI])
        )
        {
            // Get the local points and mark up the list
            const labelList& mp = boundaryMesh()[patchI].meshPoints();

            forAll (mp, i)
            {
                otherProcPoint[mp[i]] = true;
            }
        }
    }

    // Go through all cut edges and set the mask

    // Get matrix addressing
    const labelList& mp = meshPoints();

    const unallocLabelList& L = mesh.lduAddr().lowerAddr();
    const unallocLabelList& U = mesh.lduAddr().upperAddr();

    const labelList& cutOwn = cutEdgeOwnerIndices();
    const labelList& cutNei = cutEdgeNeighbourIndices();
    const labelList& doubleCut = doubleCutEdgeIndices();

    ownNeiDoubleMaskPtr_ =
        new scalarField
        (
            cutOwn.size() + cutNei.size() + doubleCut.size(),
            1.0
        );
    scalarField& mask = *ownNeiDoubleMaskPtr_;

    label coeffI = 0;

    // Owner side
    // ~~~~~~~~~~
    {
        const labelList& cutOwnStart = cutEdgeOwnerStart();

        forAll (mp, pointI)
        {
            label ownIndex = cutOwnStart[pointI];
            label endOwn = cutOwnStart[pointI + 1];

            for (; ownIndex < endOwn; ownIndex++)
            {
                if(otherProcPoint[U[cutOwn[ownIndex]]])
                {
                    mask[coeffI] = 0;
                }
                coeffI++;
            }
        }
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    {
        const labelList& cutNeiStart = cutEdgeNeighbourStart();

        forAll (mp, pointI)
        {
            label neiIndex = cutNeiStart[pointI];
            label endNei = cutNeiStart[pointI + 1];

            for (; neiIndex < endNei; neiIndex++)
            {
                if (otherProcPoint[L[cutNei[neiIndex]]])
                {
                    mask[coeffI] = 0;
                }

                coeffI++;
            }
        }
    }

    // Doubly cut coefficients
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Not needed.  HJ, 18/Apr/2002
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& processorTetPolyPatchFaceDecomp::localEdgeIndices() const
{
    if (!localEdgeIndicesPtr_)
    {
        if (isMaster())
        {
            localEdgeIndicesPtr_ =
                new labelList(calcProcLocalEdgesIndices(procPolyPatch_));
        }
        else
        {
            // Slave side.  Create the reversed patch and pick up its points
            // so that the order is correct
            const polyPatch& pp = patch();

            faceList masterFaces(pp.size());

            forAll (pp, faceI)
            {
                masterFaces[faceI] = pp[faceI].reverseFace();
            }

            localEdgeIndicesPtr_ =
                new labelList
                (
                    calcProcLocalEdgesIndices
                    (
                        primitiveFacePatch
                        (
                            masterFaces,
                            pp.points()
                        )
                    )
                );
        }
    }

    return *localEdgeIndicesPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::cutEdgeIndices() const
{
    if (!cutEdgeIndicesPtr_)
    {
        calcCutEdgeIndices();
    }

    return *cutEdgeIndicesPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::cutEdgeOwnerIndices() const
{
    if (!cutEdgeOwnerIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerIndicesPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::cutEdgeOwnerStart() const
{
    if (!cutEdgeOwnerStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerStartPtr_;
}


const labelList&
processorTetPolyPatchFaceDecomp::cutEdgeNeighbourIndices() const
{
    if (!cutEdgeNeighbourIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourIndicesPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::cutEdgeNeighbourStart() const
{
    if (!cutEdgeNeighbourStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourStartPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::doubleCutEdgeIndices() const
{
    if (!doubleCutEdgeIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutEdgeIndicesPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::doubleCutOwner() const
{
    if (!doubleCutEdgeIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutOwnerPtr_;
}


const labelList& processorTetPolyPatchFaceDecomp::doubleCutNeighbour() const
{
    if (!doubleCutEdgeIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutNeighbourPtr_;
}


const scalarField& processorTetPolyPatchFaceDecomp::ownNeiDoubleMask() const
{
    if (!ownNeiDoubleMaskPtr_)
    {
        calcOwnNeiDoubleMask();
    }

    return *ownNeiDoubleMaskPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
