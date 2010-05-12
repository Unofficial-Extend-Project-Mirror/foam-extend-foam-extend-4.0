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
    Calculate cut edge global addressing, needed for vector-matrix
    multiply on global boundaries.

\*---------------------------------------------------------------------------*/

#include "globalTetPolyPatchFaceDecomp.H"
#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void globalTetPolyPatchFaceDecomp::calcLocalEdgesIndices() const
{
    if (debug)
    {
        Info<< "labelList globalTetPolyPatchFaceDecomp::"
            << "calcLocalEdgesIndices() const : "
            << "calculating local edge indices"
            << endl;
    }

    // Get reference to the mesh
    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    // Get reference to edges
    const edgeList& patchEdges = meshEdges();

    localEdgeIndicesPtr_ = new labelList(patchEdges.size(), -1);
    labelList& localEdgeInd = *localEdgeIndicesPtr_;

    const lduAddressing& lduAddr = mesh.lduAddr();

    forAll (patchEdges, edgeI)
    {
        localEdgeInd[edgeI] =
            lduAddr.triIndex
            (
                patchEdges[edgeI].start(),
                patchEdges[edgeI].end()
            );
    }

#   ifdef DEBUGtetFemMatrix
    if (localEdgeInd.size() > 0 && min(localEdgeInd) < 0)
    {
        FatalErrorIn
        (
            "void globalTetPolyPatchFaceDecomp::"
            "calcLocalEdgesIndices() const"
        )   << "Problem in local edge addressing"
            << abort(FatalError);
    }
#   endif

    if (debug)
    {
        Info<< "void globalTetPolyPatchFaceDecomp::"
            << "calcLocalEdgesIndices() const : "
            << "finished calculating local edge indices"
            << endl;
    }
}


void globalTetPolyPatchFaceDecomp::calcCutEdgeIndices() const
{
    if (debug)
    {
        Info<< "void globalTetPolyPatchFaceDecomp::"
            << "calcCutEdgeIndices() const : "
            << "calculating cut edge indices"
            << endl;
    }

    if (cutEdgeIndicesPtr_)
    {
        FatalErrorIn
        (
            "void globalTetPolyPatchFaceDecomp::"
            "calcCutEdgesIndices() const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }

    // Get reference to the mesh
    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    // Get reference to edges
    const edgeList& patchCutEdges = meshCutEdges();

    cutEdgeIndicesPtr_ = new labelList(patchCutEdges.size(), -1);
    labelList& cutEdgeInd = *cutEdgeIndicesPtr_;

    const lduAddressing& lduAddr = mesh.lduAddr();

    forAll (patchCutEdges, edgeI)
    {
        cutEdgeInd[edgeI] =
            lduAddr.triIndex
            (
                patchCutEdges[edgeI].start(),
                patchCutEdges[edgeI].end()
            );
    }

    if (debug)
    {
        Info<< "void globalTetPolyPatchFaceDecomp::"
            << "calcCutEdgeIndices() const : "
            << "finished calculating cut edge indices"
            << endl;
    }
}


void globalTetPolyPatchFaceDecomp::calcCutEdgeAddressing() const
{
    if
    (
        cutEdgeOwnerIndicesPtr_
     || cutEdgeOwnerStartPtr_
     || cutEdgeNeighbourIndicesPtr_
     || cutEdgeNeighbourStartPtr_
     || ownNeiDoubleMaskPtr_
    )
    {
        FatalErrorIn
        (
            "void globalTetPolyPatchFaceDecomp::"
            "calcCutEdgeAddressing() "
            "const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }


    // Make a list over all edges in the mesh.  Mark the ones that are local
    // to the patch and then collect the rest.
    // For doubly cut edges, mark up the local points

    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    // Get reference to mesh points
    const labelList& mp = meshPoints();

    // Get reference to local edge indices
    const labelList& cutEdges = cutEdgeIndices();

    // Get the cut edge mask.  It will be used to create the new mask
    const scalarField& cutMask = meshCutEdgeMask();

    // Mark up the cut edges

    labelList cutIndexLookup(mesh.nEdges(), -1);

    forAll (cutEdges, edgeI)
    {
        cutIndexLookup[cutEdges[edgeI]] = edgeI;
    }

    // Get reference to addressing
    const lduAddressing& ldu = mesh.lduAddr();

    // Owner side
    //~~~~~~~~~~~

    // Allocate the array
    cutEdgeOwnerIndicesPtr_ = new labelList(cutEdges.size(), -1);
    labelList& own = *cutEdgeOwnerIndicesPtr_;
    label nOwn = 0;

    cutEdgeOwnerStartPtr_ = new labelList(mp.size() + 1, -1);
    labelList& ownStart = *cutEdgeOwnerStartPtr_;

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as cut;
    // add it to the list of cut edges and grab the mask value that goes
    // with it.
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
            if (cutIndexLookup[edgeLabel] >= 0)
            {
                // Singly cut edge
                own[nOwn] = edgeLabel;
                nOwn++;
            }
        }
    }

    // reset the size of owner edges
    own.setSize(nOwn);

    // set the last start label by hand
    ownStart[mp.size()] = nOwn;

    // Neighbour side
    //~~~~~~~~~~~~~~~

    // Allocate the array
    cutEdgeNeighbourIndicesPtr_ = new labelList(cutEdges.size(), -1);
    labelList& nei = *cutEdgeNeighbourIndicesPtr_;
    label nNei = 0;

    cutEdgeNeighbourStartPtr_ = new labelList(mp.size() + 1, -1);
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
            if (cutIndexLookup[losort[edgeLabel]] >= 0)
            {
                nei[nNei] = losort[edgeLabel];
                nNei++;
            }
        }
    }

    // reset the size of neighbour edges
    nei.setSize(nNei);

    // set the last start label by hand
    neiStart[mp.size()] = nNei;

    // Allocate the reordered mask that corresponds to the owner-neighbour
    // coefficient ordering
    ownNeiDoubleMaskPtr_ = new scalarField(own.size() + nei.size());
    scalarField& ownNeiDoubleMask = *ownNeiDoubleMaskPtr_;

    label nMask = 0;

    forAll (own, i)
    {
        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[own[i]]];
        nMask++;
    }

    forAll (nei, i)
    {
        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[nei[i]]];
        nMask++;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& globalTetPolyPatchFaceDecomp::localEdgeIndices() const
{
    if (!localEdgeIndicesPtr_)
    {
        calcLocalEdgesIndices();
    }

    return *localEdgeIndicesPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::cutEdgeIndices() const
{
    if (!cutEdgeIndicesPtr_)
    {
        calcCutEdgeIndices();
    }

    return *cutEdgeIndicesPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::cutEdgeOwnerIndices() const
{
    if (!cutEdgeOwnerIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerIndicesPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::cutEdgeOwnerStart() const
{
    if (!cutEdgeOwnerStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerStartPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::cutEdgeNeighbourIndices() const
{
    if (!cutEdgeNeighbourIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourIndicesPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::cutEdgeNeighbourStart() const
{
    if (!cutEdgeNeighbourStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourStartPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::doubleCutEdgeIndices() const
{
    if (!doubleCutEdgeIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutEdgeIndicesPtr_;
}


const labelList& globalTetPolyPatchFaceDecomp::doubleCutOwner() const
{
    if (!doubleCutOwnerPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutOwnerPtr_;
}


const labelList&
globalTetPolyPatchFaceDecomp::doubleCutNeighbour() const
{
    if (!doubleCutNeighbourPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutNeighbourPtr_;
}


const scalarField&
globalTetPolyPatchFaceDecomp::ownNeiDoubleMask() const
{
    if (!ownNeiDoubleMaskPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *ownNeiDoubleMaskPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
