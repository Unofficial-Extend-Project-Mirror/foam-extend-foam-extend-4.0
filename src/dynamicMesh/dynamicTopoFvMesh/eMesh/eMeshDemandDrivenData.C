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

Description

    Demand-driven edge mesh data

\*---------------------------------------------------------------------------*/

#include "eMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Calculate ordered edges
void eMesh::calcOrderedEdgeList()
{
    if (debug)
    {
        Info<< "void eMesh::calcOrderedEdges() const : "
            << "Calculating ordered edges" << endl;
    }

    if (edges_.size() || boundary_.size())
    {
        FatalErrorIn
        (
            "void eMesh::calcOrderedEdges() const"
        )   << "Ordered edges already allocated."
            << abort(FatalError);
    }

    // Set size first.
    nEdges_ = mesh_.nEdges();
    nInternalEdges_ = 0;
    boundary_.setSize(mesh_.boundaryMesh().size());

    // Allocate lists for re-ordering
    labelList edgePatch(nEdges_, -1);
    labelList edgePatchStarts(mesh_.boundaryMesh().size(), -1);
    labelList edgePatchSizes(mesh_.boundaryMesh().size(), 0);

    reverseEdgeMap_.setSize(nEdges_);

    // Obtain connectivity from primitive mesh
    const edgeList& edges = mesh_.edges();
    const labelListList& fEdges = mesh_.faceEdges();

    // Edge-patches are the same as faces
    for (label i = mesh_.nInternalFaces(); i < mesh_.nFaces(); i++)
    {
        const labelList& fEdge = fEdges[i];

        forAll(fEdge, edgeI)
        {
            edgePatch[fEdge[edgeI]] = mesh_.boundaryMesh().whichPatch(i);
        }
    }

    // Loop through edgePatch and renumber internal edges
    forAll(edgePatch, edgeI)
    {
        if (edgePatch[edgeI] == -1)
        {
            reverseEdgeMap_[edgeI] = nInternalEdges_++;
        }
        else
        {
            edgePatchSizes[edgePatch[edgeI]]++;
        }
    }

    // Calculate patch-starts
    label startCount = nInternalEdges_;

    forAll(edgePatchStarts, patchI)
    {
        edgePatchStarts[patchI] = startCount;
        startCount += edgePatchSizes[patchI];
    }

    // Now renumber boundary edges
    labelList patchCount(edgePatchStarts);

    forAll(edgePatch, edgeI)
    {
        if (edgePatch[edgeI] > -1)
        {
            reverseEdgeMap_[edgeI] = patchCount[edgePatch[edgeI]]++;
        }
    }

    // Renumber and fill in edges
    edges_.setSize(nEdges_, edge(-1,-1));

    forAll(edges, edgeI)
    {
        edges_[reverseEdgeMap_[edgeI]] = edges[edgeI];
    }

    // Now set the boundary, copy name. (type is default)
    forAll(boundary_, patchI)
    {
        boundary_.set
        (
            patchI,
            ePatch::New
            (
                ePatch::typeName_(),
                mesh_.boundaryMesh()[patchI].name(),
                edgePatchSizes[patchI],
                edgePatchStarts[patchI],
                patchI,
                boundary_
            )
        );
    }

    // Force calculation of other data...
    calcEdgeFaces();
    calcFaceEdges();
}


void eMesh::calcFaceEdges() const
{
    if (debug)
    {
        Info<< "void eMesh::calcFaceEdges() const : "
            << "Calculating FaceEdges" << endl;
    }

    if (fePtr_)
    {
        FatalErrorIn
        (
            "void eMesh::calcFaceEdges() const"
        )   << "fePtr_ already allocated"
            << abort(FatalError);
    }

    fePtr_ =
    (
        new labelListIOList
        (
            IOobject
            (
                "faceEdges",
                mesh_.facesInstance(),
                meshSubDir,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_.nFaces()
        )
    );

    labelListIOList& faceEdges = *fePtr_;

    // If the read was successful, return.
    if (faceEdges.headerOk())
    {
        return;
    }

    // If edgeFaces exists, use that.
    if (efPtr_)
    {
        invertManyToMany(mesh_.nFaces(), edgeFaces(), faceEdges);
    }
    else
    {
        FatalErrorIn
        (
            "void eMesh::calcFaceEdges() const"
        )   << "Cannot calculate faceEdges."
            << abort(FatalError);
    }
}

void eMesh::calcEdgeFaces() const
{
    if (debug)
    {
        Info<< "void eMesh::calcEdgeFaces() const : "
            << "Calculating EdgeFaces" << endl;
    }

    if (efPtr_)
    {
        FatalErrorIn
        (
            "void eMesh::calcEdgeFaces() const"
        )   << "efPtr_ already allocated."
            << abort(FatalError);
    }

    efPtr_ =
    (
        new labelListIOList
        (
            IOobject
            (
                "edgeFaces",
                mesh_.facesInstance(),
                meshSubDir,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            nEdges_
        )
    );

    labelListIOList& edgeFaces = *efPtr_;

    // If the read was successful, return.
    if (edgeFaces.headerOk())
    {
        return;
    }

    if (fePtr_)
    {
        // If faceEdges exists, use that.
        invertManyToMany(nEdges_, faceEdges(), edgeFaces);
    }
    else
    {
        // Obtain connectivity from primitive mesh.
        if (!reverseEdgeMap_.size())
        {
            FatalErrorIn
            (
                "void eMesh::calcEdgeFaces() const"
            )   << "reverseEdgeMap has not been allocated."
                << abort(FatalError);
        }

        const labelListList& eFaces = mesh_.edgeFaces();

        forAll(eFaces, edgeI)
        {
            edgeFaces[reverseEdgeMap_[edgeI]] = eFaces[edgeI];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
