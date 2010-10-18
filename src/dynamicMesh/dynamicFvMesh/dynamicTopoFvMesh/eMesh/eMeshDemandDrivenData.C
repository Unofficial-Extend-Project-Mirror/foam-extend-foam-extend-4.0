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

void eMesh::calcEdgePoints() const
{
    if (debug)
    {
        Info<< "void eMesh::calcEdgePoints() const : "
            << "Calculating ordered edgePoints" << endl;
    }

    if (epPtr_)
    {
        FatalErrorIn
        (
            "void eMesh::calcEdgePoints() const"
        )   << "edgePoints already allocated."
            << abort(FatalError);
    }

    if (!efPtr_->size())
    {
        FatalErrorIn
        (
            "void eMesh::calcEdgePoints() const"
        )   << "edgeFaces doesn't exist."
            << abort(FatalError);
    }

    // Size the list
    epPtr_ = new labelListList(nEdges_);
    labelListList& edgePoints = *epPtr_;

    // EdgePoints are ordered such that points appear in
    // counter-clockwise direction around point[0] of the edge.
    // If a boundary point is present, that must begin the list.
    // NOTE: Will work only on tetrahedral meshes!

    bool found;
    label faceIndex = -1, cellIndex = -1;
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const labelListList& eFaces = this->edgeFaces();

    for (label eIndex=0; eIndex < nEdges_; eIndex++)
    {
        const edge& e = edges_[eIndex];
        const labelList& eFace = eFaces[eIndex];

        // Size the list
        edgePoints[eIndex].setSize(eFace.size(), -1);

        if (eIndex < nInternalEdges_)
        {
            // Pick the first face and start with that
            faceIndex = eFace[0];
        }
        else
        {
            // Need to find a properly oriented start-face
            forAll(eFace, faceI)
            {
                if
                (
                    (!mesh_.isInternalFace(eFace[faceI]))
                 && (edgeDirection(faces[eFace[faceI]], e) == 1)
                )
                {
                    faceIndex = eFace[faceI];
                    break;
                }
            }
        }

        // Shuffle vertices to appear in CCW order
        for (label j=0; j < eFace.size(); j++)
        {
            const face& f = faces[faceIndex];

            // Add the isolated point
            edgePoints[eIndex][j] = findIsolatedPoint(f, e);

            // Figure out how this edge is oriented.
            if (edgeDirection(f, e) == 1)
            {
                // Counter-clockwise. Pick the owner.
                cellIndex = owner[faceIndex];
            }
            else
            if (mesh_.isInternalFace(faceIndex))
            {
                // Clockwise. Pick the neighbour.
                cellIndex = neighbour[faceIndex];
            }
            else
            {
                // Looks like we've hit a boundary face. Break out.
                break;
            }

            const cell& cellToCheck = cells[cellIndex];

            found = false;

            // Assuming tet-cells,
            // Loop through edgeFaces and get the next face
            forAll(eFace, faceI)
            {
                if
                (
                    eFace[faceI] != faceIndex
                 && eFace[faceI] == cellToCheck[0]
                )
                {
                    faceIndex = cellToCheck[0];
                    found = true; break;
                }

                if
                (
                    eFace[faceI] != faceIndex
                 && eFace[faceI] == cellToCheck[1]
                )
                {
                    faceIndex = cellToCheck[1];
                    found = true; break;
                }

                if
                (
                    eFace[faceI] != faceIndex
                 && eFace[faceI] == cellToCheck[2]
                )
                {
                    faceIndex = cellToCheck[2];
                    found = true; break;
                }

                if
                (
                    eFace[faceI] != faceIndex
                 && eFace[faceI] == cellToCheck[3]
                )
                {
                    faceIndex = cellToCheck[3];
                    found = true; break;
                }
            }

#           ifdef FULLDEBUG
            if (!found)
            {
                // Something's terribly wrong
                FatalErrorIn
                (
                    "void eMesh::calcEdgePoints() const"
                )
                << " Failed to determine a vertex ring. " << nl
                << " edgeFaces connectivity is inconsistent. " << nl
                << "Edge: " << eIndex << ":: " << e << nl
                << "edgeFaces: " << eFace
                << abort(FatalError);
            }
#           endif
        }
    }
}


void eMesh::calcPointEdges() const
{
    if (debug)
    {
        Info<< "void eMesh::calcPointEdges() const : "
            << "Calculating PointEdges" << endl;
    }

    if (pePtr_)
    {
        FatalErrorIn
        (
            "void eMesh::calcPointEdges() const"
        )   << "pePtr_ already allocated"
            << abort(FatalError);
    }

    pePtr_ = new labelListList(mesh_.nPoints());
    invertManyToMany(mesh_.nPoints(), edges(), *pePtr_);
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
