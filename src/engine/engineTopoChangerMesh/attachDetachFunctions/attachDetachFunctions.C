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

#include "attachDetachFunctions.H"
#include "HashSet.H"

using namespace Foam;

// Find edge between points v0 and v1.
label Foam::findEdge
(
    const primitiveMesh& mesh,
    const label v0,
    const label v1
)
{
    const labelList& pEdges = mesh.pointEdges()[v0];

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        const edge& e = mesh.edges()[edgeI];

        if (e.otherVertex(v0) == v1)
        {
            return edgeI;
        }
    }

    FatalErrorIn
    (
        "findEdge(const primitiveMesh&, const label, const label)"
    )   << "Cannot find edge between mesh points " << v0 << " and " << v1
        << abort(FatalError);

    return -1;
}

// Checks whether patch present
void checkPatch
(
    const polyBoundaryMesh& bMesh,
    const word& name
)
{
    label patchI = bMesh.findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Cannot find patch " << name << endl
            << "It should be present but of zero size" << endl
            << "Valid patches are " << bMesh.names()
            << exit(FatalError);
    }

    if (bMesh[patchI].size() != 0)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but not of zero size"
            << exit(FatalError);
    }
}

void changePatchID
(
    const polyMesh& mesh,
    const label faceID,
    const label patchID,
    directTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}
