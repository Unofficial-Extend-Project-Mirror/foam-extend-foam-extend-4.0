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

#include "faMesh.H"
#include "mapPolyMesh.H"
#include "MapFaFields.H"
#include "faMeshMapper.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faMesh::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        Info<< "bool faMesh::updateMesh(const mapPolyMesh& mpm) const : "
            << "Updating mesh" << endl;
    }

    if (!mpm.morphing())
    {
        // No topo change
        return false;
    }

    // Create fa mesh mapper, using the old mesh
    const faMeshMapper mapper(*this, mpm);


    // Rebuild mesh

    // Cast away const for interface reasons.  HJ, 12/Aug/2011
    faMesh& m = const_cast<faMesh&>(*this);


    // Clear existing mesh data
    clearOut();

    // Set new labels
    m.faceLabels_ = mapper.areaMap().newFaceLabels();

    const indirectPrimitivePatch& bp = patch();

    // Collect patch data
    const label nTotalEdges = bp.nEdges();
    const label nInternalEdges = bp.nInternalEdges();
    const labelListList& edgeFaces = bp.edgeFaces();

    labelListList patchEdges(boundary_.size());

    // Special handling required for faces that have more than one edge
    // Each patch will be visited separately

    labelList edgeToPatch(nTotalEdges - nInternalEdges, -1);
    const labelList& newFaceLabelsMap = mapper.areaMap().newFaceLabelsMap();

    const labelListList& oldPatchEdgeFaces = mapper.oldPatchEdgeFaces();

    forAll (oldPatchEdgeFaces, patchI)
    {
        labelList& curPatchEdges = patchEdges[patchI];
        curPatchEdges.setSize(nTotalEdges - nInternalEdges);
        label nCurPatchEdges = 0;

        // Note: it is possible to pick up the old-to-new boundary patch
        // mapping, but currently this is not done.  HJ, 13/Aug/2011

        // Make a fast lookup
        labelHashSet oldFaceLookup(oldPatchEdgeFaces[patchI]);

        for  (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
        {
            if (edgeToPatch[edgeI - nInternalEdges] > -1)
            {
                // Edge already found; continue with the next one
                continue;
            }

            // Boundary edges will only have one face next to them
            const label oldFaceIndex = newFaceLabelsMap[edgeFaces[edgeI][0]];

            if (oldFaceIndex > -1)
            {
                // Old face exists.  See if it has got an edge in this patch
                if (oldFaceLookup.found(oldFaceIndex))
                {
                    // Face found, add it to the patch
                    curPatchEdges[nCurPatchEdges] = edgeI;
                    nCurPatchEdges++;

                    edgeToPatch[edgeI - nInternalEdges] = patchI;
                }
            }
        }

        // Collected all faces for the current patch
        curPatchEdges.setSize(nCurPatchEdges);
    }

    // Set new edges for all patches
    forAll (m.boundary_, patchI)
    {
        m.boundary_[patchI].resetEdges(patchEdges[patchI]);
    }

    m.setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    m.boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    m.boundary_.calcGeometry();


    // Map fields
    mapFields(mapper);

    // Map old areas
    mapOldAreas(mapper);

    // Update edge interpolation
    edgeInterpolation::movePoints();

    return true;
}


void Foam::faMesh::mapFields(const faMeshMapper& mapper) const
{
    // Map all the areaFields in the objectRegistry
    MapGeometricFields<scalar, faPatchField, faMeshMapper, areaMesh>(mapper);
    MapGeometricFields<vector, faPatchField, faMeshMapper, areaMesh>(mapper);
    MapGeometricFields<sphericalTensor, faPatchField, faMeshMapper, areaMesh>
        (mapper);
    MapGeometricFields<symmTensor, faPatchField, faMeshMapper, areaMesh>
        (mapper);
    MapGeometricFields<tensor, faPatchField, faMeshMapper, areaMesh>(mapper);

    // Map all the edgeFields in the objectRegistry
    MapGeometricFields<scalar, faePatchField, faMeshMapper, edgeMesh>(mapper);
    MapGeometricFields<vector, faePatchField, faMeshMapper, edgeMesh>(mapper);
    MapGeometricFields<sphericalTensor, faePatchField, faMeshMapper, edgeMesh>
        (mapper);
    MapGeometricFields<symmTensor, faePatchField, faMeshMapper, edgeMesh>
        (mapper);
    MapGeometricFields<tensor, faePatchField, faMeshMapper, edgeMesh>(mapper);
}


void Foam::faMesh::mapOldAreas(const faMeshMapper& mapper) const
{
    if (S0Ptr_)
    {
        if (debug)
        {
            InfoIn("void faMesh::mapOldAreas(const faMeshMapper& mapper)")
                << "Mapping old face areas." << endl;
        }

        scalarField& S0 = *S0Ptr_;

        scalarField savedS0(S0);
        S0.setSize(nFaces());

        const labelList& faceMap = mapper.areaMap().newFaceLabelsMap();

        // Map existing old areas; for new faces set area to zero
        forAll (faceMap, faceI)
        {
            if (faceMap[faceI] > -1)
            {
                S0[faceI] = savedS0[faceMap[faceI]];
            }
            else
            {
                S0[faceI] = 0;
            }
        }
    }

    if (S00Ptr_)
    {
        if (debug)
        {
            InfoIn("void faMesh::mapOldAreas(const faMeshMapper& mapper)")
                << "Mapping old-old face areas." << endl;
        }

        scalarField& S00 = *S00Ptr_;

        scalarField savedS00(S00);
        S00.setSize(nFaces());

        const labelList& faceMap = mapper.areaMap().newFaceLabelsMap();

        // Map old areas for existing faces; for new faces, set area to zero
        forAll (faceMap, faceI)
        {
            if (faceMap[faceI] > -1)
            {
                S00[faceI] = savedS00[faceMap[faceI]];
            }
            else
            {
                S00[faceI] = 0;
            }
        }
    }

}


// ************************************************************************* //
