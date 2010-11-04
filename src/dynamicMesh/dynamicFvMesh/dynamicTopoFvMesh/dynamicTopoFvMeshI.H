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

Class
    dynamicTopoFvMesh

Description
    An implementation of dynamic changes to mesh-topology

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Stack.H"
#include "meshOps.H"
#include "tetPointRef.H"
#include "linePointRef.H"
#include "lengthScaleEstimator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a non-const reference to the lengthScaleEstimator
inline lengthScaleEstimator&
dynamicTopoFvMesh::lengthEstimator()
{
    if (!lengthEstimator_.valid())
    {
        FatalErrorIn
        (
            "inline lengthScaleEstimator& "
            "dynamicTopoFvMesh::lengthEstimator()"
        ) << nl
          << " Invalid request. Length scale estimator was not loaded. "
          << abort(FatalError);
    }

    return lengthEstimator_();
}


// Return a const reference to the lengthScaleEstimator
inline const lengthScaleEstimator&
dynamicTopoFvMesh::lengthEstimator() const
{
    if (!lengthEstimator_.valid())
    {
        FatalErrorIn
        (
            "inline const lengthScaleEstimator& "
            "dynamicTopoFvMesh::lengthEstimator() const"
        ) << nl
          << " Invalid request. Length scale estimator was not loaded. "
          << abort(FatalError);
    }

    return lengthEstimator_();
}


// Return a const reference to the multiThreader
inline const multiThreader&
dynamicTopoFvMesh::threader() const
{
    if (!threader_.valid())
    {
        FatalErrorIn
        (
            "inline const multiThreader& "
            "dynamicTopoFvMesh::threader() const"
        ) << nl
          << " Invalid request. multiThreader was not loaded. "
          << abort(FatalError);
    }

    return threader_();
}


// Return a reference to the entity mutexes.
// The index 'entity' ranges from 0 to 3 for point/edge/face/cell.
inline const Mutex& dynamicTopoFvMesh::entityMutex
(
    const label entity
) const
{
    return entityMutex_[entity];
}


// Return the edge index for a provided edge
inline label dynamicTopoFvMesh::getEdgeIndex
(
    const edge& edgeToCheck
) const
{
    if (twoDMesh_)
    {
        // No efficient search method. Loop through all edges.
        forAll(edges_, edgeI)
        {
            if (edges_[edgeI] == edgeToCheck)
            {
                return edgeI;
            }
        }
    }
    else
    {
        // Look througg pointEdges list
        const labelList& pEdges = pointEdges_[edgeToCheck.start()];

        forAll(pEdges, edgeI)
        {
            if (edges_[pEdges[edgeI]] == edgeToCheck)
            {
                return pEdges[edgeI];
            }
        }
    }

    // Probably not a valid edge.
    FatalErrorIn
    (
        "inline label dynamicTopoFvMesh::getEdgeIndex"
        "(const edge& edgeToCheck) const"
    )
        << "Could not find an appropriate edge index for edge:"
        << edgeToCheck
        << abort(FatalError);

    return -1;
}


// Given any quad face, pick out a boundary edge that
// contains a triangular face. For 2D simplical meshes only.
inline label dynamicTopoFvMesh::getTriBoundaryEdge
(
    const label fIndex
) const
{
    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        // Obtain edgeFaces for this edge
        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            if (faces_[eFaces[faceI]].size() == 3)
            {
                // Found a triangular face. Return this edge.
                return fEdges[edgeI];
            }
        }
    }

    // This bit should never happen.
    FatalErrorIn
    (
        "inline label dynamicTopoFvMesh::getTriBoundaryEdge"
        "(const label fIndex) const"
    )
        << "Cannot find a triangular face bordering face: "
        << fIndex << " :: " << faces_[fIndex]
        << abort(FatalError);

    return -1;
}


// Return length-scale at an face-location in the mesh [2D]
inline scalar dynamicTopoFvMesh::faceLengthScale
(
    const label fIndex
) const
{
    // Reset the scale first
    scalar scale = 0.0;

    label facePatch = whichPatch(fIndex);

    // Determine whether the face is internal
    if (facePatch < 0)
    {
#       ifdef FULLDEBUG
        // Check whether neighbour is valid
        if (neighbour_[fIndex] == -1)
        {
            FatalErrorIn
            (
                "inline scalar dynamicTopoFvMesh::faceLengthScale"
                "(const label fIndex) const"
            )
                << nl << "Face: " << fIndex
                << ": " << faces_[fIndex]
                << " is not internal."
                << abort(FatalError);
        }
#       endif

        scale =
        (
            0.5 *
            (
                lengthScale_[owner_[fIndex]]
              + lengthScale_[neighbour_[fIndex]]
            )
        );
    }
    else
    {
        // Fetch the fixed-length scale
        scale = lengthEstimator().fixedLengthScale(fIndex, facePatch);

        // If this is a floating face, pick the owner length-scale
        if (lengthEstimator().isFreePatch(facePatch))
        {
            scale = lengthScale_[owner_[fIndex]];
        }

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (lengthEstimator().isProximityPatch(facePatch))
        {
            label proximityFace = -1;

            // Perform a proximity-check.
            scalar distance = testProximity(fIndex, proximityFace);

            if (debug > 3 && self() == 0)
            {
                if
                (
                    (proximityFace > -1) &&
                    ((distance / 5.0) < scale)
                )
                {
                    Info << " Closest opposing face detected for face: " << nl
                         << '\t' << fIndex
                         << " :: " << faces_[fIndex]
                         << " was face:\n"
                         << '\t' << proximityFace
                         << " :: " << polyMesh::faces()[proximityFace] << nl
                         << " with distance: " << distance
                         << endl;
                }
            }

            scale =
            (
                Foam::min
                (
                    scale,
                    ((distance / 3.0) - SMALL)/lengthEstimator().ratioMax()
                )
            );
        }

        // Limit scales if necessary
        lengthEstimator().limitScale(scale);
    }

    return scale;
}


// Compute length-scale at an edge-location in the mesh [3D]
inline scalar dynamicTopoFvMesh::edgeLengthScale
(
    const label eIndex
) const
{
    // Reset the scale first
    scalar scale = 0.0;

    const labelList& eFaces = edgeFaces_[eIndex];

    label edgePatch = whichEdgePatch(eIndex);

    // Determine whether the edge is internal
    if (edgePatch < 0)
    {
        forAll(eFaces, faceI)
        {
#           ifdef FULLDEBUG
            // Check whether neighbour is valid
            if (neighbour_[eFaces[faceI]] == -1)
            {
                FatalErrorIn
                (
                    "inline scalar dynamicTopoFvMesh::edgeLengthScale"
                    "(const label eIndex) const"
                )
                    << nl << "Face: " << eFaces[faceI]
                    << ": " << faces_[eFaces[faceI]]
                    << " is not internal, while edge: "
                    << eIndex << ": " << edges_[eIndex] << " is."
                    << abort(FatalError);
            }
#           endif

            scale += lengthScale_[owner_[eFaces[faceI]]];
            scale += lengthScale_[neighbour_[eFaces[faceI]]];
        }

        scale /= (2.0*eFaces.size());
    }
    else
    {
        // Search for boundary faces, and average their scale
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                scale +=
                (
                    lengthEstimator().fixedLengthScale
                    (
                        eFaces[faceI],
                        edgePatch
                    )
                );
            }
        }

        scale *= 0.5;

        // If proximity-based refinement is requested,
        // test the proximity to the nearest non-neighbour.
        if (lengthEstimator().isProximityPatch(edgePatch))
        {
            label proximityFace = -1;

            // Perform a proximity-check.
            scalar distance = testProximity(eIndex, proximityFace);

            if (debug > 3 && self() == 0)
            {
                if
                (
                    (proximityFace > -1) &&
                    ((distance / 5.0) < scale)
                )
                {
                    Info << " Closest opposing face detected for edge: " << nl
                         << '\t' << eIndex
                         << " :: " << edges_[eIndex]
                         << " was face:\n"
                         << '\t' << proximityFace
                         << " :: " << polyMesh::faces()[proximityFace] << nl
                         << " with distance: " << distance
                         << endl;
                }
            }

            scale =
            (
                Foam::min
                (
                    scale,
                    ((distance / 3.0) - SMALL)/lengthEstimator().ratioMax()
                )
            );
        }

        // If curvature-based refinement is requested,
        // test the variation in face-normal directions.
        if (lengthEstimator().isCurvaturePatch(edgePatch))
        {
            // Obtain face-normals for both faces.
            label count = 0;
            FixedList<vector, 2> fNorm;

            forAll(eFaces, faceI)
            {
                if (neighbour_[eFaces[faceI]] == -1)
                {
                    // Obtain the normal.
                    fNorm[count] = faces_[eFaces[faceI]].normal(points_);

                    // Normalize it.
                    fNorm[count] /= mag(fNorm[count]);

                    count++;
                }
            }

            scalar deviation = (fNorm[0] & fNorm[1]);
            scalar refDeviation = lengthEstimator().curvatureDeviation();

            if (mag(deviation) < refDeviation)
            {
                // Fetch the edge
                const edge& edgeToCheck = edges_[eIndex];

                // Get the edge-length.
                scalar length =
                (
                    linePointRef
                    (
                        points_[edgeToCheck.start()],
                        points_[edgeToCheck.end()]
                    ).mag()
                );

                if (debug > 3 && self() == 0)
                {
                    Info << "Deviation: " << deviation << nl
                         << "curvatureDeviation: " << refDeviation
                         << ", Edge: " << eIndex << ", Length: " << length
                         << ", Scale: " << scale << nl
                         << " Half-length: " << (0.5*length) << nl
                         << " MinRatio: "
                         << (lengthEstimator().ratioMin()*scale)
                         << endl;
                }

                scale =
                (
                    Foam::min
                    (
                        scale,
                        ((length - SMALL)/lengthEstimator().ratioMax())
                    )
                );
            }
        }

        // Limit scales if necessary
        lengthEstimator().limitScale(scale);
    }

    return scale;
}


// Check for edge bisection
inline bool dynamicTopoFvMesh::checkBisection
(
    const label index
) const
{
    scalar length = 0.0, scale = 0.0;

    if (twoDMesh_)
    {
        // Fetch the edge
        const edge& edgeToCheck = edges_[getTriBoundaryEdge(index)];

        // Measure the boundary edge-length of the face in question
        length =
        (
            linePointRef
            (
                points_[edgeToCheck.start()],
                points_[edgeToCheck.end()]
            ).mag()
        );

        // Determine the length-scale at this face
        scale = faceLengthScale(index);
    }
    else
    {
        // Fetch the edge
        const edge& edgeToCheck = edges_[index];

        // Measure the edge-length
        length =
        (
            linePointRef
            (
                points_[edgeToCheck.start()],
                points_[edgeToCheck.end()]
            ).mag()
        );

        // Determine the length-scale at this point in the mesh
        scale = edgeLengthScale(index);
    }

    if (length > lengthEstimator().ratioMax() * scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Check for edge collapse
inline bool dynamicTopoFvMesh::checkCollapse
(
    const label index
) const
{
    scalar length = 0.0, scale = 0.0;

    if (twoDMesh_)
    {
        // Fetch the edge
        const edge& edgeToCheck = edges_[getTriBoundaryEdge(index)];

        // Measure the boundary edge-length of the face in question
        length =
        (
            linePointRef
            (
                points_[edgeToCheck.start()],
                points_[edgeToCheck.end()]
            ).mag()
        );

        // Determine the length-scale at this face
        scale = faceLengthScale(index);
    }
    else
    {
        // Fetch the edge
        const edge& edgeToCheck = edges_[index];

        // Measure the edge-length
        length =
        (
            linePointRef
            (
                points_[edgeToCheck.start()],
                points_[edgeToCheck.end()]
            ).mag()
        );

        // Determine the length-scale at this point in the mesh
        scale = edgeLengthScale(index);
    }

    if (length < lengthEstimator().ratioMin() * scale)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Return the entity stack
inline Stack& dynamicTopoFvMesh::stack
(
    const label threadID
)
{
    return entityStack_[threadID];
}


// Return the integer ID for a given thread
// Return zero for single-threaded operation
inline label dynamicTopoFvMesh::self() const
{
    if (threader_->multiThreaded())
    {
        for (label i = 1; i <= threader_->getNumThreads(); i++)
        {
            if (handlerPtr_[i].self())
            {
                return i;
            }
        }
    }

    return 0;
}


// Initialize edge-stacks
inline void dynamicTopoFvMesh::initStacks
(
    const labelHashSet& entities
)
{
    forAll(entityStack_, stackI)
    {
        entityStack_[stackI].clear();
    }

    // Prepare a filling sequence based on threading operation
    label tIndex = 0;
    labelList tID(threader().getNumThreads());

    if (threader_->multiThreaded())
    {
        forAll(tID, tI)
        {
            tID[tI] = (tI + 1);
        }
    }
    else
    {
        tID = 0;
    }

    if (twoDMesh_)
    {
        forAll(faces_, faceI)
        {
            if (faces_[faceI].size() == 4)
            {
                stack(tID[tIndex]).insert(faceI);

                tIndex = tID.fcIndex(tIndex);
            }
        }
    }
    else
    {
        forAll(edges_, edgeI)
        {
            if (edgeFaces_[edgeI].size())
            {
                stack(tID[tIndex]).insert(edgeI);

                tIndex = tID.fcIndex(tIndex);
            }
        }
    }
}


// Method to determine the old boundary patch index for a given face
// Similar to the polyBoundaryMesh routine, but works on local information
inline label dynamicTopoFvMesh::whichPatch
(
    const label index
) const
{
    if (index < nOldInternalFaces_)
    {
        return -1;
    }

    for (label i = 0; i < boundaryMesh().size(); i++)
    {
        if
        (
            index >= oldPatchStarts_[i]
         && index < oldPatchStarts_[i] + oldPatchSizes_[i]
        )
        {
            return i;
        }
    }

    // If not in any of the above, it's possible that the face was added
    // at the end of the list. Check addedFacePatches_ for the patch info
    if (addedFacePatches_.found(index))
    {
        return addedFacePatches_[index];
    }
    else
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::whichPatch"
            "(const label index) const"
        )
            << "Cannot find patch information for face index: "
            << index << nl
            << " It appears that face ordering is"
            << " inconsistent with patch information." << nl
            << " Mesh info: " << nl
            << "    nOldInternalFaces: " << nOldInternalFaces_ << nl
            << "    oldPatchStarts: " << oldPatchStarts_ << nl
            << "    oldPatchSizes: " << oldPatchSizes_ << nl
            << abort(FatalError);
    }

    return -2;
}


// Method to determine the old boundary patch index for a given edge
inline label dynamicTopoFvMesh::whichEdgePatch
(
    const label index
) const
{
    if (index < nOldInternalEdges_)
    {
        return -1;
    }

    for (label i = 0; i < boundaryMesh().size(); i++)
    {
        if
        (
            index >= oldEdgePatchStarts_[i]
         && index < oldEdgePatchStarts_[i] + oldEdgePatchSizes_[i]
        )
        {
            return i;
        }
    }

    // If not in any of the above, it's possible that the edge was added
    // at the end of the list. Check addedEdgePatches_ for the patch info
    if (addedEdgePatches_.found(index))
    {
        return addedEdgePatches_[index];
    }
    else
    {
        FatalErrorIn
        (
            "label dynamicTopoFvMesh::whichEdgePatch"
            "(const label index) const"
        )
            << "Cannot find patch information for edge index: "
            << index << nl
            << " It appears that edge ordering is"
            << " inconsistent with patch information." << nl
            << " Mesh info: " << nl
            << "    nOldInternalEdges: " << nOldInternalEdges_ << nl
            << "    oldEdgePatchStarts: " << oldEdgePatchStarts_ << nl
            << "    oldEdgePatchSizes: " << oldEdgePatchSizes_ << nl
            << abort(FatalError);
    }

    return -2;
}


// Report topo-modification status
//  - Enumerants are as follows
//    [0] nModifications
//    [1] Total swaps (internal and surface)
//    [2] Surface swaps
//    [3] Total bisections
//    [4] Total collapses
//    [5] Surface bisections
//    [6] Surface collapses
//    [7] Slivers removed
inline label dynamicTopoFvMesh::status(const label type) const
{
    if (type < statistics_.size())
    {
        return statistics_[type];
    }
    else
    {
        FatalErrorIn
        (
            "inline label dynamicTopoFvMesh::status"
            "(const label type) const"
        )
            << " Unknown type: " << type << nl
            << abort(FatalError);
    }

    return 0;
}


// Set a particular face index as flipped.
inline void dynamicTopoFvMesh::setFlip(const label fIndex)
{
    if (fIndex < nOldFaces_)
    {
        if (flipFaces_.found(fIndex))
        {
            flipFaces_.erase(fIndex);
        }
        else
        {
            flipFaces_.insert(fIndex);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
