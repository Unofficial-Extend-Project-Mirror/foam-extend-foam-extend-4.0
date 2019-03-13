/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "ImmersedFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Distance>
void Foam::ImmersedFace<Distance>::createSubfaces
(
    const face& localFace,
    const scalarField& depth
)
{
    // Cut edges that cross the surface at the surface and add to
    // points and intersections

    // Make a copy of starting face points
    pointField localPoints(facePointsAndIntersections_);

    // Note: depth corresponds to local points

    // Expand the list for additional points. This leaves sufficient
    // space for intersection at every edge
    facePointsAndIntersections_.setSize(2*localPoints.size());
    scalarField newDepth(2*localPoints.size());

    // Get list of edges
    const edgeList edges = localFace.edges();

    // Count the number of newly created points, including original points
    label nNewPoints = 0;

    // Loop through all edges
    forAll (edges, edgeI)
    {
        // Take reference to currentEdge
        const edge& curEdge = edges[edgeI];

        const label start = curEdge.start();
        const label end = curEdge.end();

        const scalar edgeLength = curEdge.mag(localPoints);

        // Check if there is a legitimate cut to be found
        // Note: synced tolerances in ImmersedCell and ImmersedFace
        // HJ, 13/Mar/2019
        if
        (
            depth[start]*depth[end] < 0
         && edgeLength > SMALL
         && mag(depth[start]) > edgeLength*immersedPoly::tolerance_()
         && mag(depth[end]) > edgeLength*immersedPoly::tolerance_()
        )
        {
            // Prepare a new point to insert and determine its location
            point cutPoint;
            scalar depthAtCut = 0;

            if (!dist_.iterateDistance())
            {
                // Intersection is along the edge length (pf[end] - pf[start])
                // times the ratio of the depth at start and the difference
                // between depth at start and end; add to this the start point
                // and you have the location
                cutPoint =
                    localPoints[start]
                  + depth[start]/(depth[start] - depth[end])*
                    (localPoints[end] - localPoints[start]);
            }
            else
            {
                // Initialize bisection starting points
                point p0 = localPoints[start];
                point p1 = localPoints[end];

                // Depth at starting points
                scalar d0 = depth[start];
                scalar d1 = depth[end];

                // Initial guess of starting point same
                // as in non-iterative approach
                cutPoint = p0 + mag(d0)/(mag(d0) + mag(d1))*(p1 - p0);

                // Convergence criterion is the depth at newP
                depthAtCut = dist_.distance(cutPoint);

                // Initialize loop counter
                label iters = 0;

                while
                (
                    (mag(depthAtCut) > immersedPoly::tolerance_())
                 && (iters < immersedPoly::nIter_())
                )
                {
                    // Is the guessed point on the same side of the surface
                    // as p0? If yes, move p0 to the guessed point and thus
                    // shorten the interval
                    if (sign(depthAtCut) == sign(d0))
                    {
                        d0 = depthAtCut;
                        p0 = cutPoint;
                    }
                    // otherwise, shorten the other side
                    else
                    {
                        d1 = depthAtCut;
                        p1 = cutPoint;
                    }

                    // determine new intersection point
                    cutPoint =  p0 + mag(d0)/(mag(d0) + mag(d1))*(p1 - p0);

                    // and calculate its depth
                    depthAtCut = dist_.distance(cutPoint);

                    iters++;
                }

            }

            // Store first point of edge
            facePointsAndIntersections_[nNewPoints] =
                localPoints[curEdge.start()];

            // Store first point depth
            newDepth[nNewPoints] = depth[curEdge.start()];

            nNewPoints++;

            // Store the newly found cut point
            facePointsAndIntersections_[nNewPoints] = cutPoint;

            // Store newly found cut depth
            newDepth[nNewPoints] = depthAtCut;

            nNewPoints++;
        }
        else
        {
            // No intersection: just copy first point of edge
            facePointsAndIntersections_[nNewPoints] =
                localPoints[curEdge.start()];

            // Store first point depth
            newDepth[nNewPoints] = depth[curEdge.start()];

            nNewPoints++;
        }
    }

    // Point list should now be complete because last point of last edge should
    // be the starting point of the first edge
    facePointsAndIntersections_.setSize(nNewPoints);
    newDepth.setSize(nNewPoints);

    // Analyse new depth

    // For each point, determine if it is submerged( = -1), dry( = 1) or
    // on the surface ( = 0)
    labelField isSubmerged(newDepth.size());

    forAll (newDepth, pointI)
    {
        if (mag(newDepth[pointI]) < immersedPoly::tolerance_())
        {
            isSubmerged[pointI] = 0;
        }
        else
        {
            isSubmerged[pointI] = sign(newDepth[pointI]);
        }
    }

    // Determine if face is on surface, fully dry, fully submerged
    // or intersected by surface
    label sumMagSubmerged = sum(mag(isSubmerged));
    label sumSubmerged = sum(isSubmerged);

    if (sumSubmerged == sumMagSubmerged)
    {
        // Face fully dry
        drySubface_ = localFace;
        wetSubface_ = face();
    }
    else if (sumSubmerged == -sumMagSubmerged)
    {
        // Face fully wet
        drySubface_ = face();
        wetSubface_ = localFace;
    }
    else if (sumMagSubmerged == 0)
    {
        // Face fully on surface:
        // set both wet and dry face to the originial face
        drySubface_ = localFace;
        wetSubface_ = localFace;
    }
    else
    {
        // Face is intersected by surface

        // Initialise both faces to full size of intesection points
        // to be truncated after completion

        drySubface_.setSize(facePointsAndIntersections_.size());
        label nDry = 0;

        wetSubface_.setSize(facePointsAndIntersections_.size());
        label nWet = 0;

        forAll (facePointsAndIntersections_, pointI)
        {
            if (isSubmerged[pointI] == 1)
            {
                // Point is dry, add to dry sub-face
                drySubface_[nDry] = pointI;
                nDry++;
            }
            else if (isSubmerged[pointI] == -1)
            {
                // Point is submerged, add to wet sub-face
                wetSubface_[nWet] = pointI;
                nWet++;
            }
            else
            {
                // Point is on surface, add to both dry and wet sub-face
                drySubface_[nDry] = pointI;
                nDry++;

                wetSubface_[nWet] = pointI;
                nWet++;
            }
        }

        // Check if surface is merely touching the face
        // in that case, either dry or wet sub-face have less
        // than 3 points
        if (nDry < 3)
        {
            // The face is wet
            isAllWet_ = true;
            isAllDry_ = false;

            drySubface_.clear();
        }
        else
        {
            drySubface_.setSize(nDry);

            // // Check area: if it is very small, reset the face
            // if
            // (
            //     drySubface_.mag(facePointsAndIntersections_)
            //   < immersedPoly::liveFactor_()*
            //     localFace.mag(facePointsAndIntersections_)
            // )
            // {
            //     // The face is practically wet
            //     isAllWet_ = true;

            //     drySubface_.clear();
            // }
        }

        if (nWet < 3)
        {
            // The face is dry
            isAllWet_ = false;
            isAllDry_ = true;

            wetSubface_.clear();
        }
        else
        {
            wetSubface_.setSize(nWet);

            // // Check area: if it is very small, reset the face
            // // Note: must use relative tolerance, ie divide with original face
            // // area magnitude.  HJ, 30/Aug/2018
            // if
            // (
            //     wetSubface_.mag(facePointsAndIntersections_)
            //   < immersedPoly::liveFactor_()*
            //     localFace.mag(facePointsAndIntersections_)
            // )
            // {
            //     // The face is practically dry
            //     isAllDry_ = true;

            //     wetSubface_.clear();
            // }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Distance>
Foam::ImmersedFace<Distance>::ImmersedFace
(
    const label faceID,
    const polyMesh& mesh,
    const Distance& dist
)
:
    dist_(dist),
    wetSubface_(),
    drySubface_(),
    isAllWet_(false),
    isAllDry_(false)
{
    const face& origFace = mesh.faces()[faceID];

    // Store face points locally
    facePointsAndIntersections_ = origFace.points(mesh.points());
    face localFace(origFace.size());

    // Local face addresses into local points
    forAll (origFace, pointI)
    {
        localFace[pointI]  = pointI;
    }

    // Distance from the surface for every point of face
    scalarField depth = dist_.distance(facePointsAndIntersections_);

    // Check if all points are wet or dry, using absolute tolerance
    if (max(depth) < immersedPoly::tolerance_())
    {
        // All points are wet within a tolerance: face is wet
        isAllWet_ = true;
        isAllDry_ = false;

        wetSubface_ = localFace;
    }
    else if (min(depth) > -immersedPoly::tolerance_())
    {
        // All points are dry within a tolerance: face is dry
        isAllWet_ = false;
        isAllDry_ = true;

        drySubface_ = localFace;
    }
    else
    {
        // Face appears to be cut by the free surface.
        // Perform detailed analysis to create dry and wet sub-face
        createSubfaces(localFace, depth);
    }
}


// ************************************************************************* //
