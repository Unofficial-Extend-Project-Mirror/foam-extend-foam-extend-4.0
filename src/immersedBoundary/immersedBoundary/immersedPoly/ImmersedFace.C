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

    // For each point, determine if it is submerged( = -1), dry( = 1) or
    // on the surface ( = 0)
    // This is done during cutting to avoid using another tolerance check later
    // to dermine which points are on the surface, below or above it. By
    // definition, points that are a result of cutting are on the surface. (IG
    // 14/May/2019)
    labelList isSubmerged(facePointsAndIntersections_.size());

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

            // Determine whether it is above or below the surface.
            // NOTE: it must be one or the other since this is an original point
            // of the edge, and it passed the if statement above (IG
            // 14/May/2019)
            isSubmerged[nNewPoints] = sign(depth[curEdge.start()]);

            nNewPoints++;

            // Store the newly found cut point
            facePointsAndIntersections_[nNewPoints] = cutPoint;

            // Store newly found cut depth
            newDepth[nNewPoints] = depthAtCut;

            // The cut point is by definition on the surface and therefore
            // shared by the dry and wet face (IG 14/May/2019)
            isSubmerged[nNewPoints] = 0;

            nNewPoints++;
        }
        else
        {
            // No intersection: just copy first point of edge
            facePointsAndIntersections_[nNewPoints] =
                localPoints[curEdge.start()];

            // Store first point depth
            newDepth[nNewPoints] = depth[curEdge.start()];

            // Determine whether it is above, below or on the surface.
            // NOTE: now it can be any of the options since end or start is
            // sitting on the surface, othervise the if statement above would
            // have been true.(IG 14/May/2019)
            if
            (
                mag(depth[curEdge.start()])
              < edgeLength*immersedPoly::tolerance_()
            )
            {
                isSubmerged[nNewPoints] = 0;
            }
            else
            {
                isSubmerged[nNewPoints] = sign(depth[curEdge.start()]);
            }

            nNewPoints++;
        }
    }

    // Point list should now be complete because last point of last edge should
    // be the starting point of the first edge
    facePointsAndIntersections_.setSize(nNewPoints);
    newDepth.setSize(nNewPoints);
    isSubmerged.setSize(nNewPoints);

    // Count the number of points on wet and dry parts of the face and create
    // the faces
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

            // Since cell cut is adjusted, face cut cannot be.
            // HJ, 5/Apr/2019
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

            // Since cell cut is adjusted, face cut cannot be.
            // HJ, 5/Apr/2019
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

    // Calculating absolute tolerances based on minimum edge length
    scalar absTol = 0.0;
    {
        // Use local edges
        const edgeList edges = localFace.edges();

        // Calculate min edge length for a quick check
        scalar minEdgeLength = GREAT;

        // Note: expensive calculation of min length.  HJ, 28/May/2015
        forAll (edges, edgeI)
        {
            minEdgeLength =
                Foam::min
                (
                    minEdgeLength,
                    edges[edgeI].mag(facePointsAndIntersections_)
                );
        }

        absTol = minEdgeLength*immersedPoly::tolerance_();
    }

    // Check if all points are wet or dry, using absolute tolerance
    if (max(depth) < absTol)
    {
        // All points are wet within a tolerance: face is wet
        isAllWet_ = true;
        isAllDry_ = false;

        wetSubface_ = localFace;
    }
    else if (min(depth) > -absTol)
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
