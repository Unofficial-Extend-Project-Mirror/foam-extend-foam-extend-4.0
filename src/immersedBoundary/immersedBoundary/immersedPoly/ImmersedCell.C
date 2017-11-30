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

\*---------------------------------------------------------------------------*/

#include "ImmersedCell.H"
#include "plane.H"
#include "transform.H"
#include "SortableList.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Distance>
void Foam::ImmersedCell<Distance>::getBase
(
    const vector& n,
    vector& e0,
    vector& e1
) const
{
    // Copy from class: geomCellLooper

    // Guess for vector normal to n.
    vector base(1, 0, 0);

    scalar nComp = n & base;

    if (mag(nComp) > 0.8)
    {
        // Was bad guess. Try with different vector.

        base.x() = 0;
        base.y() = 1;

        nComp = n & base;

        if (mag(nComp) > 0.8)
        {
            base.y() = 0;
            base.z() = 1;

            nComp = n & base;
        }
    }

    // Use component normal to n as base vector.
    e0 = base - nComp*n;

    e0 /= mag(e0) + VSMALL;

    e1 = n ^ e0;
}


template<class Distance>
Foam::label Foam::ImmersedCell<Distance>::insertIntersectionPoints
(
    scalarField& depth
)
{
    // Get list of edges
    edgeList edges = this->edges();

    // Expandable list with additional points and depths
    DynamicList<point> extraPoints(edges.size());
    DynamicList<scalar> extraDepths(edges.size());

    // Loop through all edges
    forAll (edges, edgeI)
    {
        // Get reference to currentEdge
        const edge& curEdge = edges[edgeI];

        const label start = curEdge.start();
        const label end = curEdge.end();

        const scalar edgeLength = curEdge.mag(points_);

        // Check if there is a legitimate cut to be found
        if
        (
            depth[start]*depth[end] < 0
         && edgeLength > SMALL
         && mag(depth[start]) > absTol_
         && mag(depth[end]) > absTol_
        )
        {
            // Prepare a new point to insert
            point cutPoint;
            scalar depthAtCut = 0;

            // and determine its location
            if (!dist_.iterateDistance())
            {
                // Intersection is along the edge length (pf[end] - pf[start])
                // times the ratio of the depth at start and the difference
                // between depth at start and end; add to this the start point
                // and you have the location
                cutPoint =
                    points_[start]
                  + depth[start]/(depth[start] - depth[end])*
                    (points_[end] - points_[start]);
            }
            else
            {
                // Initialize bisection starting points
                point p0 = points_[start];
                point p1 = points_[end];

                // Depth at starting points
                scalar d0 = depth[start];
                scalar d1 = depth[end];

                // initial guess of starting point same
                // as in non-iterative approach
                cutPoint = p0 + mag(d0)/(mag(d0) + mag(d1))*(p1 - p0);

                // convergence criterion is the depth at newP
                depthAtCut = dist_.distance(cutPoint);

                // initialize loop counter
                label iters = 0;

                while
                (
                    (mag(depthAtCut) > absTol_)
                 && (iters < immersedPoly::nIter_())
                )
                {
                    // is the guessed point on the same side of the surface
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

                    // Determine new intersection point
                    cutPoint =  p0 + mag(d0)/(mag(d0) + mag(d1))*(p1 - p0);

                    // and calculate its depth
                    depthAtCut = dist_.distance(cutPoint);

                    iters++;
                }
            }

            // and store the newly found cut point
            extraPoints.append(cutPoint);
            extraDepths.append(depthAtCut);

            // Index of last inserted point in expanded point list
            label cutPointID = points_.size() + extraPoints.size() - 1;

            // Find faces connected to edge
            labelList edgeFaceIDs = this->edgeFaces()[edgeI];

            // Add the new point to each connected face,
            // but at the right position!
            forAll (edgeFaceIDs, edgeFaceI)
            {
                // Get reference edgeFace
                face& edgeFace = faces_[edgeFaceIDs[edgeFaceI]];

                // Find edge that is identical to the one we selected to cut
                edgeList edgeFaceEdges = edgeFace.edges();

                // New face point list will have one more point
                DynamicList<label> newFacePointLabels(edgeFace.size()+1);

                // therefore, loop through all edges of the current face
                forAll (edgeFaceEdges, edgeFaceEgdeI)
                {
                    if (edgeFaceEdges[edgeFaceEgdeI] == curEdge)
                    {
                        //if that edge is identical to the current edge
                        //insert both intersection and old edge end point
                        newFacePointLabels.append(cutPointID);
                        newFacePointLabels.append
                        (
                            edgeFaceEdges[edgeFaceEgdeI].end()
                        );

                    }
                    else
                    {
                        //if not, just insert the old edge end point
                        newFacePointLabels.append
                        (
                            edgeFaceEdges[edgeFaceEgdeI].end()
                        );
                    }
                }

                newFacePointLabels.shrink();

                // and replace with the face with an extra point added
                edgeFace = face(newFacePointLabels);
            }
        }
    }

    // Reduce memory usage of extraPoints to min required amount
    extraPoints.shrink();
    extraDepths.shrink();

    // Update points list and depths
    points_.append(extraPoints);
    depth.append(extraDepths);

    // Check for successful intersection: more than 3 added points
    // can form an internal face
    return extraPoints.size();
}


template<class Distance>
void Foam::ImmersedCell<Distance>::createInternalFace
(
    const label nIntersections
)
{
    // Sanity check: Do we have at least 3 intersection points?
    if (nIntersections < 3)
    {
         FatalErrorIn
         (
             "ImmersedCell::createInternalFace(const label nIntersections)"
         )   << "Less than 3 intersection points between cell and free surface"
             << abort(FatalError);
    }

    // Declare internal face with mixed-up point ordering
    face unorderedInternalFace(nIntersections);

    // Make a local list of intersection points
    pointField intersectionPoints(nIntersections);

    forAll (unorderedInternalFace, i)
    {
        label pointID = points_.size() - nIntersections + i;

        unorderedInternalFace[i] = pointID;

        intersectionPoints[i] = points_[pointID];
    }

    // Order points, so that they form a polygon
    // Algorithm in analogy to geomCellLooper.C

    // Calculate centre
    point centre = average(intersectionPoints);

    // Get base vectors of coordinate system normal
    // define plane that approximates the surface from 3 points

    // Line segment between points 0 and 1
    // Note: face orientation is unknown and needs to be adjusted
    // after the face has been created
    // HJ, 28/Nov/2017
    vector S0 = intersectionPoints[1] - intersectionPoints[0];
    S0 /= mag(S0) + SMALL;

    label pointID = -1;
    scalar minDotProd = 1 - SMALL;

    // Take best non-colinear value
    for (label pI = 2; pI < intersectionPoints.size(); pI++)
    {
        // Create second line segment
        vector S1 = intersectionPoints[pI] - intersectionPoints[0];
        S1 /= mag(S1) + SMALL;
        scalar curDotProd = mag(S0 & S1);

        if (curDotProd < minDotProd)
        {
            pointID = pI;
            minDotProd = curDotProd;
        }
    }

    if (pointID == -1)
    {
        // All intersection points are colinear
        return;
    }

    // Now create surface
    plane surface
    (
        intersectionPoints[0],
        intersectionPoints[1],
        intersectionPoints[pointID]
    );

    vector e0, e1;
    getBase(surface.normal(), e0, e1);

    // Get sorted angles from point on loop to centre of loop.
    SortableList<scalar> sortedAngles(intersectionPoints.size());

    forAll (sortedAngles, angleI)
    {
        vector toCentre(intersectionPoints[angleI] - centre);
        toCentre /= mag(toCentre);

        sortedAngles[angleI] = pseudoAngle(e0, e1, toCentre);
    }
    sortedAngles.sort();

    // Re-order points
    const labelList& indices = sortedAngles.indices();

    face orderedInternalFace(intersectionPoints.size());

    forAll (indices, i)
    {
        orderedInternalFace[i] = unorderedInternalFace[indices[i]];
    }

    // Put cut face on front of faces_ since it is the only internal face

    faceList externalFaces = faces_;

    faces_.setSize(externalFaces.size() + 1);
    faces_[0] = orderedInternalFace;

    // Copy the outside faces into the list
    for (label i = 0; i < externalFaces.size(); i++)
    {
        faces_[i + 1] = externalFaces[i];
    }
}


template<class Distance>
void Foam::ImmersedCell<Distance>::splitFace
(
    const label faceID,
    const scalarField& depth
)
{
    // Get old face with inserted point(s)
    const face& oldFace = faces_[faceID];

    // now make two faces: wet and dry
    // Wet face: wet points and intersection points
    face wetFace(oldFace.size());
    label nWet = 0;

    // Dry face: dry points and intersection points
    face dryFace(oldFace.size());
    label nDry = 0;

    forAll (oldFace, pointI)
    {
        if (mag(depth[oldFace[pointI]]) < absTol_)
        {
            // Point close to surface => intersection point, add to both lists
            wetFace[nWet] = oldFace[pointI];
            nWet++;

            dryFace[nDry] = oldFace[pointI];
            nDry++;
        }
        else if (depth[oldFace[pointI]] < 0)
        {
            // Point is submerged, add to wetFace
            wetFace[nWet] = oldFace[pointI];
            nWet++;
        }
        else
        {
            // Otherwise point must be dry, add to dryFace
            dryFace[nDry] = oldFace[pointI];
            nDry++;
        }
    }

    // Sanity check
    if (nWet >= 3 && nDry >= 3)
    {
        // Two valid faces are produced.  Insert them into the list

        // Shrink faces to minimum size
        wetFace.setSize(nWet);
        dryFace.setSize(nDry);

        // Replace original face by wet face
        faces_[faceID] = wetFace;

        // Append dry face as new face to faces
        faces_.setSize(faces_.size() + 1);
        faces_[faces_.size() - 1] = dryFace;
    }
    // else face cut has failed.  Do nothing
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Distance>
Foam::ImmersedCell<Distance>::ImmersedCell
(
    const label cellID,
    const polyMesh& mesh,
    const Distance& dist
)
:
    primitiveMesh
    (
        mesh.cells()[cellID].labels(mesh.faces()).size(), // nPoints
        0,                              // nInternalFaces (init to zero)
        mesh.cells()[cellID].size(),    // nFaces
        1                               // nCells
    ),
    cellID_(cellID),
    mesh_(mesh),
    dist_(dist),
    absTol_(0),
    isAllWet_(false),
    isAllDry_(false),
    // Initialize points_ with  points from cell
    points_(mesh_.cells()[cellID_].points(mesh_.faces(), mesh_.points())),
    faces_(),
    // We start with single cell with ID = 0, so it owns all faces
    faceOwner_(faces_.size(), 0),
    faceNeighbour_()
{
    const cell& origCell = mesh.cells()[cellID];

    // Distance from the surface for every point of cell
    // Depth will be modified with distance for added points
    // HJ, 28/Nov/2017
    scalarField depth = dist_.distance(points_);

    // Calculating absolute tolerances based on minimum edge length
    // Note: use original cell for calculation because the mesh is not
    // complete until faces are set.  HJ, 21/Jan/2014
    const edgeList cellEdges = origCell.edges(mesh.faces());

    // Calculate min edge length for a quick check
    scalar minEdgeLength = GREAT;

    // Note: expensive calculation of min length.  HJ, 28/May/2015
    forAll (cellEdges, edgeI)
    {
        minEdgeLength =
            Foam::min(minEdgeLength, cellEdges[edgeI].mag(mesh.points()));
    }

    absTol_ = minEdgeLength*immersedPoly::tolerance_();

    // Check if we have to perform cut at all
    if (max(depth) < absTol_)
    {
        // All points of cell are below water surface
        isAllWet_ = true;

        return;
    }
    else if (min(depth) > -absTol_)
    {
        // All points are above water surface
        isAllDry_ = true;

        return;
    }

    // Create hash table that maps points on global mesh to local point list
    HashTable<label, label, Hash<label> > pointMapTable(points_.size());

    labelList origCellPointLabels = origCell.labels(mesh.faces());

    forAll (points_, pointI)
    {
        // Insert globalID and localID
        pointMapTable.insert(origCellPointLabels[pointI], pointI);
    }

    // Make local face list by remapping the faces of the cell
    faces_ = faceList(origCell.size());

    forAll (origCell, faceI)
    {
        // Get old point list of faceI
        face origFace(mesh.faces()[origCell[faceI]]);

        // Make sure that all faces point outward,
        // since they are going to be outside cells
        if (!(mesh.faceOwner()[origCell[faceI]] == cellID))
        {
            // Cell is not owner of face, revert face orientation
            origFace = origFace.reverseFace();
        }

        // Make list to store new points
        labelList newLabels(origFace.size());

        // Map labels
        forAll (origFace, facePointI)
        {
            newLabels[facePointI] =
                pointMapTable.find(origFace[facePointI])();
        }

        // Create face from new point labels
        faces_[faceI] = face(newLabels);
    }

    /*********************************************************************/
    // Starting to modify the 1-cell primitiveMesh,
    // beyond this point be sure to know what points_, faces_, etc. contain,
    // before calling inherited primitiveMesh functions.
    // Here be dragons!
    /*********************************************************************/

    // Cut all edges that are intersected by the zero distance surface:
    // Add cutting points to points_
    // Add cutting points to faces connected to edge (will be reordered later)

#   ifdef WET_DEBUG
    Info << "Cell ID: " << cellID << "  BEFORE" << nl
        << "points: " << points_ << nl
        << "faces: " << faces_ << nl
        << "depth: " << depth << endl;
#   endif

    // Insert intersection points and adjust depth for intersections
    label nIntersections = insertIntersectionPoints(depth);

    // Recheck, if there has been a successful cut at all
    if (nIntersections < 3)
    {
        // Check if improvised cell centre is wet or dry
        if (dist_.distance(average(points_)) < 0)
        {
            // All points of cell are below water surface
            isAllWet_ = true;

            return;
        }
        else
        {
            // All points are above water surface
            isAllDry_ = true;

            return;
        }
    }

    // From here on, there exists a valid intersection

#   ifdef WET_DEBUG
    Info << "Cell ID: " << cellID << "  AFTER" << nl
        << "points: " << points_ << nl
        << "faces: " << faces_ << nl
        << "depth: " << depth << endl;
#   endif

    // For all faces with inserted points, do face splitting
    forAll (origCell, oldFaceI)
    {
        const face& oldFace = mesh.faces()[origCell[oldFaceI]];
        const face& newFace = faces_[oldFaceI];

        if (newFace.size() != oldFace.size())
        {
            splitFace(oldFaceI, depth);
        }
    }

    // Assign owners and neighbors to faces
    faceOwner_.setSize(faces_.size());

    forAll (faces_, faceI)
    {
        if (dist_.distance(faces_[faceI].centre(points_)) > 0)
        {
            // Face is dry
            faceOwner_[faceI] = DRY;
        }
        else
        {
            // Face is wet
            faceOwner_[faceI] = WET;
        }
    }

    // If we are not merely touching the water surface
    // with one point or edge, insert internal face that
    // connects all intersection points
    // create internal face, which gets inserted at front of faces_
    createInternalFace(nIntersections);

    // ... and add owner of new face at beginning of faceOwner_
    // Additional faces will be added after the internal face
    labelList faceOwnerExternal = faceOwner_;
    faceOwner_.setSize(faceOwnerExternal.size() + 1);
    faceNeighbour_.setSize(1);

    // Face needs to point out of the wet cell.  Make the wet cell its owner
    faceOwner_[0] = WET;
    faceNeighbour_[0] = DRY;

    forAll (faceOwnerExternal, i)
    {
        faceOwner_[i + 1] = faceOwnerExternal[i];
    }

    // Set ownership depending on internal face orientation
    // Create a point just below the face and find out if it is wet or dry
    const vector cutFaceCentre = faces_[0].centre(points_);
    const vector cutFaceNormal = faces_[0].normal(points_);

    // Calculate sampling point below the face centre using normal
    // Normal is sized with face area.  Possible tolerance issue for
    // very small faces.  Reconsider
    // HJ, 28/Nov/2017
    const vector samplePoint = cutFaceCentre - cutFaceNormal;

    if (dist_.distance(samplePoint) > 0)
    {
        // Face pointing out of dry cell.  Turn it
        faces_[0] = faces_[0].reverseFace();
    }

    // Update primitiveMesh parameters
    this->reset
    (
        points_.size(),                 // nPoints
        faceNeighbour_.size(),          // nInternalFaces
        faces_.size(),                  // nFaces
        faceNeighbour_.size() + 1       // nCells
    );

    // Recheck, if the cut cell has significant volume.  If not, reset it
    scalar wetCut = cellVolumes()[WET]/mesh_.cellVolumes()[cellID_];

    if (wetCut > 1 - immersedPoly::tolerance_())
    {
        // The cell is practically wet
        isAllWet_ = true;
    }

    if (wetCut < immersedPoly::tolerance_())
    {
        // The cell is practically dry
        isAllDry_ = true;
    }
}


// ************************************************************************* //
