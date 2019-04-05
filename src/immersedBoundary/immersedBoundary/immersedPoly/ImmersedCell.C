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
void Foam::ImmersedCell<Distance>::insertIntersectionPoints()
{
    // Get list of edges
    const edgeList& edges = this->edges();

    // Get edge-face addressing
    const labelListList& edgeFaces = this->edgeFaces();

    // There may be an extra point on every edge.  Resize the list of points
    const label oldSize = points_.size();
    points_.setSize(oldSize + edges.size());
    label nPoints = oldSize;

    // Loop through all edges
    forAll (edges, edgeI)
    {
        // Get reference to currentEdge
        const edge& curEdge = edges[edgeI];

        const label start = curEdge.start();
        const label end = curEdge.end();

        const scalar edgeLength = mag(points_[end] - points_[start]);

        // Check if there is a legitimate cut to be found
        // Note: synced tolerances in ImmersedCell and ImmersedFace
        // HJ, 13/Mar/2019
        if
        (
            depth_[start]*depth_[end] < 0
         && edgeLength > SMALL
         && mag(depth_[start]) > edgeLength*immersedPoly::tolerance_()
         && mag(depth_[end]) > edgeLength*immersedPoly::tolerance_()
        )
        {
            // Prepare a new point to insert
            point cutPoint;
            scalar depthAtCut = 0;

            // Intersection is along the edge length (pf[end] - pf[start])
            // times the ratio of the depth at start and the difference
            // between depth at start and end; add to this the start point
            // and you have the location
            cutPoint =
                points_[start]
              + depth_[start]/(depth_[start] - depth_[end])*
                (points_[end] - points_[start]);

            // Execute iterative cut if necessary
            if (dist_.iterateDistance())
            {
                // Initialize bisection starting points
                point p0 = points_[start];
                point p1 = points_[end];

                // Depth at starting points
                scalar d0 = depth_[start];
                scalar d1 = depth_[end];

                // Convergence criterion is the depth at newP
                depthAtCut = dist_.distance(cutPoint);

                // initialize loop counter
                label iters = 0;

                while
                (
                    (mag(depthAtCut) > immersedPoly::tolerance_())
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
                    // Otherwise, shorten the other side
                    else
                    {
                        d1 = depthAtCut;
                        p1 = cutPoint;
                    }

                    // Determine new intersection point and its depth
                    cutPoint = p0 + mag(d0)/(mag(d0) + mag(d1))*(p1 - p0);

                    depthAtCut = dist_.distance(cutPoint);

                    iters++;
                }
            }

            // Store the newly found cut point
            points_[nPoints] = cutPoint;

            // Find faces connected to edge
            const labelList& edgeFaceIDs = edgeFaces[edgeI];

            // Add the new point to each connected face at the right position!
            forAll (edgeFaceIDs, edgeFaceI)
            {
                // Get old face
                const face& oldFace = faces_[edgeFaceIDs[edgeFaceI]];

                // Make new face with one extra label
                face newFace(oldFace.size() + 1);

                // Count points added to new face
                label nfp = 0;

                // Loop through old face.  If this edge is found, add the
                // cut point label into the edge
                forAll (oldFace, fpI)
                {
                    // Add the point
                    newFace[nfp] = oldFace[fpI];
                    nfp++;

                    const label curPoint = oldFace[fpI];
                    const label nextPoint = oldFace.nextLabel(fpI);

                    if
                    (
                        (curPoint == start && nextPoint == end)
                     || (curPoint == end && nextPoint == start)
                    )
                    {
                        // Found the edge.  Inser the point
                        newFace[nfp] = nPoints;
                        nfp++;
                    }
                }

                // Debug: check if point insertion was successful
                if (nfp < newFace.size())
                {
                    FatalErrorInFunction
                        << "badInsertion"
                        << abort(FatalError);
                }

                faces_[edgeFaceIDs[edgeFaceI]] = newFace;
            }

            // Finished point insertion
            nPoints++;
        }
    }

    // Resize the points list
    points_.setSize(nPoints);

    // Extra depths are all zero
    depth_.setSize(nPoints);

    // For all cut points set depth to exactly zero
    for (label i = oldSize; i < depth_.size(); i++)
    {
        depth_[i] = 0;
    }
}


template<class Distance>
Foam::face Foam::ImmersedCell<Distance>::createInternalFace() const
{
    // Declare internal face with mixed-up point ordering
    face unorderedInternalFace(points_.size());

    // Collect all points with zero distance to surface
    label nPif = 0;

    forAll (depth_, pointI)
    {
        if (mag(depth_[pointI]) < absTol_)
        {
            // Found point on zero plane
            unorderedInternalFace[nPif] = pointI;
            nPif++;
        }
    }

    unorderedInternalFace.setSize(nPif);

    // Sanity check: Do we have at least 3 points at zero distance?
    if (nPif < 3)
    {
         FatalErrorInFunction
             << "Less than 3 intersection points in cell on free surface." << nl
             << "depth: " << depth_
             << abort(FatalError);
    }

    // Order points, so that they form a polygon
    // Algorithm in analogy to geomCellLooper.C

    // Calculate centre
    point centre = average(unorderedInternalFace.points(points_));

    // Get base vectors of coordinate system normal
    // define plane that approximates the surface from 3 points

    // Line segment between points 0 and 1
    // Note: face orientation is unknown and needs to be adjusted
    // after the face has been created
    // HJ, 28/Nov/2017
    vector S0 =
        points_[unorderedInternalFace[1]]
      - points_[unorderedInternalFace[0]];

    S0 /= mag(S0) + SMALL;

    label pointID = -1;
    scalar minDotProd = 1 - SMALL;

    // Take best non-colinear value
    for (label pI = 2; pI < unorderedInternalFace.size(); pI++)
    {
        // Create second line segment
        vector S1 =
            points_[unorderedInternalFace[pI]]
          - points_[unorderedInternalFace[0]];

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
        FatalErrorInFunction
            << "Colinear points in cut"
            << abort(FatalError);
    }

    // Now create surface
    plane surface
    (
        points_[unorderedInternalFace[0]],
        points_[unorderedInternalFace[1]],
        points_[unorderedInternalFace[pointID]]
    );

    vector e0, e1;
    getBase(surface.normal(), e0, e1);

    // Get sorted angles from point on loop to centre of loop.
    SortableList<scalar> sortedAngles(unorderedInternalFace.size());

    forAll (sortedAngles, angleI)
    {
        vector toCentre(points_[unorderedInternalFace[angleI]] - centre);
        toCentre /= mag(toCentre);

        sortedAngles[angleI] = pseudoAngle(e0, e1, toCentre);
    }
    sortedAngles.sort();

    // Re-order points
    const labelList& indices = sortedAngles.indices();

    face orderedInternalFace(unorderedInternalFace.size());

    forAll (indices, i)
    {
        orderedInternalFace[i] = unorderedInternalFace[indices[i]];
    }

    // Check direction of the new face using average wet and dry point
    // HJ, 5/Dec/2017
    point wetPoint = vector::zero;
    label nWet = 0;

    point dryPoint = vector::zero;
    label nDry = 0;

    label nUndecided = 0;

    forAll (depth_, i)
    {
        if (depth_[i] > absTol_)
        {
            dryPoint += points_[i];
            nDry++;
        }
        else if  (depth_[i] < -absTol_)
        {
            wetPoint += points_[i];
            nWet++;
        }
        else
        {
            nUndecided++;
        }
    }

    if (nUndecided == depth_.size())
    {
         FatalErrorInFunction
             << "All points lay on the tri surface, zero volume cell?"
             << nl << "Points: " << points_
             << abort(FatalError);
    }

    wetPoint /= nWet;
    dryPoint /= nDry;

    // Good direction points out of the wet cell
    vector dir = dryPoint - wetPoint;
    dir /= mag(dir) + SMALL;

    vector n = orderedInternalFace.normal(points_);
    n /= mag(n);

    if ((dir & n) < 0)
    {
        orderedInternalFace = orderedInternalFace.reverseFace();
    }

    // Note: the face may have wrong orientation here.  It is corrected later
    // HJ, 5/Dec/2017
    return orderedInternalFace;
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
    isBadCut_(false),
    // Initialize points_ with  points from cell
    points_(mesh_.cells()[cellID_].points(mesh_.faces(), mesh_.points())),
    faces_(),
    // We start with single cell with ID = 0, so it owns all faces
    faceOwner_(faces_.size(), 0),
    faceNeighbour_(),
    depth_(dist.distance(points_))
{
    const cell& origCell = mesh_.cells()[cellID];

    // Build a valid 1-cell mesh in local addressing

    // Create hash table that maps points on global mesh to local point list
    HashTable<label, label, Hash<label> > pointMapTable(points_.size());

    labelList origCellPointLabels = origCell.labels(mesh_.faces());

    forAll (points_, pointI)
    {
        // Insert globalID and localID
        pointMapTable.insert(origCellPointLabels[pointI], pointI);
    }

    // Make local face list by remapping the faces of the cell
    // Maximum number of new faces is twice the number of original faces
    // plus one internal face
    faces_ = faceList(origCell.size());

    forAll (origCell, faceI)
    {
        // Get old point list of faceI
        face origFace(mesh_.faces()[origCell[faceI]]);

        // Make sure that all faces point outward,
        // since they are going to be outside cells
        if (!(mesh_.faceOwner()[origCell[faceI]] == cellID))
        {
            // Cell is not owner of face, revert face orientation
            // for the use in a 1-cell mesh
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

    // At this point, a 1-cell mesh is valid

    // Calculating absolute tolerances based on minimum edge length
    {
        // Use local edges
        const edgeList& cellEdges = edges();

        // Calculate min edge length for a quick check
        scalar minEdgeLength = GREAT;

        // Note: expensive calculation of min length.  HJ, 28/May/2015
        forAll (cellEdges, edgeI)
        {
            minEdgeLength =
                Foam::min(minEdgeLength, cellEdges[edgeI].mag(points_));
        }

        absTol_ = minEdgeLength*immersedPoly::tolerance_();
    }

    // Check if we have to perform cut at all
    if (max(depth_) < absTol_)
    {
        // All points of cell are below water surface
        isAllWet_ = true;

        return;
    }
    else if (min(depth_) > -absTol_)
    {
        // All points are above water surface
        isAllDry_ = true;

        return;
    }

#   ifdef WET_DEBUG
    Info << "Cell ID: " << cellID << "  BEFORE" << nl
        << "points: " << points_ << nl
        << "faces: " << faces_ << nl
        << "depth: " << depth_ << endl;
#   endif

    /*********************************************************************/
    // Starting to modify the 1-cell primitiveMesh.
    // Beyond this point be sure to know what points_, faces_, etc. contain,
    // before calling inherited primitiveMesh functions of this class.
    // Here be dragons!
    /*********************************************************************/

    // Created expanded point and face lists

    // Insert intersection points and adjust depth for intersections
    // This will add further points into the intersected face if needed
    // Depth at intersection will be zero.  HJ, 5/Dec/2017
    // Note that it is possible to have the cut face even if no new points
    // have been introduced.  HJ, 13/Mar/2019
    insertIntersectionPoints();

    // Update primitiveMesh parameters
    this->reset
    (
        points_.size(),         // nPoints
        0,                      // nInternalFaces
        faces_.size(),          // nFaces
        1                       // nCells
    );

#   ifdef WET_DEBUG
    Info << "Cell ID: " << cellID << "  ENRICHED" << nl
        << "points: " << points_ << nl
        << "faces: " << faces_ << nl
        << "depth: " << depth_ << endl;
#   endif
    // At this point, a 1-cell mesh with faces enriched for intersections
    // is valid.  HJ, 5/Dec/2017

    // Check if there has been a successful cut at all
    // For a good cut there should be at least 3 points at zero level
    label nIntersections = 0;

    forAll (depth_, pointI)
    {
        if (mag(depth_[pointI]) < absTol_)
        {
            nIntersections++;
        }
    }

    if (nIntersections < 3)
    {
        // Check if cell centre is wet or dry, depending on greatest distance
        // away from the cutting surface
        // Note: cannot measure  distance geometrically because of
        // the unknown resolution of the immersed surface
        // HJ, 5/Dec/2017
        if (mag(min(depth_)) > mag(max(depth_)))
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

    // Resize the face list.  Each face can be split into two, with one
    // extra internal face.  HJ, 5/Dec/2017

    // Make a copy of enriched faces, on which the cutting is performed
    faceList enrichedFaces = faces_;

    // Reset face lists, preserving existing faces
    faces_.setSize(2*faces_.size() + 1);
    faceOwner_.setSize(2*faces_.size() + 1);
    faceNeighbour_.setSize(1);

    // If we are not merely touching the water surface
    // with one point or edge, insert internal face that
    // connects all intersection points
    // create internal face, which gets inserted at front of faces_ list
    faces_[0] = createInternalFace();

    // Internal face points out of the wet cell.  Make the wet cell its owner
    faceOwner_[0] = WET;
    faceNeighbour_[0] = DRY;

    // Count new faces
    label nFaces = 1;

    // For all faces with inserted points, do face splitting
    forAll (enrichedFaces, oldFaceI)
    {
        const face& oldFace = mesh_.faces()[origCell[oldFaceI]];
        const face& newFace = enrichedFaces[oldFaceI];

        // Calculate old face area locally to avoid triggering polyMesh
        const scalar oldFaceArea = mag(mesh_.faceAreas()[origCell[oldFaceI]]);

        // If a face has been modified, it will have extra points
        if (newFace.size() != oldFace.size())
        {
            // Make two faces: wet and dry
            // Wet face: wet points and intersection points
            face wetFace(newFace.size());
            label nWet = 0;

            // Dry face: dry points and intersection points
            face dryFace(newFace.size());
            label nDry = 0;

            forAll (newFace, pointI)
            {
                if (mag(depth_[newFace[pointI]]) < absTol_)
                {
                    // Intersection point. Add to both faces
                    wetFace[nWet] = newFace[pointI];
                    nWet++;

                    dryFace[nDry] = newFace[pointI];
                    nDry++;
                }
                else if (depth_[newFace[pointI]] < -absTol_)
                {
                    // Point is submerged, add to wetFace
                    wetFace[nWet] = newFace[pointI];
                    nWet++;
                }
                else // depth_[newFace[pointI]] > absTol_
                {
                    // Otherwise point must be dry, add to dryFace
                    dryFace[nDry] = newFace[pointI];
                    nDry++;
                }
            }

            // Check for a successful cut
            if (nWet >= 3)
            {
                // Insert wet face
                wetFace.setSize(nWet);
                faces_[nFaces] = wetFace;
                faceOwner_[nFaces] = WET;

                nFaces++;

                // Check for bad wet face cut
                if
                (
                    wetFace.mag(points_)
                  > (1 + immersedPoly::badCutFactor_())*oldFaceArea
                )
                {
                    // Wet face area is greater than original face area
                    // This is a bad cut
                    Info<< "Bad cell face cut: wet = ("
                        << wetFace.mag(points_) << " "
                        << oldFaceArea
                        << ")" << endl;

                    isBadCut_ = true;
                }
            }

            if (nDry >= 3)
            {
                // Insert dry face
                dryFace.setSize(nDry);
                faces_[nFaces] = dryFace;
                faceOwner_[nFaces] = DRY;

                nFaces++;

                // Check for bad dry face cut
                if
                (
                    dryFace.mag(points_)
                  > (1 + immersedPoly::badCutFactor_())*oldFaceArea
                )
                {
                    // Dry face area is greater than original face area
                    // This is a bad cut
                    Info<< "Bad cell face cut: dry = ("
                        << dryFace.mag(points_) << " "
                        << oldFaceArea
                        << ")" << endl;

                    isBadCut_ = true;
                }
            }
        }
        else
        {
            // Face cut has failed.  Insert original face and owner
            faces_[nFaces] = newFace;

            // Determine wet/dry based on distance to face centre
            // Note: cannot measure  distance geometrically because of
            // the unknown resolution of the immersed surface
            // HJ, 5/Dec/2017

            // Create face depth distance as a subset
            scalarField faceDepth(depth_, newFace);

            // Since the face has not been cut, all faceDepth should have the
            // same sign.  Otherwise, the face should straddle the immersed
            // surface. Check on minimum.
            // Note: this is a very precise check on purpose: there is no cut
            // and the face belongs either to a wet cell or a dry cell
            // HJ, 12/Mar/2019
            if (min(faceDepth) < scalar(0))
            {
                // Negative distance: wet face
                faceOwner_[nFaces] = WET;
            }
            else
            {
                // Positive distance: dry face
                faceOwner_[nFaces] = DRY;
            }

            nFaces++;
        }
    }

    faces_.setSize(nFaces);
    faceOwner_.setSize(nFaces);

    // Update primitiveMesh parameters
    this->reset
    (
        points_.size(),                 // nPoints
        faceNeighbour_.size(),          // nInternalFaces
        faces_.size(),                  // nFaces
        faceNeighbour_.size() + 1       // nCells
    );


#   ifdef WET_DEBUG
    this->checkMesh();
    Info << "Cell ID: " << cellID << "  AFTER" << nl
        << "points: " << points_ << nl
        << "faces: " << faces_ << nl
        << "depth: " << depth_ << endl;
#   endif

    const scalar oldCellVolume = mesh_.cellVolumes()[cellID_];

    // Note: is it legal to cut a zero volume cell?  HJ, 11/Mar/2019

    scalar wetCut = cellVolumes()[WET]/oldCellVolume;

    scalar dryCut = cellVolumes()[DRY]/oldCellVolume;

    // Check for bad cell cut based on volume
    if
    (
        wetCut < -immersedPoly::badCutFactor_()
     || wetCut > (1 + immersedPoly::badCutFactor_())
     || dryCut < -immersedPoly::badCutFactor_()
     || dryCut > (1 + immersedPoly::badCutFactor_())
    )
    {
        isBadCut_ = true;
    }

    //  If the cut is not bad, adjust the cell for thin cell cut
    if (!isBadCut_)
    {
        if (mag(wetCut) < immersedPoly::liveFactor_())
        {
            // Cell is dry; reset
            isAllDry_ = true;
        }

        if (mag(dryCut) < immersedPoly::liveFactor_())
        {
            // Cell is wet; reset
            isAllWet_ = true;
        }
    }
    else
    {
        Info<< "Bad cell cut: volume = (" << wetCut << " " << dryCut
            << ") = " << wetCut + dryCut << nl
            // << "Points: " << nl << this->points() << nl
            // << "Faces: " << nl << this->faces() << nl
            // << "Owner: " << nl << this->faceOwner() << nl
            // << "Neighbour: " << nl << this->faceNeighbour() << nl
            // << "Cut (wet dry) = (" << isAllWet_ << " " << isAllDry_ << ")"
            << endl;
    }

    // Correction on cutting is not allowed, as it results in an open cell
    // if faces are cut and the cell is not.
    // Previous check confirmed more than 3 valid cut points in the cell,
    // which means that some of the faces were cut.
    // Cutting tolerances for the cell and face have been adjusted to make sure
    // identical cut has been produced.
    // HJ, 11/Mar/2019
}


// ************************************************************************* //
