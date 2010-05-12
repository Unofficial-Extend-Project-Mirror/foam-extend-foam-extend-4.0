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

#include "polyMesh.H"
#include "meshSearch.H"
#include "triPointRef.H"
#include "octree.H"
#include "pointIndexHit.H"
#include "DynamicList.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::meshSearch, 0);

Foam::scalar Foam::meshSearch::tol_ = 1E-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// tree based searching
Foam::label Foam::meshSearch::findNearestCellTree(const point& location) const
{
    treeBoundBox tightest(treeBoundBox::greatBox);

    scalar tightestDist = GREAT;

    return cellCentreTree().findNearest(location, tightest, tightestDist);
}


// linear searching
Foam::label Foam::meshSearch::findNearestCellLinear(const point& location) const
{
    const vectorField& centres = mesh_.cellCentres();

    label nearestCelli = 0;
    scalar minProximity = magSqr(centres[0] - location);

    forAll(centres, celli)
    {
        scalar proximity = magSqr(centres[celli] - location);

        if (proximity < minProximity)
        {
            nearestCelli = celli;
            minProximity = proximity;
        }
    }

    return nearestCelli;
}


// walking from seed
Foam::label Foam::meshSearch::findNearestCellWalk
(
    const point& location,
    const label seedCellI
) const
{
    if (seedCellI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestCellWalk(const point&, const label)"
        )   << "illegal seedCell:" << seedCellI << exit(FatalError);
    }

    const vectorField& centres = mesh_.cellCentres();
    const labelListList& cc = mesh_.cellCells();


    // Walk in direction of face that decreases distance

    label curCell = seedCellI;
    scalar distanceSqr = magSqr(centres[curCell] - location);

    bool closer;

    do
    {
        closer = false;

        // set the current list of neighbouring cells
        const labelList& neighbours = cc[curCell];

        forAll (neighbours, nI)
        {
            scalar curDistSqr = magSqr(centres[neighbours[nI]] - location);

            // search through all the neighbours.
            // If the cell is closer, reset current cell and distance
            if (curDistSqr < distanceSqr)
            {
                distanceSqr = curDistSqr;
                curCell = neighbours[nI];
                closer = true;    // a closer neighbour has been found
            }
        }
    } while (closer);

    return curCell;
}


Foam::label Foam::meshSearch::findCellLinear(const point& location) const
{
    bool cellFound = false;
    label n = 0;

    label cellI = -1;

    while ((!cellFound) && (n < mesh_.nCells()))
    {
        if (pointInCell(location, n))
        {
            cellFound = true;
            cellI = n;
        }
        else
        {
            n++;
        }
    }
    if (cellFound)
    {
        return cellI;
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::meshSearch::findNearestBoundaryFaceWalk
(
    const point& location,
    const label seedFaceI
) const
{
    if (seedFaceI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestBoundaryFaceWalk"
            "(const point&, const label)"
        )   << "illegal seedFace:" << seedFaceI << exit(FatalError);
    }

    // Start off from seedFaceI

    label curFaceI = seedFaceI;

    const face& f =  mesh_.faces()[curFaceI];

    scalar minDist =
        f.nearestPoint
        (
            location,
            mesh_.points()
        ).distance();

    bool closer;

    do
    {
        closer = false;

        // Search through all neighbouring boundary faces by going
        // across edges

        label lastFaceI = curFaceI;

        const labelList& myEdges = mesh_.faceEdges()[curFaceI];

        forAll (myEdges, myEdgeI)
        {
            const labelList& neighbours =
                mesh_.edgeFaces()[myEdges[myEdgeI]];

            // Check any face which uses edge, is boundary face and
            // is not curFaceI itself.

            forAll(neighbours, nI)
            {
                label faceI = neighbours[nI];

                if
                (
                    (faceI >= mesh_.nInternalFaces())
                 && (faceI != lastFaceI)
                )
                {
                    const face& f =  mesh_.faces()[faceI];

                    pointHit curHit =
                        f.nearestPoint
                        (
                            location,
                            mesh_.points()
                        );

                    // If the face is closer, reset current face and distance
                    if (curHit.distance() < minDist)
                    {
                        minDist = curHit.distance();
                        curFaceI = faceI;
                        closer = true;  // a closer neighbour has been found
                    }
                }
            }
        }
    } while (closer);

    return curFaceI;
}


Foam::vector Foam::meshSearch::offset
(
    const point& bPoint,
    const label bFaceI,
    const vector& dir
) const
{
    // Get the neighbouring cell
    label ownerCellI = mesh_.faceOwner()[bFaceI];

    const point& c = mesh_.cellCentres()[ownerCellI];

    // Typical dimension: distance from point on face to cell centre
    scalar typDim = mag(c - bPoint);

    return tol_*typDim*dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshSearch::meshSearch(const polyMesh& mesh, const bool faceDecomp)
:
    mesh_(mesh),
    faceDecomp_(faceDecomp),
    cloud_(mesh_, IDLList<passiveParticle>()),
    boundaryTreePtr_(NULL),
    cellTreePtr_(NULL),
    cellCentreTreePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSearch::~meshSearch()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::octree<Foam::octreeDataFace>& Foam::meshSearch::boundaryTree() const
{
    if (!boundaryTreePtr_)
    {
        //
        // Construct tree
        //

        // all boundary faces (not just walls)
        octreeDataFace shapes(mesh_);

        treeBoundBox overallBb(mesh_.points());

        boundaryTreePtr_ = new octree<octreeDataFace>
        (
            overallBb,  // overall search domain
            shapes,     // all information needed to do checks on cells
            1,          // min levels
            20.0,       // maximum ratio of cubes v.s. cells
            10.0
        );
    }

    return *boundaryTreePtr_;
}


const Foam::octree<Foam::octreeDataCell>& Foam::meshSearch::cellTree() const
{
    if (!cellTreePtr_)
    {
        //
        // Construct tree
        //

        octreeDataCell shapes(mesh_);

        treeBoundBox overallBb(mesh_.points());

        cellTreePtr_ = new octree<octreeDataCell>
        (
            overallBb,  // overall search domain
            shapes,     // all information needed to do checks on cells
            1,          // min levels
            20.0,       // maximum ratio of cubes v.s. cells
            2.0
        );
    }

    return *cellTreePtr_;
    
}


const Foam::octree<Foam::octreeDataPoint>& Foam::meshSearch::cellCentreTree()
    const
{
    if (!cellCentreTreePtr_)
    {
        //
        // Construct tree
        //

        // shapes holds reference to cellCentres!
        octreeDataPoint shapes(mesh_.cellCentres());

        treeBoundBox overallBb(mesh_.cellCentres());

        cellCentreTreePtr_ = new octree<octreeDataPoint>
        (
            overallBb,  // overall search domain
            shapes,     // all information needed to do checks on cells
            1,          // min levels
            20.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        );
    }

    return *cellCentreTreePtr_;
}


// Is the point in the cell
// Works by checking if there is a face inbetween the point and the cell
// centre.
// Check for internal uses proper face decomposition or just average normal.
bool Foam::meshSearch::pointInCell(const point& p, label cellI) const
{
    if (faceDecomp_)
    {
        const point& ctr = mesh_.cellCentres()[cellI];

        vector dir(p - ctr);
        scalar magDir = mag(dir);

        // Check if any faces are hit by ray from cell centre to p.
        // If none -> p is in cell.
        const labelList& cFaces = mesh_.cells()[cellI];

        // Make sure half_ray does not pick up any faces on the wrong
        // side of the ray.
        scalar oldTol = intersection::setPlanarTol(0.0);

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            const face& f = mesh_.faces()[faceI];

            forAll(f, fp)
            {
                pointHit inter = 
                    f.ray
                    (
                        ctr,
                        dir,
                        mesh_.points(),
                        intersection::HALF_RAY,
                        intersection::VECTOR
                    );

                if (inter.hit())
                {
                    scalar dist = inter.distance();

                    if (dist < magDir)
                    {
                        // Valid hit. Hit face so point is not in cell.
                        intersection::setPlanarTol(oldTol);

                        return false;
                    }
                }
            }
        }

        intersection::setPlanarTol(oldTol);

        // No face in between point and cell centre so point is inside.
        return true;
    }
    else
    {
        const labelList& f = mesh_.cells()[cellI];
        const labelList& owner = mesh_.faceOwner();
        const vectorField& cf = mesh_.faceCentres();
        const vectorField& Sf = mesh_.faceAreas();

        forAll(f, facei)
        {
            label nFace = f[facei];
            vector proj = p - cf[nFace];
            vector normal = Sf[nFace];
            if (owner[nFace] == cellI)
            {
                if ((normal & proj) > 0)
                {
                    return false;
                }
            }
            else
            {
                if ((normal & proj) < 0)
                {
                    return false;
                }
            }
        }

        return true;
    }
}


Foam::label Foam::meshSearch::findNearestCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    if (seedCellI == -1)
    {
        if (useTreeSearch)
        {
            return findNearestCellTree(location);
        }
        else
        {
            return findNearestCellLinear(location);
        }
    }
    else
    {
        return findNearestCellWalk(location, seedCellI);
    }
}


Foam::label Foam::meshSearch::findCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    // Find the nearest cell centre to this location
    label nearCellI = findNearestCell(location, seedCellI, useTreeSearch);

    if (debug)
    {
        Pout<< "findCell : nearCellI:" << nearCellI
            << " ctr:" << mesh_.cellCentres()[nearCellI]
            << endl;
    }

    if (pointInCell(location, nearCellI))
    {
        return nearCellI;
    }
    else
    {
        if (useTreeSearch)
        {
            // Start searching from cell centre of nearCell
            point curPoint = mesh_.cellCentres()[nearCellI];

            vector edgeVec = location - curPoint;
            edgeVec /= mag(edgeVec);

            do
            {
                // Walk neighbours (using tracking) to get closer
                passiveParticle tracker(cloud_, curPoint, nearCellI);

                if (debug)
                {
                    Pout<< "findCell : tracked from curPoint:" << curPoint
                        << " nearCellI:" << nearCellI;
                }

                tracker.track(location);

                if (debug)
                {
                    Pout<< " to " << tracker.position()
                        << " need:" << location
                        << " onB:" << tracker.onBoundary()
                        << " cell:" << tracker.cell()
                        << " face:" << tracker.face() << endl;
                }

                if (!tracker.onBoundary())
                {
                    // stopped not on boundary -> reached location
                    return tracker.cell();
                }

                // stopped on boundary face. Compare positions
                scalar typDim = sqrt(mag(mesh_.faceAreas()[tracker.face()]));

                if ((mag(tracker.position() - location)/typDim) < SMALL)
                {
                    return tracker.cell();
                }

                // tracking stopped at boundary. Find next boundary
                // intersection from current point (shifted out a little bit)
                // in the direction of the destination

                curPoint =
                    tracker.position()
                  + offset(tracker.position(), tracker.face(), edgeVec);

                if (debug)
                {
                    Pout<< "Searching for next boundary from curPoint:"
                        << curPoint
                        << " to location:" << location  << endl;
                }
                pointIndexHit curHit = intersection(curPoint, location);
                if (debug)
                {
                    Pout<< "Returned from line search with ";
                    curHit.write(Pout);
                    Pout<< endl;
                }

                if (!curHit.hit())
                {
                    return -1;
                }
                else
                {
                    nearCellI = mesh_.faceOwner()[curHit.index()];
                    curPoint =
                        curHit.hitPoint() 
                      + offset(curHit.hitPoint(), curHit.index(), edgeVec);
                }
            }
            while(true);
        }
        else
        {
             return findCellLinear(location);
        }
    }
    return -1;
}


Foam::label Foam::meshSearch::findNearestBoundaryFace
(
    const point& location,
    const label seedFaceI,
    const bool useTreeSearch
) const
{
    if (seedFaceI == -1)
    {
        if (useTreeSearch)
        {
            treeBoundBox tightest(treeBoundBox::greatBox);

            scalar tightestDist = GREAT;

            label index =
                boundaryTree().findNearest
                (
                    location,
                    tightest,
                    tightestDist
                );

            return boundaryTree().shapes().meshFaces()[index];
        }
        else
        {
            scalar minDist = GREAT;

            label minFaceI = -1;

            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                const face& f =  mesh_.faces()[faceI];

                pointHit curHit =
                    f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minFaceI = faceI;
                }
            }
            return minFaceI;
        }
    }
    else
    {
        return findNearestBoundaryFaceWalk(location, seedFaceI);
    }
}


Foam::pointIndexHit Foam::meshSearch::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    pointIndexHit curHit = boundaryTree().findLine(pStart, pEnd);

    if (curHit.hit())
    {
        // Change index into octreeData into face label
        curHit.setIndex(boundaryTree().shapes().meshFaces()[curHit.index()]);
    }
    return curHit;
}


Foam::List<Foam::pointIndexHit> Foam::meshSearch::intersections
(
    const point& pStart,
    const point& pEnd
) const
{
    DynamicList<pointIndexHit> hits;

    vector edgeVec = pEnd - pStart;
    edgeVec /= mag(edgeVec);

    point pt = pStart;

    pointIndexHit bHit;
    do
    {
        bHit = intersection(pt, pEnd);

        if (bHit.hit())
        {
            // Verify if this hit point has not been visited before
            // Otherwise, we are entering an infinite loop.
            // HJ, bug fix  10/Nov/2008
            bool alreadyVisitedPoint = false;
            forAll(hits, hI)
            {
                if(bHit.hitPoint() == hits[hI].hitPoint())
                {
                    alreadyVisitedPoint = true;
                    break;
                }
            }
            //- Start of infinite loop?
            if(alreadyVisitedPoint) 
                break;

            hits.append(bHit);

            const vector& area = mesh_.faceAreas()[bHit.index()];

            scalar typDim = Foam::sqrt(mag(area));

            if ((mag(bHit.hitPoint() - pEnd)/typDim) < SMALL)
            {
                break;
            }

            // Restart from hitPoint shifted a little bit in the direction
            // of the destination

            pt =
                bHit.hitPoint()
              + offset(bHit.hitPoint(), bHit.index(), edgeVec);
        }

    } while(bHit.hit());


    hits.shrink();

    // Copy into straight list
    List<pointIndexHit> allHits(hits.size());

    forAll(hits, hitI)
    {
        allHits[hitI] = hits[hitI];
    }
    return allHits;
}


bool Foam::meshSearch::isInside(const point& p) const
{
    return boundaryTree().getSampleType(p) == octree<octreeDataFace>::INSIDE;
}


// Delete all storage
void Foam::meshSearch::clearOut()
{
    deleteDemandDrivenData(boundaryTreePtr_);
    deleteDemandDrivenData(cellTreePtr_);
    deleteDemandDrivenData(cellCentreTreePtr_);
}


void Foam::meshSearch::correct()
{
    clearOut();
}


// ************************************************************************* //
